"""E115: bounded contact-region, polar-exposure and feasible-pose diagnostics."""
import argparse
from pathlib import Path
import time
import numpy as np
import audit_native_state_poses as native
from check_refined_shape_poses import assess
from diagnose_contact_features import controlled_candidate
from aidd_agent import contact_region_review as review
base = native.base


def export_primary(out, candidate, chemistry, state, query, matrices):
    from aidd_agent.consensus_admission import read_structure
    from aidd_agent.predocking_qc import _chemical_sdf
    from aidd_agent.anchor_extraction import protein_hbond_roles
    record = state["record"]
    structure = read_structure(Path(record["native_mmcif"]), query["target"])
    atoms = [a for a in structure["chains"][record["chain"]] if not a["label_alt_id"]
             and float(a["occupancy"]) >= .9 and not a.get("pdbx_PDB_ins_code")]
    receptor = native.move([a["xyz"] for a in atoms], native.rigid(record["transform"]))
    base.require(np.allclose(receptor, state["tree"].data), "Receptor export frame mismatch")
    lines = []
    for i, (atom, xyz) in enumerate(zip(atoms, receptor), 1):
        lines.append(f"ATOM  {i:5d} {atom['auth_atom_id']:>4s} {atom['auth_comp_id']:>3s} "
                     f"A{int(atom['auth_seq_id']):4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                     f"  1.00  0.00          {atom['type_symbol']:>2s}")
    (out / "primary_receptor.pdb").write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    base.save(out / "primary_receptor_mapping.json", [
        dict(export_serial=i, native_chain=a["auth_asym_id"], native_residue=a["auth_seq_id"],
             canonical_residue=a["canonical_residue"], atom=a["auth_atom_id"])
        for i, a in enumerate(atoms, 1)])
    details = {}
    for name, matrix in matrices.items():
        points = native.move(candidate.shape_points, matrix)
        (out / f"primary_{name}.sdf").write_text(
            _chemical_sdf(points, chemistry, "6Q9L:HTZ:A:201"), encoding="utf-8")
        attribution = review.color_attribution(
            state["aligned"]["feature_points"], state["aligned"]["feature_types"],
            native.move(candidate.feature_points, matrix), candidate.feature_types)
        for row in attribution["candidate_features"]:
            row["member_atoms"] = list(map(int, chemistry.feature_members[row["index"]]))
        exposure = review.polar_exposure(points, chemistry, candidate.feature_types,
                                         receptor, [a["type_symbol"] for a in atoms])
        polar_neighbors = []
        for i, kind in enumerate(candidate.feature_types):
            if int(kind) not in (1, 2):
                continue
            center = native.move(candidate.feature_points[i:i+1], matrix)[0]
            for j in np.flatnonzero(np.linalg.norm(receptor-center, axis=1) <= 3.5):
                atom = atoms[j]
                donor, acceptor = protein_hbond_roles(atom["auth_comp_id"], atom["auth_atom_id"])
                if (int(kind) == 1 and acceptor) or (int(kind) == 2 and donor):
                    polar_neighbors.append(dict(candidate_feature=i,
                        protein_residue=atom["canonical_residue"], atom=atom["auth_atom_id"],
                        distance=float(np.linalg.norm(receptor[j]-center))))
        details[name] = dict(color=attribution, polar_exposure=exposure,
                             complementary_polar_neighbors_distance_only=polar_neighbors)
    base.save(out / "primary_review.json", dict(molecule=candidate.molecule_id,
        template="6Q9L:HTZ:A:201", poses=details,
        limitations=["SDF preserves companion connectivity/charges but does not reconstruct stereo tags.",
                     "Feature assignments are geometric hypotheses, not proven hydrogen bonds.",
                     "SASA excludes waters, hydrogens and other chains; no permeability prediction."]))
    (out / "view_primary.pml").write_text(
        "load primary_receptor.pdb, receptor\nload primary_selected.sdf, ligand\n"
        "hide everything\nshow surface, receptor\ncolor lightblue, receptor\n"
        "set transparency, 0.45, receptor\nshow sticks, ligand\n"
        "color indigo, ligand\ncolor red, ligand and elem O\n"
        "color blue, ligand and elem N\nbg_color white\nzoom ligand, 8\n", encoding="utf-8")


def pilot(args, selected, states, query, batch, out, controls, started):
    from aidd_agent.gaussian_batch import ArtifactCatalogReader
    from aidd_agent.interaction_fast import InteractionFeatureReader
    from aidd_agent.chemical_companion import ChemicalCompanionReader
    bundle = base.read(args.proposals)
    receipt = args.receipt.resolve()
    base.verify(receipt / "report.json", bundle["source_report_sha256"])
    base.verify(receipt / "protocol.json", bundle["source_protocol_sha256"])
    report = base.read(receipt / "report.json")
    protocol = base.read(receipt / "protocol.json")
    base.require(report["status"] == "complete", "Incomplete E113 receipt")
    for name, digest in report["output_hashes"].items():
        base.verify(receipt / name, digest)
    base.require(protocol["query_sha256"] == base.read(out / "protocol.json")["query_sha256"],
                 "Query lineage mismatch")
    base.require(protocol["previous_report_sha256"] == base.sha(args.previous / "report.json"),
                 "Sample lineage mismatch")
    evidence = {(r["kind"], r["id"]): r for r in base.read(receipt / "pose_diagnostics.json")}
    sample = {str(r["gid"]): r for r in selected}
    groups = {key: review.reference_groups(state) for key, state in states.items()}
    base.save(out / "reference_groups.json", groups)
    reader = ArtifactCatalogReader(batch / "artifacts/catalog.json")
    features = InteractionFeatureReader(batch / "artifacts/catalog.json", batch / "chemical/catalog.json")
    chemistry_reader = ChemicalCompanionReader(batch / "chemical/catalog.json")
    rows, identities = [], []
    deadline = time.monotonic() + args.seconds
    try:
        # Endpoint replay is mandatory even when the optimization budget is exhausted.
        for p in bundle["proposals"]:
            source = evidence[(p["kind"], p["id"])]["best"][p["objective"]]
            base.require(source["template"] == p["template"] and np.allclose(
                source["transform"], p["original_transform"], atol=1e-12), "Proposal source mismatch")
            state = states[p["template"]]
            if p["kind"] == "library":
                candidate = reader.get(int(p["id"])); cf = features.get(int(p["id"]))
                item = sample[p["id"]]
                base.require(candidate.conformer_id == cf.conformer_id == item["cid"] and
                             candidate.molecule_id == cf.molecule_id == item["mid"] == p["mid"],
                             "Candidate identity mismatch")
            else:
                base.require(p["kind"] == "leave_self_out", "Unexpected control kind")
                candidate, cf = controlled_candidate(states[p["id"]], p["id"])
            endpoints = [native.rigid(p[k]) for k in ("original_transform", "transform")]
            def evaluate(matrix):
                score = assess(candidate, cf, state, matrix)
                score["grouped"] = review.grouped_score(state, cf, matrix, groups[p["template"]])
                return score
            result = review.refine_feasible(candidate.shape_points, endpoints, evaluate, deadline)
            deadline += result["endpoint_replay_seconds"]
            for endpoint, field in zip(result["archive"][:2], ("original_shape", "optimized_shape")):
                base.require(abs(endpoint["score"]["shape"] - p[field]) < 1e-6, "Endpoint replay mismatch")
            row = dict(kind=p["kind"], id=p["id"], mid=p["mid"], template=p["template"], **result)
            rows.append(row)
            if p["kind"] == "library" and str(p["id"]) == "4382601":
                chemistry = chemistry_reader.get(int(p["id"]))
                base.require(chemistry.conformer_id == candidate.conformer_id and
                             chemistry.molecule_id == candidate.molecule_id, "Chemistry identity mismatch")
                chosen = result["selected"]
                base.require(chosen is not None, "Primary candidate lost feasible endpoint")
                export_started = time.monotonic()
                export_primary(out, candidate, chemistry, state, query,
                    dict(original=endpoints[0], shape_only=endpoints[1],
                         selected=native.rigid(chosen["transform"])))
                deadline += time.monotonic() - export_started
            print(f"E115 reviewed {len(rows)}/{len(bundle['proposals'])}", flush=True)
        for key, state in states.items():
            if not state["record"]["supported_for_requested_contact_rule"]:
                continue
            data = state["aligned"]
            from types import SimpleNamespace
            cf = SimpleNamespace(feature_points=data["feature_points"], feature_types=data["feature_types"],
                feature_directions=data["feature_directions"], feature_kinds=data["feature_direction_kinds"])
            score = review.grouped_score(state, cf, np.eye(4), groups[key])
            base.require(abs(score["mean"] - 1) < 1e-6, "Grouped identity control failed")
            identities.append(dict(template=key, **score))
    finally:
        features.close()
        if hasattr(reader, "close"):
            reader.close()
    base.save(out / "reviewed_poses.json", rows)
    base.save(out / "identity_controls.json", identities)
    base.write_csv(out / "pose_summary.csv", [
        dict(kind=r["kind"], id=r["id"], mid=r["mid"], template=r["template"],
             original_feasible=review.feasible(r["archive"][0]["score"]),
             original_shape=r["archive"][0]["score"]["shape"],
             original_color=r["archive"][0]["score"]["color"],
             original_grouped=r["archive"][0]["score"]["grouped"]["mean"],
             selected_shape=r["selected"]["score"]["shape"] if r["selected"] else None,
             selected_color=r["selected"]["score"]["color"] if r["selected"] else None,
             selected_grouped=r["selected"]["score"]["grouped"]["mean"] if r["selected"] else None,
             original_policy_preserved=r["original_policy_preserved"],
             trials=r["trial_evaluations"], budget_exhausted=r["budget_exhausted"])
        for r in rows])
    settings = base.read(out / "protocol.json")
    settings.update(mode="E115_no_library_search", runner_sha256=base.sha(__file__),
        module_sha256=base.sha(review.__file__), proposals_sha256=base.sha(args.proposals),
        e113_report_sha256=bundle["source_report_sha256"], optimization_seconds=args.seconds,
        protocol_sha256=base.sha(Path(__file__).resolve().parents[1] /
                                 "experiments/E115-contact-region-review.md"))
    settings["diagnostic_dependencies_sha256"] = {
        name: base.sha(Path(__file__).with_name(name)) for name in (
            "audit_native_state_poses.py", "audit_strict_cached_poses.py",
            "check_refined_shape_poses.py", "diagnose_contact_features.py")}
    base.save(out / "protocol.json", settings)
    result = dict(status="complete", optimization_budget_exhausted=any(r["budget_exhausted"] for r in rows),
        total_seconds=time.monotonic()-started, poses=len(rows), identity_controls=len(identities),
        original_feasible=sum(review.feasible(r["archive"][0]["score"]) for r in rows),
        selected_feasible=sum(r["selected"] is not None for r in rows),
        retained_original_policy_positive=sum(r["original_policy_preserved"] for r in rows),
        limitations=["Outcome-selected panel, not a pass-rate estimate.",
                     "No production gate changes; grouped score has no calibrated acceptance threshold.",
                     "Finite local rigid proposals, no conformer or global search.",
                     "Complete means all endpoint reviews finished; inspect optimization budget flag."])
    result["output_hashes"] = {p.name: base.sha(p) for p in out.iterdir()
                               if p.is_file() and p.name != "report.json"}
    base.save(out / "report.json", result)
    print(base.json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("previous", "prepared", "reference", "output", "receipt"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--proposals", type=Path, default=Path(__file__).resolve().parents[1] /
                        "to_human/E114_SHAPE_PROPOSALS.json")
    parser.add_argument("--seconds", type=float, default=180)
    args = parser.parse_args()
    if not np.isfinite(args.seconds) or args.seconds <= 0:
        parser.error("--seconds must be finite and positive")
    args.fresh = True; args.sample = 10000; args.workers = 1
    native.fresh_pilot = pilot
    native.run(args)


if __name__ == "__main__":
    main()
