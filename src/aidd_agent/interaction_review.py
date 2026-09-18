"""Molecule-level diagnostic ranks and immutable native-frame pose exports."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil

import numpy as np

from .gaussian_batch import ArtifactCatalogReader
from .chemical_companion import ChemicalCompanionReader
from .gaussian_overlay import apply_transform
from .interaction_matching import _average_ranks
from .predocking_qc import _chemical_sdf

OBJECTIVE = "atomcentered_anchored_joint"


def molecule_comparison(ids, molecules, gaussian, interaction, *, top=100, outside=500, per_group=10):
    ids, molecules, gaussian, interaction = map(np.asarray, (ids, molecules, gaussian, interaction))
    if not (ids.ndim == 1 and ids.shape == molecules.shape == gaussian.shape == interaction.shape
            and len(np.unique(ids)) == len(ids) and len(ids) > 0
            and np.isfinite(gaussian).all() and np.isfinite(interaction).all()):
        raise ValueError("invalid molecule comparison arrays")
    representatives = []
    for scores in (gaussian, interaction):
        best = {}
        for index in np.lexsort((ids, -scores)):
            best.setdefault(str(molecules[index]), int(index))
        representatives.append(best)
    gbest, ibest = representatives
    rows = []
    for mid in sorted(gbest):
        g, i = gbest[mid], ibest[mid]
        rows.append(dict(molecule_id=mid, gaussian_index=g, interaction_index=i,
                         gaussian_global_id=int(ids[g]), interaction_global_id=int(ids[i]),
                         gaussian_score=float(gaussian[g]), interaction_score=float(interaction[i]),
                         interaction_at_gaussian_pose=float(interaction[g]), gaussian_at_interaction_pose=float(gaussian[i])))
    for metric in ("gaussian", "interaction"):
        order = sorted(range(len(rows)), key=lambda j: (-rows[j][metric + "_score"], rows[j][metric + "_global_id"]))
        for rank, index in enumerate(order, 1): rows[index][metric + "_rank"] = rank
    g = np.asarray([r["gaussian_score"] for r in rows]); i = np.asarray([r["interaction_score"] for r in rows])
    rho = float(np.corrcoef(_average_ranks(g), _average_ranks(i))[0, 1]) if np.ptp(g) and np.ptp(i) else None
    overlap = {}
    for k in (100, 500, 1000):
        actual = min(k, len(rows))
        shared = sum(r["gaussian_rank"] <= actual and r["interaction_rank"] <= actual for r in rows)
        overlap[str(k)] = dict(denominator=actual, count=shared, fraction=shared / actual)
    groups = {
        "gaussian_high_interaction_low": sorted([r for r in rows if r["gaussian_rank"] <= top and r["interaction_rank"] > outside],
                                                 key=lambda r: (r["gaussian_rank"], -r["interaction_rank"])),
        "interaction_high_gaussian_low": sorted([r for r in rows if r["interaction_rank"] <= top and r["gaussian_rank"] > outside],
                                                 key=lambda r: (r["interaction_rank"], -r["gaussian_rank"])),
        "both_high": sorted([r for r in rows if max(r["gaussian_rank"], r["interaction_rank"]) <= top],
                            key=lambda r: (max(r["gaussian_rank"], r["interaction_rank"]), r["molecule_id"]))}
    selected = {name: values[:per_group] for name, values in groups.items()}
    return rows, dict(molecules=len(rows), spearman=rho, top_k_overlap=overlap,
                     different_representatives=sum(r["gaussian_global_id"] != r["interaction_global_id"] for r in rows),
                     eligible_groups={name: len(values) for name, values in groups.items()},
                     selected_groups={name: len(values) for name, values in selected.items()},
                     unit="source-grouped molecule; diagnostic ranks only"), selected


def write_review(artifact, chemical, query_dir, rigid, sidecar, output, query_id):
    output.mkdir(parents=True, exist_ok=True)
    with np.load(rigid, allow_pickle=False) as archive:
        r = {k: archive[k] for k in archive.files}
    with np.load(sidecar, allow_pickle=False) as archive:
        s = {k: archive[k] for k in archive.files}
    for key in ("global_ids", "molecule_ids", "conformer_ids"):
        if not np.array_equal(r[key], s[key]): raise ValueError("review input identity mismatch")
    rows, summary, selected = molecule_comparison(r["global_ids"], r["molecule_ids"], r[f"{OBJECTIVE}__objective"],
                                                 s[f"{OBJECTIVE}__interaction_match_score"])
    with (output / "molecule-ranks.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    pdb_id = query_id.split(":")[0]
    shutil.copyfile(query_dir / f"{pdb_id}.cif", output / "receptor.cif")
    with np.load(query_dir / "gaussian-query.npz", allow_pickle=False) as q:
        anchor_indices = q["anchor_feature_indices"]; anchor_points = q["feature_points"][anchor_indices]
    reader, chem_reader = ArtifactCatalogReader(artifact), ChemicalCompanionReader(chemical)
    exports = []
    for group, molecules in selected.items():
        for number, row in enumerate(molecules, 1):
            # Export each distinct winning conformer using the same anchored objective pose.
            for index in sorted({row["gaussian_index"], row["interaction_index"]}):
                gid = int(r["global_ids"][index]); candidate = reader.get(gid); chemistry = chem_reader.get(gid)
                if (candidate.molecule_id != str(r["molecule_ids"][index])
                        or candidate.conformer_id != str(r["conformer_ids"][index])
                        or chemistry.molecule_id != candidate.molecule_id or chemistry.conformer_id != candidate.conformer_id):
                    raise ValueError("pose export identity mismatch")
                transform = r[f"{OBJECTIVE}__transform"][index]
                points = apply_transform(candidate.shape_points, transform)
                filename = f"{group}-{number:02d}-gid{gid}.sdf"
                (output / filename).write_text(_chemical_sdf(points, chemistry, query_id), encoding="utf-8")
                assignments = s[f"{OBJECTIVE}__anchor_assignments"][index]
                contributions = s[f"{OBJECTIVE}__anchor_scores"][index]
                moved_features = apply_transform(candidate.feature_points, transform)
                pairs = [dict(anchor_feature_index=int(anchor_indices[a]), candidate_feature_index=int(c),
                              anchor_score=float(contributions[a]), query_xyz=anchor_points[a].tolist(),
                              candidate_xyz=moved_features[c].tolist() if c >= 0 else None)
                         for a, c in enumerate(assignments)]
                exports.append(dict(group=group, file=filename, global_id=gid, molecule_id=candidate.molecule_id,
                                    conformer_id=candidate.conformer_id, transform=transform.tolist(),
                                    gaussian_score=float(r[f"{OBJECTIVE}__objective"][index]),
                                    interaction_score=float(s[f"{OBJECTIVE}__interaction_match_score"][index]),
                                    anchors=pairs, molecule_comparison=row))
    (output / "poses.json").write_text(json.dumps(exports, indent=2), encoding="utf-8")
    (output / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (output / "review-context.json").write_text(json.dumps(dict(query_id=query_id)), encoding="utf-8")
    (output / "review.py").write_text('''from pathlib import Path
import json
from pymol import cmd
root = Path(__file__).resolve().parent
cmd.load(str(root / "receptor.cif"), "crystal")
_, ccd, chain, residue = json.loads((root / "review-context.json").read_text())["query_id"].split(":")
cmd.select("query_ligand", "crystal and resn %s and chain %s and resi %s" % (ccd, chain, residue))
cmd.hide("everything", "all")
cmd.show("cartoon", "crystal and polymer")
cmd.show("sticks", "query_ligand")
for n, row in enumerate(json.loads((root / "poses.json").read_text())):
    name = "review_%03d_gid%d" % (n, row["global_id"])
    cmd.load(str(root / row["file"]), name)
    cmd.show("sticks", name)
    cmd.group(row["group"], name)
    for a, pair in enumerate(row["anchors"]):
        if pair["candidate_xyz"] is None:
            continue
        q, c, d = name + "_q%d" % a, name + "_c%d" % a, name + "_match%d" % a
        cmd.pseudoatom(q, pos=pair["query_xyz"])
        cmd.pseudoatom(c, pos=pair["candidate_xyz"])
        cmd.distance(d, q, c)
        cmd.group(name + "_anchors", "%s %s %s" % (q, c, d))
        cmd.group(row["group"], name + "_anchors")
        cmd.disable(name + "_anchors")
    cmd.disable(name)
cmd.zoom("query_ligand")
print("Enable one review object at a time. Match lines are scoring assignments, not validated hydrogen bonds.")
''', encoding="utf-8")
    (output / "README.md").write_text(
        "# Diagnostic pose review\n\nRun `run /absolute/path/to/review.py` in PyMOL. Enable one pose at a time.\n"
        "Both winning conformers are exported when different. Empty groups are not backfilled.\n"
        "Coordinates are in the original crystal frame; no docking or minimization was performed.\n"
        "SDF contains heavy atoms/bonds/formal charges; full source stereochemical annotation is not reconstructed.\n"
        "Assignments are not validated hydrogen bonds. Inspect clashes, type/direction and per-anchor scores.\n",
        encoding="utf-8")
    return dict(**summary, exported_poses=len(exports), output=str(output))
