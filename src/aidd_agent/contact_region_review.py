"""Exploratory reference grouping and feasibility-preserving rigid refinement."""
import time
import numpy as np
from scipy.spatial.transform import Rotation
from .chemical_geometry import DIRECTION_AXIAL, DIRECTION_NONE
from .interaction_matching import interaction_match, _maximum_weight_assignment


def reference_groups(state):
    data = state["aligned"]
    labels = {}
    for anchor in state["anchors"]:
        labels.setdefault(int(anchor["feature_index"]), set()).add(
            (int(anchor["target_residue"]), anchor["protein_atom"], anchor["feature_class"]))
    groups = []
    for index in sorted(labels):
        def compatible(other):
            if labels[index] != labels[other]:
                return False
            if any(data[k][index] != data[k][other] for k in
                   ("feature_types", "feature_direction_kinds")):
                return False
            if np.linalg.norm(data["feature_points"][index] - data["feature_points"][other]) > 2:
                return False
            kind = data["feature_direction_kinds"][index]
            cosine = np.dot(data["feature_directions"][index], data["feature_directions"][other])
            return kind == DIRECTION_NONE or (
                abs(cosine) if kind == DIRECTION_AXIAL else cosine) >= np.cos(np.pi / 6)
        destination = next((g for g in groups if all(compatible(j) for j in g)), None)
        if destination is None:
            groups.append([index])
        else:
            destination.append(index)
    return [dict(reference_features=g, labels=[list(x) for x in sorted(labels[g[0]])])
            for g in groups]


def grouped_score(state, features, matrix, groups):
    if not groups:
        return dict(mean=0.0, groups=[])
    data = state["aligned"]
    indices = sorted({i for g in groups for i in g["reference_features"]})
    pairs = interaction_match(
        data["feature_points"][indices], data["feature_types"][indices],
        data["feature_directions"][indices], data["feature_direction_kinds"][indices],
        np.ones(len(indices)),
        features.feature_points @ matrix[:3, :3].T + matrix[:3, 3],
        features.feature_types, features.feature_directions @ matrix[:3, :3].T,
        features.feature_kinds)["pair_scores"]
    lookup = {j: i for i, j in enumerate(indices)}
    values = np.asarray([pairs[[lookup[i] for i in g["reference_features"]]].max(axis=0)
                         for g in groups])
    assignments = _maximum_weight_assignment(values, np.ones(len(groups)))
    rows = []
    for n, (group, candidate) in enumerate(zip(groups, assignments)):
        member = None if candidate < 0 else max(
            group["reference_features"], key=lambda i: pairs[lookup[i], candidate])
        rows.append(dict(**group, candidate_feature=int(candidate), matched_reference_feature=member,
                         score=0.0 if candidate < 0 else float(values[n, candidate])))
    return dict(mean=float(np.mean([r["score"] for r in rows])), groups=rows)


def color_attribution(reference, reference_types, candidate, candidate_types):
    def kernel(a, b):
        return np.exp(-np.sum((a[:, None] - b[None, :]) ** 2, axis=2) / 2)
    cross = kernel(reference, candidate) * (reference_types[:, None] == candidate_types)
    left = kernel(reference, reference) * (reference_types[:, None] == reference_types)
    right = kernel(candidate, candidate) * (candidate_types[:, None] == candidate_types)
    denominator = float(left.sum() + right.sum() - cross.sum())
    return dict(tanimoto=float(cross.sum() / denominator) if denominator else 0.0,
                reference_self=float(left.sum()), candidate_self=float(right.sum()),
                cross=float(cross.sum()),
                types=[dict(type=int(kind),
                    reference_self=float(left[reference_types == kind].sum()),
                    candidate_self=float(right[candidate_types == kind].sum()),
                    cross=float(cross[reference_types == kind].sum()))
                    for kind in sorted(set(reference_types) | set(candidate_types))],
                candidate_features=[dict(index=i, type=int(candidate_types[i]),
                    cross=float(cross[:, i].sum()), self=float(right[i].sum()),
                    best_reference_overlap=float(cross[:, i].max()) if len(reference) else 0.0)
                    for i in range(len(candidate))],
                reference_features=[dict(index=i, type=int(reference_types[i]),
                    cross=float(cross[i].sum()), self=float(left[i].sum()))
                    for i in range(len(reference))])


def feasible(score):
    return all(score["flags"][k] for k in ("contacts", "steric", "envelope", "regions"))


def refine_feasible(points, endpoints, evaluate, deadline, max_trials=96):
    """Keep endpoints unconditionally; never replace a feasible pose by an invalid one."""
    replay_started = time.monotonic()
    archive = [dict(transform=np.asarray(m).reshape(4, 4).ravel().tolist(), score=evaluate(
        np.asarray(m).reshape(4, 4)), origin="endpoint") for m in endpoints]
    replay_seconds = time.monotonic() - replay_started
    deadline += replay_seconds
    feasible_rows = [r for r in archive if feasible(r["score"])]
    incumbent = max(feasible_rows, key=lambda r: r["score"]["shape"]) if feasible_rows else None
    trials = 0
    for translation, rotation in ((.20, .08), (.08, .03)):
        improved = True
        while incumbent is not None and improved and trials < max_trials and time.monotonic() < deadline:
            improved = False
            matrix = np.array(incumbent["transform"]).reshape(4, 4)
            center = (points @ matrix[:3, :3].T + matrix[:3, 3]).mean(axis=0)
            for axis in range(6):
                for sign in (-1, 1):
                    if trials >= max_trials or time.monotonic() >= deadline:
                        break
                    delta = np.eye(4)
                    if axis < 3:
                        delta[axis, 3] = sign * translation
                    else:
                        vector = np.zeros(3); vector[axis - 3] = sign * rotation
                        delta[:3, :3] = Rotation.from_rotvec(vector).as_matrix()
                        delta[:3, 3] = center - delta[:3, :3] @ center
                    trial = delta @ matrix
                    score = evaluate(trial); trials += 1
                    if feasible(score):
                        row = dict(transform=trial.ravel().tolist(), score=score, origin="feasible_trial")
                        archive.append(row)
                        if score["shape"] > incumbent["score"]["shape"] + 1e-9:
                            incumbent = row; improved = True
    return dict(selected=incumbent, archive=archive, trial_evaluations=trials,
                endpoint_replay_seconds=replay_seconds,
                budget_exhausted=time.monotonic() >= deadline,
                original_policy_preserved=any(r["score"]["original_policy_pass"] for r in archive))


def polar_exposure(points, chemistry, feature_types, receptor, receptor_elements):
    from rdkit import Chem
    from .anchor_extraction import atomic_sasa
    from .conformer_artifacts import FEATURE_TYPES
    polar_types = {FEATURE_TYPES[name] for name in ("Donor", "Acceptor", "PosIonizable", "NegIonizable")}
    members = sorted({int(j) for i, kind in enumerate(feature_types) if kind in polar_types
                      for j in chemistry.feature_members[i]})
    elements = [Chem.GetPeriodicTable().GetElementSymbol(int(z)) for z in chemistry.atomic_numbers]
    isolated = atomic_sasa(points, elements, points, elements, target_occluder_indices=range(len(points)))
    bound = atomic_sasa(points, elements, np.vstack((points, receptor)),
                        elements + list(receptor_elements), target_occluder_indices=range(len(points)))
    return dict(method="Heavy-atom Shrake-Rupley; 960 points; 1.4 A probe; no waters",
                isolated_polar_sasa=float(isolated[members].sum()),
                bound_polar_sasa=float(bound[members].sum()),
                atoms=[dict(index=i, element=elements[i], isolated_sasa=float(isolated[i]),
                            bound_sasa=float(bound[i]), protein_buried_sasa=float(isolated[i]-bound[i]))
                       for i in members],
                limitation="Geometric polar exposure is not permeability or hydrogen-bond satisfaction.")
