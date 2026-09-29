"""Experimental receptor advice; representative sets are not physical states."""
from collections import Counter
import hashlib
import json

import numpy as np


DEFAULTS = dict(minimum_pdbs=4, minimum_series=2, minimum_repeat_entries=3,
                coverage_fraction=.95, fallback_tolerance=.15,
                maximum_representatives=32, minimum_plateau_width=.08,
                minimum_cluster_pdbs=3)


def fuse(overlap, local, chemical, chemical_weight=.15):
    """Freeze fitted scales in each report; never rank-transform small differences."""
    arrays = [1-np.asarray(overlap, float), np.asarray(local, float)]
    scales = []
    for values in arrays:
        pairs = values[np.triu_indices(len(values), 1)]
        scales.append(max(.1, float(np.quantile(pairs, .9))) if len(pairs) else 1.)
    # Bounded, monotone, zero-preserving; the fitted 90th percentile maps to 0.5.
    normalized = [a/(a+s) for a, s in zip(arrays, scales)]
    shape = .5*normalized[0] + .5*normalized[1]
    distance = (1-chemical_weight)*shape + chemical_weight*np.asarray(chemical)
    np.fill_diagonal(distance, 0)
    return distance, dict(version='robust-weighted-v1', scales=scales,
                         shape_weights=[.5, .5], chemical_weight=chemical_weight,
                         normalization='x/(x+max(0.1,pair_p90)); frozen per report',
                         scope='Exploratory cohort calibration; not transferable cutoffs')


def ligand_series(rows):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    series, identities, missing, mapping = set(), set(), set(), {}
    queries = {q['query_id']: q for r in rows for q in r.get('queries', [])}
    for qid, query in sorted(queries.items()):
        mol = Chem.MolFromSmiles(query.get('smiles') or '')
        if mol is None or not mol.GetNumAtoms():
            missing.add(qid)
            continue
        identities.add(Chem.MolToSmiles(mol, isomericSmiles=True))
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if not scaffold.GetNumAtoms():
            missing.add(qid)
            continue
        key = Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(scaffold))
        sid = 'series-'+hashlib.sha256(key.encode()).hexdigest()[:12]
        series.add(sid)
        mapping[qid] = sid
    return dict(independent_series_proxy=len(series), unique_ligand_identities=len(identities),
                unclassified_ligand_instances=len(missing), query_series=mapping,
                method='Generic Murcko ring-scaffold proxy; acyclic/missing structures unclassified',
                limitation='Scaffold count is not verified medicinal-chemistry series independence')


def partition_signature(groups):
    return tuple(sorted(tuple(sorted(g['members'])) for g in groups))


def stable_partition(distance, rows, tolerance, settings):
    from .pocket_states import clusters
    if len(rows) < 2*settings['minimum_cluster_pdbs']:
        return [], []
    upper = float(distance.max())
    cuts = np.linspace(0, upper, 41)
    scans = [(float(t), clusters(distance, rows, float(t))) for t in cuts]
    intervals = []
    start = 0
    for end in range(1, len(scans)+1):
        if end < len(scans) and partition_signature(scans[end][1]) == partition_signature(scans[start][1]):
            continue
        lo, groups = scans[start]
        hi = scans[end-1][0]
        if 1 < len(groups) < len(rows):
            supported = all(g['distinct_pdb_support'] >= settings['minimum_cluster_pdbs'] for g in groups)
            intervals.append(dict(lower=lo, upper=hi, width=hi-lo, k=len(groups), supported=supported,
                                  groups=groups, cutoff=(lo+hi)/2))
        start = end
    candidates = sorted(intervals, key=lambda x: (-x['width'], x['k']))
    for interval in candidates:
        interval['leave_entry_out_stable'] = False
        if not interval['supported'] or interval['width'] < max(settings['minimum_plateau_width'], tolerance*.5):
            continue
        groups = interval['groups']
        cross = [distance[i,j] for g in groups for i,r in enumerate(rows) if r['id'] in g['members']
                 for j,s in enumerate(rows) if s['id'] not in g['members']]
        if not cross or min(cross) <= tolerance:
            continue
        stable = True
        # Entry-level controls: do not treat chains of one crystal as independent.
        for entry in sorted({r['pdb_id'] for r in rows}):
            ix = [i for i,r in enumerate(rows) if r['pdb_id'] != entry]
            ids = {rows[i]['id'] for i in ix}
            expected = tuple(sorted(tuple(sorted(set(g['members']) & ids)) for g in groups
                                    if set(g['members']) & ids))
            actual = clusters(distance[np.ix_(ix,ix)], [rows[i] for i in ix], interval['cutoff'])
            if partition_signature(actual) != expected:
                stable = False
                break
        interval['leave_entry_out_stable'] = stable
        if stable:
            return groups, [{k:v for k,v in x.items() if k != 'groups'} for x in intervals]
    return [], [{k:v for k,v in x.items() if k != 'groups'} for x in intervals]


def recommend(rows, distance, settings=None):
    p = dict(DEFAULTS, **(settings or {}))
    if set(p) != set(DEFAULTS):
        raise ValueError('Unknown receptor advice setting')
    if any(type(v) not in (int,float) or not np.isfinite(v) or v <= 0 for v in p.values()):
        raise ValueError('Invalid receptor advice settings')
    if p['coverage_fraction'] > 1 or any(type(p[k]) is not int for k in (
            'minimum_pdbs','minimum_series','minimum_repeat_entries','maximum_representatives','minimum_cluster_pdbs')):
        raise ValueError('Invalid receptor advice bounds')
    d = np.asarray(distance, float)
    n = len(rows)
    if not n or d.shape != (n,n) or not np.isfinite(d).all() or (d < 0).any() or not np.allclose(d,d.T) or not np.allclose(d.diagonal(),0):
        raise ValueError('Invalid receptor advice matrix')
    if len({r['id'] for r in rows}) != n:
        raise ValueError('Duplicate receptor IDs')
    counts = Counter(r['pdb_id'] for r in rows)
    weights = np.array([1/counts[r['pdb_id']] for r in rows]); weights /= weights.sum()
    chemistry = ligand_series(rows)
    repeats = {entry: [float(d[i,j]) for i in range(n) for j in range(i)
                       if rows[i]['pdb_id'] == rows[j]['pdb_id'] == entry]
               for entry in counts if counts[entry] > 1}
    # Each entry has one contribution, irrespective of the number of chains.
    repeat_values = [float(np.quantile(v,.9)) for v in repeats.values()]
    sufficient_repeats = len(repeat_values) >= p['minimum_repeat_entries']
    repeat_q90 = float(np.quantile(repeat_values,.9)) if repeat_values else None
    tolerance = max(.02, repeat_q90) if sufficient_repeats else p['fallback_tolerance']
    costs = d @ weights
    reference = min(range(n), key=lambda i:(not bool(rows[i].get('queries')), float(costs[i]),
                                          rows[i].get('resolution') or 99, rows[i]['id']))
    pairs = d[np.triu_indices(n,1)]
    pair_q90 = float(np.quantile(pairs,.9)) if len(pairs) else 0.
    sparse = len(counts) < p['minimum_pdbs'] or chemistry['independent_series_proxy'] < p['minimum_series']
    # Require all observations to lie within the provisional repeat envelope;
    # a rare but distinct pocket must not disappear behind a bulk quantile.
    unresolved = sufficient_repeats and float(d.max()) <= tolerance
    partitions, plateau = ([], []) if sparse or unresolved else stable_partition(d, rows, tolerance, p)
    partition_over_budget = len(partitions) > p['maximum_representatives']
    if partition_over_budget:
        partitions = []
    if sparse:
        mode, reason = 'single_insufficient_evidence', 'Limited independent structure or ligand-series evidence; start with one reference.'
    elif unresolved:
        mode, reason = 'single_unresolved_variation', 'Observed differences fall within the provisional near-repeat envelope; start with one receptor.'
    elif partitions:
        mode, reason = 'partition_representatives', 'A supported partition is stable across cutoffs and leave-entry-out checks; review its representatives.'
    else:
        mode, reason = 'coverage_representatives', 'No robust partition was resolved; use a representative coverage set without asserting physical continuity.'
    if partition_over_budget:
        reason = 'A stable partition exceeds the requested receptor budget; provide a limited coverage proposal for user review.'
    selected = [next(i for i,r in enumerate(rows) if r['id'] == g['representative']) for g in partitions] if partitions else [reference]
    curve = []
    while True:
        residual = d[:,selected].min(axis=1)
        covered = float(weights[residual <= tolerance].sum())
        curve.append(dict(k=len(selected), weighted_coverage=covered, maximum_remaining_distance=float(residual.max()),
                          weighted_mean_remaining_distance=float(residual @ weights)))
        if mode != 'coverage_representatives' or covered >= p['coverage_fraction']-1e-12 or len(selected) >= min(n,p['maximum_representatives']):
            break
        options = [i for i in range(n) if i not in selected]
        # Maximize newly covered PDB-weighted support; residual improvement breaks ties.
        best = min(options, key=lambda i:(-float(weights[(residual > tolerance)&(d[:,i] <= tolerance)].sum()),
                        float(np.minimum(residual,d[:,i]) @ weights), not bool(rows[i].get('queries')), rows[i]['id']))
        selected.append(best)
    residual = d[:,selected].min(axis=1)
    proposals = []
    for position,i in enumerate(selected):
        if partitions:
            members = partitions[position]['members']
        else:
            # Overlapping local neighborhoods, not Voronoi bins passed off as states.
            members = [rows[j]['id'] for j in range(n) if d[j,i] <= tolerance]
        rid = 'pocket-'+hashlib.sha256(json.dumps([mode,rows[i]['id'],sorted(members)]).encode()).hexdigest()[:12]
        newly_covered = d[:,i] <= tolerance
        if position:
            newly_covered &= d[:,selected[:position]].min(axis=1) > tolerance
        gain = float(weights[newly_covered].sum())
        proposals.append(dict(id=rid, representative=rows[i]['id'], members=members,
            distinct_pdb_support=len({r['pdb_id'] for r in rows if r['id'] in members}),
            small_support=len({r['pdb_id'] for r in rows if r['id'] in members})<3,
            role='partition_representative' if partitions else 'local_coverage_reference',
            consensus_ready=bool(rows[i].get('queries')),
            added_weighted_coverage=gain,
            reason='Represents a stable structural partition.' if partitions else
                   f'Covers {len(members)} nearby chains; adds {gain:.1%} PDB-weighted coverage.'))
    return dict(version=1, mode=mode, recommendation=reason, recommended_k=len(selected),
        evidence=dict(qualified_chains=n, independent_pdbs=len(counts), **chemistry),
        repeat_baseline=dict(entries=len(repeats), pairs=sum(map(len,repeats.values())), entry_balanced_q90=repeat_q90,
                             sufficient=sufficient_repeats, pair_distance_q90=pair_q90,
                             interpretation='Experimental near-repeat variation, not a pure noise floor'),
        coverage=dict(tolerance=tolerance, tolerance_source='near_repeat_proxy' if sufficient_repeats else 'uncalibrated_engineering_default',
                      target=p['coverage_fraction'], achieved=curve[-1]['weighted_coverage'],
                      target_met=curve[-1]['weighted_coverage'] >= p['coverage_fraction']-1e-12,
                      curve=curve, uncovered=[rows[i]['id'] for i in range(n) if residual[i] > tolerance]),
        representatives=proposals, policy=p, plateau_diagnostics=plateau,
        literature=dict(status='not_reviewed', instruction='Use the domain agent literature_search tool for target-specific publications before adding cited ideas.'),
        user_decision='Accept, select a subset, or revise the representative budget; no automatic adoption.',
        limitations=['Engineering recommendation, not proof of rigidity or physical states.',
                     'Grid/alignment perturbation calibration remains unperformed.',
                     'Scaffold series are proxies; unknown series reduce evidence, not receptor eligibility.',
                     'PDB weighting corrects chain multiplicity, not sampling bias or equilibrium occupancy.',
                     'Uncovered structures remain available for explicit user selection.'])


def attach(report, settings=None):
    """Attach advice to a fresh report; preserve the old partition as diagnostics."""
    if report.get('receptor_advice'):
        raise ValueError('Advice is already attached; use the original source report')
    d, calibration = fuse(report['pairwise']['overlap'], report['pairwise']['local_difference'],
                          report['pairwise']['chemical_distance'], report['policy']['chemical_weight'])
    report['pairwise']['legacy_distance'] = report['pairwise']['distance']
    report['pairwise']['distance'] = d.tolist()
    report['metric_calibration'] = calibration
    report['diagnostic_clusters'] = report['clusters']
    report['diagnostic_sensitivity'] = report.pop('sensitivity', [])
    advice = recommend(report['structures'], d, settings)
    report['receptor_advice'] = advice
    report['clusters'] = advice['representatives']
    # Every qualified structure remains selectable even when outside the recommendation.
    report['receptor_options'] = [dict(id='pocket-'+hashlib.sha256(('individual:'+r['id']).encode()).hexdigest()[:12],
        representative=r['id'], members=[r['id']], role='user_selected_single_structure',
        distinct_pdb_support=1, small_support=True, consensus_ready=bool(r.get('queries')))
        for r in report['structures']]
    report['descriptor_version'] = 3
    return report


def summary(advice):
    e = advice['evidence']; coverage = advice['coverage']
    lines = ['Receptor selection recommendation', advice['recommendation'],
             f"Evidence: {e['qualified_chains']} qualified chains, {e['independent_pdbs']} independent PDB entries, "
             f"{e['independent_series_proxy']} scaffold-series proxies ({e['unclassified_ligand_instances']} ligand instances unclassified).",
             f"Recommended receptors: {advice['recommended_k']}."]
    if not coverage['target_met']:
        lines.append('The coverage target was not reached; this is a provisional limited set, not an adequate-coverage claim.')
    lines += [f"- {r['representative']}: {r['reason']} Selection ID: {r['id']}."+
              (' Ligand preparation is required before consensus.' if not r['consensus_ready'] else '') for r in advice['representatives'][:5]]
    if len(advice['representatives']) > 5:
        lines.append(f"{len(advice['representatives'])-5} additional representatives are listed in the detailed report.")
    budget = [f"{point['k']} receptors: {point['weighted_coverage']:.1%}" for point in coverage['curve']
              if point['k'] in {1,4,8,advice['recommended_k']}]
    if len(budget) > 1:
        lines.append('Coverage alternatives: '+ '; '.join(budget)+'.')
    lines += [f"PDB-weighted coverage: {coverage['achieved']:.1%}; {len(coverage['uncovered'])} structures remain outside the selected neighborhoods.",
              'Coverage tolerance: '+coverage['tolerance_source'].replace('_',' ')+'.',
              'Literature has not yet been reviewed; request target-specific cited interpretation from the research agent.',
              'Choose the proposed set, a subset, or another qualified structure before proceeding to local interaction consensus.']
    return '\n\n'.join(lines)


def compact(advice):
    """Keep the default agent context small; full diagnostics remain in the report."""
    value = dict(advice)
    value.pop('plateau_diagnostics', None)
    value['evidence'] = {k:v for k,v in advice['evidence'].items() if k != 'query_series'}
    return value
