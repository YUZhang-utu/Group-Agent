"""Crystal-derived 3D query hypotheses; no candidate complex or docking analysis."""
from __future__ import annotations

from collections import defaultdict
import csv
import html
from pathlib import Path
import time
import tempfile

import numpy as np

from . import library_acceptance as ev
from .anchor_extraction import _protein_atoms
from .chemistry_prep import _column, _mmcif_dict
from .expanded_wee1 import fingerprint, ensure_file_descriptor_limit
from .gaussian_batch import _load_query, _atomic_json, prepare_gaussian_query
from .interaction_fast import score_batched
from .interaction_review import OBJECTIVE
from .screening_selection import archive, check_hashes

CLASSES = ("hydrogen_bond", "hydrophobic", "pi_stacking", "cation_pi", "salt_bridge",
           "water_bridge", "metal_coordination", "halogen_bond")
RINGS = {"PHE": ("CG CD1 CE1 CZ CE2 CD2",), "TYR": ("CG CD1 CE1 CZ CE2 CD2",),
         "HIS": ("CG ND1 CE1 NE2 CD2",), "TRP": ("CG CD1 NE1 CE2 CD2", "CD2 CE2 CZ2 CH2 CZ3 CE3")}
HYDROPHOBIC = {"ALA":"CB", "VAL":"CB CG1 CG2", "LEU":"CB CG CD1 CD2", "ILE":"CB CG1 CG2 CD1",
              "MET":"CB", "PRO":"CB CG", "PHE":"CB CG CD1 CD2 CE1 CE2 CZ", "TYR":"CB CG CD1 CD2 CE1 CE2",
              "TRP":"CB CG CD2 CE3 CZ2 CZ3 CH2", "LYS":"CB CG CD", "ARG":"CB CG", "GLU":"CB", "GLN":"CB"}
THRESHOLDS = dict(hydrophobic_distance=4.0, salt_distance=5.5, pi_distance=5.5,
                  pi_angle_deviation=30.0, pi_offset=2.0, cation_pi_distance=6.0,
                  polar_water_distance=3.5, water_protein_distance=3.5, metal_distance=3.0)


def partner(atom):
    return {k: atom[k] for k in ("residue", "chain", "residue_number", "atom_name")}


def crystal_environment(mmcif, query_id):
    from .chemistry_prep import enumerate_ligand_instances
    _, ccd, chain, residue = query_id.split(":")
    instances = [x for x in enumerate_ligand_instances(mmcif,[ccd]) if x["chain_id"] == chain and x["residue_number"] == residue]
    if len(instances) != 1: raise ValueError("Ambiguous crystal ligand instance")
    model = instances[0]["model"]
    atoms = _protein_atoms(mmcif,model,(ccd,chain,residue))
    data = _mmcif_dict(mmcif); n = len(data["_atom_site.group_PDB"])
    names = ["auth_comp_id","auth_asym_id","auth_seq_id","auth_atom_id","type_symbol","Cartn_x","Cartn_y","Cartn_z","pdbx_PDB_model_num"]
    columns = {k:_column(data,"_atom_site."+k,n) for k in names}
    waters, metals = [], []
    for i in range(n):
        if columns["pdbx_PDB_model_num"][i] != model: continue
        symbol, res = columns["type_symbol"][i].upper(), columns["auth_comp_id"][i]
        if not (res in {"HOH","WAT","DOD"} and symbol == "O") and symbol not in {"ZN","MG","CA","MN","FE","CO","NI","CU"}: continue
        row = dict(residue=res,chain=columns["auth_asym_id"][i],residue_number=columns["auth_seq_id"][i],
                   atom_name=columns["auth_atom_id"][i],xyz=np.array([float(columns[k][i]) for k in ("Cartn_x","Cartn_y","Cartn_z")]))
        (waters if res in {"HOH","WAT","DOD"} else metals).append(row)
    return atoms, waters, metals


def extra_hypotheses(query_id, query, atoms, waters=(), metals=()):
    """Transparent geometric hypotheses, not PLIP-equivalent assignments or energies."""
    groups = defaultdict(dict)
    for a in atoms: groups[(a["chain"],a["residue_number"],a["residue"])][a["atom_name"]] = a
    rings, charges = [], []
    for (_, _, res), group in groups.items():
        for names in RINGS.get(res, ()):
            names = names.split()
            if all(n in group for n in names):
                points = np.array([group[n]["xyz"] for n in names]); center = points.mean(axis=0)
                _, _, vh = np.linalg.svd(points-center,full_matrices=False)
                rings.append((center,vh[-1],dict(partner(group[names[0]]),atom_name="/".join(names))))
        definition = {"LYS":(1,["NZ"]), "ARG":(1,["NE","NH1","NH2"]), "ASP":(-1,["OD1","OD2"]), "GLU":(-1,["OE1","OE2"])}.get(res)
        if definition and all(n in group for n in definition[1]):
            sign, names = definition
            charges.append((sign,np.array([group[n]["xyz"] for n in names]).mean(axis=0),dict(partner(group[names[0]]),atom_name="/".join(names))))
    rows = []
    def add(kind, index, target, distance, **geometry):
        key = f"{query_id}/F{index}:{kind}:{target['chain']}:{target['residue_number']}:{target['atom_name']}"
        rows.append(dict(anchor_id=key,feature_index=index,feature_class=kind,query_id=query_id,
            ligand_atom_indices=[],atom_center=query["feature_points"][index].tolist(),weight=1.0,
            evidence_class="crystal_derived_geometric_hypothesis",
            evidence=dict(interaction=kind,protein_partner=target,distance_angstrom=float(distance),
                          angle_degrees=geometry.pop("angle_degrees",None),angle_status="not_a_validated_hydrogen_bond",**geometry),
            interpretation="Match the crystal ligand feature, not a candidate-protein contact; shared feature columns are not independent evidence"))
    for index, (point, kind) in enumerate(zip(query["feature_points"],query["feature_types"])):
        if kind == 5:
            nearest = {}
            for atom in atoms:
                if atom["atom_name"] not in HYDROPHOBIC.get(atom["residue"],"").split(): continue
                d = np.linalg.norm(point-atom["xyz"])
                key = (atom["chain"],atom["residue_number"],atom["residue"])
                if d <= THRESHOLDS["hydrophobic_distance"] and (key not in nearest or d < nearest[key][0]): nearest[key] = (d,atom)
            for d, atom in nearest.values(): add("hydrophobic",index,partner(atom),d)
        if kind in (3,4):
            for sign, center, target in charges:
                d = np.linalg.norm(point-center)
                if sign == (-1 if kind == 3 else 1) and d <= THRESHOLDS["salt_distance"]:
                    add("salt_bridge",index,target,d,charge_assumption="RDKit ionizable ligand feature and residue template; protonation unvalidated")
        if kind == 6:
            normal = query["feature_directions"][index]
            for center, pn, target in rings:
                delta = point-center; d = np.linalg.norm(delta)
                angle = float(np.degrees(np.arccos(np.clip(abs(np.dot(normal,pn)),0,1))))
                offset = min(np.linalg.norm(delta-np.dot(delta,pn)*pn),np.linalg.norm(delta-np.dot(delta,normal)*normal))
                if query["feature_direction_kinds"][index] == 2 and d <= 5.5 and (angle <= 30 or angle >= 60) and offset <= 2:
                    add("pi_stacking",index,target,d,angle_degrees=angle,offset_angstrom=float(offset))
            for sign, center, target in charges:
                d = np.linalg.norm(point-center)
                if sign == 1 and d <= 6: add("cation_pi",index,target,d,scope="distance hypothesis; protonation and face geometry need review")
        if kind == 3:
            for center, pn, target in rings:
                d = np.linalg.norm(point-center)
                if d <= 6: add("cation_pi",index,target,d,scope="distance hypothesis; protonation and face geometry need review")
        if kind in (1,2):
            for water in waters:
                dw = np.linalg.norm(point-water["xyz"])
                if dw > 3.5: continue
                for atom in atoms:
                    dp = np.linalg.norm(atom["xyz"]-water["xyz"])
                    if dp <= 3.5 and (atom["donor"] or atom["acceptor"]):
                        add("water_bridge",index,partner(atom),dw,water=partner(water),water_protein_distance=float(dp),scope="two-distance hypothesis; water orientation unresolved")
            if kind == 2:
                for metal in metals:
                    d = np.linalg.norm(point-metal["xyz"])
                    if d <= 3: add("metal_coordination",index,partner(metal),d,scope="proximity hypothesis; metal-specific coordination chemistry unvalidated")
    # Multiple waters can support the same feature/residue hypothesis. Keep the nearest.
    unique = {}
    for row in sorted(rows,key=lambda r:r["evidence"]["distance_angstrom"]): unique.setdefault(row["anchor_id"],row)
    return list(unique.values())


def write_tables(report, output):
    records = []
    for q in report["queries"]:
        for a in q["anchors"]:
            e=a["evidence"]; p=e.get("protein_partner",{})
            records.append(dict(query=q["query_id"],category=a.get("interaction_class","hydrogen_bond"),
                anchor_id=a["anchor_id"],feature_index=a["feature_index"],
                protein_partner=" ".join(str(p.get(k,"")) for k in ("chain","residue","residue_number","atom_name")),
                distance_angstrom=e.get("distance_angstrom"),angle_degrees=e.get("angle_degrees"),
                ligand_atoms=str(a['ligand_atom_indices']),molecules_at_score_0_5=next((d['molecules'] for d in a['diagnostic_counts'] if d['minimum_score']==.5),None),
                evidence=a["evidence_class"],interpretation=a["interpretation"]))
    fields=["query","category","anchor_id","feature_index","ligand_atoms","protein_partner","distance_angstrom","angle_degrees","molecules_at_score_0_5","evidence","interpretation"]
    with (output/"interaction-classes.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(records)
    pieces=['<!doctype html><html lang="en"><meta charset="utf-8"><title>Crystal-derived 3D search evidence</title>',
        '<style>body{font:15px system-ui;margin:32px;color:#16324a}table{border-collapse:collapse;width:100%}td,th{padding:9px;border:1px solid #ccd6df;text-align:left;overflow-wrap:anywhere}th{background:#eaf1f7}h2{margin-top:32px}</style>',
        '<h1>Crystal-derived 3D search evidence</h1><p>Reference geometry and candidate feature matching. No docking or candidate-complex interaction validation.</p>',
        '<p>Choose an anchor ID, then request a same-pose selection preview in chat. A feature shared by several contacts is not independent evidence for each contact.</p>']
    if report.get('kind') == 'full_library_conditions':
        pieces.append('<p><strong>Coverage status: '+html.escape(report['status'])+'</strong>. '
            'Counts concern one heuristic pose per evaluated conformer; no Top-K or Top-N membership limit.</p>')
        for q in report['queries']:
            pieces.append('<p>'+html.escape(q['query_id'])+': '+str(q['counts']['evaluated_conformers'])+
                ' / '+str(q['counts']['total_conformers'])+' conformers evaluated.</p>')
    for kind in CLASSES:
        matching=[r for r in records if r['category']==kind]
        pieces.append('<h2>'+html.escape(kind.replace('_',' ').title())+'</h2>')
        if not matching:
            pieces.append('<p>'+('Not indexed: halogen-specific ligand features require a separate validated extension.' if kind=='halogen_bond' else 'No hypothesis detected by these rules in the supplied crystal. This does not establish chemical absence.')+'</p>');continue
        pieces.append('<table><thead><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in fields)+'</tr></thead><tbody>')
        for r in matching: pieces.append('<tr>'+''.join('<td>'+html.escape(str(r[k]))+'</td>' for k in fields)+'</tr>')
        pieces.append('</tbody></table>')
    pieces.append('</html>');(output/'interactions.html').write_text('\n'.join(pieces),encoding='utf-8')


def classify(review_path, output):
    started=time.perf_counter(); original=ev.read(review_path); check_hashes(original['sources'])
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    sources=dict(original['sources']);sources.update(fingerprint([review_path]));queries=[]
    for q in original['queries']:
        q=dict(q)
        query_path=Path(q.get('query_npz') or ev.read(Path(q['sidecar']).with_suffix('.manifest.json'))['inputs']['query']['path'])
        manifest, arrays=_load_query(query_path)
        members=manifest['source'].get('feature_atom_indices')
        if members is None:
            src=manifest['source']
            with tempfile.TemporaryDirectory() as tmp:
                rebuilt=Path(tmp)/'query.npz'
                rebuilt_manifest=prepare_gaussian_query(Path(src['mmcif']),Path(src['ccd']),Path(src['query_manifest']),rebuilt)
                rebuilt_arrays=archive(rebuilt)
                if set(arrays)!=set(rebuilt_arrays) or any(not np.array_equal(v,rebuilt_arrays[k]) for k,v in arrays.items()):
                    raise ValueError('Legacy feature atom mapping reconstruction differs')
                members=rebuilt_manifest['source']['feature_atom_indices']
        atoms,waters,metals=crystal_environment(Path(manifest['source']['mmcif']),q['query_id'])
        anchors=[dict(a,interaction_class='hydrogen_bond') for a in q['anchors']]
        extras=extra_hypotheses(q['query_id'],arrays,atoms,waters,metals)
        for a in extras: a['ligand_atom_indices']=members[a['feature_index']]
        anchors.extend(dict(a,interaction_class=a['feature_class']) for a in extras)
        indices=sorted({a['feature_index'] for a in anchors})
        for a in anchors: a['score_column']=indices.index(a['feature_index'])
        expanded=dict(arrays,anchor_feature_indices=np.array(indices,dtype=np.int64))
        r=archive(q['rigid'])
        import os
        if os.name=='posix': ensure_file_descriptor_limit(len(ev.read(q['artifact_catalog'])['shards']))
        ev.log('Classifying 3D feature matches for '+q['query_id'])
        scores, assignments, contributions=score_batched(Path(q['artifact_catalog']),Path(q['chemical_companion']),expanded,
            r['global_ids'],r['molecule_ids'],r['conformer_ids'],[OBJECTIVE],{OBJECTIVE:r[OBJECTIVE+'__transform']},
            sigma=1.,cutoff=4.5,angular_power=2.)
        side=output/(q['query_id'].split(':')[0].lower()+'-classified-matches.npz')
        np.savez(side,**{k:r[k] for k in ('global_ids','molecule_ids','conformer_ids')},query_anchor_feature_indices=np.array(indices),
            **{OBJECTIVE+'__anchor_scores':contributions[OBJECTIVE],OBJECTIVE+'__anchor_assignments':assignments[OBJECTIVE],
               OBJECTIVE+'__interaction_match_score':scores[OBJECTIVE]})
        for a in anchors:
            values=contributions[OBJECTIVE][:,a['score_column']];assigned=assignments[OBJECTIVE][:,a['score_column']]>=0
            a['diagnostic_counts']=[dict(minimum_score=t,conformers=int(np.sum(assigned & (values>=t))),
                molecules=len(set(map(str,r['molecule_ids'][assigned & (values>=t)])))) for t in (.25,.5,.75)]
        q.update(anchors=anchors,sidecar=str(side),query_npz=str(query_path));queries.append(q)
        sources.update(fingerprint([side,query_path,query_path.with_suffix('.manifest.json')]))
    result=dict(original,queries=queries,sources=sources,classification='native-crystal-feature-hypotheses-v1',
        unsupported_classes=['halogen_bond'],thresholds=THRESHOLDS,wall_seconds=time.perf_counter()-started,
        policy='Expanded query feature matching on original rigid poses. Original E031 outputs unchanged; new scores are not claimed equivalent to E031.')
    result['classes']={kind:[a['anchor_id'] for q in queries for a in q['anchors'] if a['interaction_class']==kind] for kind in CLASSES}
    check_hashes(sources);_atomic_json(output/'report.json',result);write_tables(result,output)
    return result
