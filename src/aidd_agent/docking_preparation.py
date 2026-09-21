"""Crystal-pocket-derived preparation and gated Glide execution, after 3D selection."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

from . import library_acceptance as ev
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json, _load_query
from .screening_selection import check_hashes


def pocket_geometry(points, padding=5.):
    points=np.asarray(points,dtype=float)
    if points.ndim!=2 or points.shape[1]!=3 or not len(points) or not np.isfinite(points).all():
        raise ValueError('Invalid crystal ligand coordinates')
    center=(points.min(axis=0)+points.max(axis=0))/2
    return dict(center_angstrom=center.tolist(),
                radius_angstrom=float(np.ceil(np.linalg.norm(points-center,axis=1).max()+padding)),
                outer_box_angstrom=np.maximum(20.,np.ceil(np.ptp(points,axis=0)+2*padding)).tolist(),
                inner_box_angstrom=[10.,10.,10.],padding_angstrom=padding,
                policy='Axis-aligned reference-ligand bounding box midpoint; dimensions require review for larger candidates')


def profile_path():
    return Path(os.environ.get('AIDD_DOCKING_PROFILE',
        str(Path(os.environ.get('AIDD_CONFIG_DIR','/mnt/local/hand/yuzhang/aidd/config'))/'docking.local.json'))).expanduser().resolve()


def prepare(handoff_path, output):
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from rdkit import Chem
    from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances
    handoff=ev.read(handoff_path)
    if handoff.get('kind')!='docking_handoff' or handoff['counts']['selected_molecules']<1:
        raise ValueError('Docking preparation needs a nonempty exported selection')
    check_hashes(handoff['sources']);check_hashes(handoff['outputs'])
    q=handoff['query'];query_path=Path(q['query_npz'])
    manifest,arrays=_load_query(query_path);source=manifest['source']
    for key in ('mmcif','ccd','query_manifest'):
        if ev.sha(source[key])!=source[key+'_sha256']:raise ValueError('Crystal source changed')
    geometry=pocket_geometry(arrays['shape_points']);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    _,ccd,chain,residue=q['query_id'].split(':')
    structure=MMCIFParser(QUIET=True).get_structure('crystal',source['mmcif'])
    first=next(structure.get_models()); removed=[]
    class Receptor(Select):
        def accept_model(self, model): return int(model.id==first.id)
        def accept_residue(self, item):
            reject=item.get_parent().id==chain and str(item.id[1])==residue and item.resname==ccd
            if reject:removed.append(item)
            return int(not reject)
    writer=PDBIO();writer.set_structure(structure);writer.save(str(output/'receptor-input.pdb'),Receptor())
    if len(removed)!=1:raise ValueError('Reference ligand removal was ambiguous')
    instances=[r for r in enumerate_ligand_instances(Path(source['mmcif']),[ccd]) if r['chain_id']==chain and r['residue_number']==residue]
    if len(instances)!=1:raise ValueError('Reference ligand is ambiguous')
    reference=_ccd_molecule(Path(source['ccd']),instances[0]['atoms'])
    with Chem.SDWriter(str(output/'reference-native.sdf')) as sd:sd.write(reference)
    shutil.copyfile(handoff['sdf'],output/'selected-poses.sdf')
    shutil.copyfile(handoff['ids'],output/'selected-ids.json')
    _atomic_json(output/'pocket.json',dict(query_id=q['query_id'],**geometry))
    config=ev.read(profile_path())
    root=Path(config['schrodinger_root']).expanduser()
    grid='\n'.join(['RECEP_FILE prepared-receptor.mae','GRIDFILE receptor-grid.zip',
        'GRID_CENTER '+', '.join(map(str,geometry['center_angstrom'])),
        'INNERBOX 10, 10, 10','OUTERBOX '+', '.join(map(str,geometry['outer_box_angstrom']))])+'\n'
    (output/'grid.in').write_text(grid,encoding='utf-8')
    for name,ligand in (('reference','prepared-reference.maegz'),('candidates','prepared-ligands.maegz')):
        (output/(name+'.in')).write_text('GRIDFILE receptor-grid.zip\nLIGANDFILE '+ligand+'\nPRECISION SP\nPOSES_PER_LIG 5\n',encoding='utf-8')
    commands=[dict(name='protein_preparation',argv=[str(root/'utilities/prepwizard'),'receptor-input.pdb','prepared-receptor.mae','-WAIT'],expected='prepared-receptor.mae'),
        dict(name='reference_preparation',argv=[str(root/'ligprep'),'-isd','reference-native.sdf','-omae','prepared-reference.maegz','-WAIT'],expected='prepared-reference.maegz'),
        dict(name='grid_generation',argv=[str(root/'glide'),'-WAIT','grid.in'],expected='receptor-grid.zip'),
        dict(name='reference_redocking',argv=[str(root/'glide'),'-WAIT','reference.in'],expected='reference_pv.maegz'),
        dict(name='reference_conversion',argv=[str(root/'utilities/structconvert'),'-imae','reference_pv.maegz','-osd','reference-poses.sdf'],expected='reference-poses.sdf'),
        dict(name='ligand_preparation',argv=[str(root/'ligprep'),'-isd','selected-poses.sdf','-omae','prepared-ligands.maegz','-WAIT'],expected='prepared-ligands.maegz'),
        dict(name='candidate_docking',argv=[str(root/'glide'),'-WAIT','candidates.in'],expected='candidates_pv.maegz')]
    result=dict(status='complete',kind='docking_preparation',engine='glide',query_id=q['query_id'],geometry=geometry,
        commands=commands,redocking_rmsd_threshold_angstrom=2.,
        threshold_scope='Provisional reference-pose engineering gate, not enrichment validation',
        preparation_status='inputs_and_commands_prepared_not_executed',
        sources={**handoff['sources'],**handoff['outputs'],**fingerprint([handoff_path,profile_path()])},
        policy='Selected ligand removed; other residues, waters and cofactors retained for PrepWizard review. Default installed preparation chemistry; inspect protonation, stereochemistry and retained cofactors before execution.')
    result['outputs']=fingerprint(p for p in output.iterdir() if p.is_file() and p.name!='report.json')
    _atomic_json(output/'report.json',result)
    return result


def reference_rmsd(native, predicted):
    """Symmetry-aware RMSD in receptor coordinates, without pose alignment."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    refs=[m for m in Chem.SDMolSupplier(str(native),removeHs=True) if m is not None]
    if len(refs)!=1:raise ValueError('Expected one native reference ligand')
    ref=refs[0];values=[]
    for mol in Chem.SDMolSupplier(str(predicted),removeHs=True):
        if mol is None or mol.GetNumAtoms()!=ref.GetNumAtoms() or not mol.HasSubstructMatch(ref):continue
        values.append(float(rdMolAlign.CalcRMS(mol,ref)))
    if not values or not np.isfinite(values).all():raise ValueError('No chemically compatible reference pose for RMSD validation')
    return dict(best_in_place_heavy_atom_rmsd=min(values),compatible_poses=len(values),
                policy='Best of returned poses, symmetry-aware, no coordinate alignment; not a top-ranked-pose acceptance claim')


def run(preparation_path, output):
    plan=ev.read(preparation_path)
    if plan.get('kind')!='docking_preparation':raise ValueError('Expected a docking preparation')
    check_hashes(plan['sources']);check_hashes(plan['outputs'])
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    if any(output.iterdir()):raise ValueError('Use a fresh docking execution task; partial engine jobs are not silently reused')
    for name in plan['outputs']:shutil.copyfile(name,output/Path(name).name)
    result=dict(status='running',kind='docking_execution',steps=[],sources={**plan['sources'],**plan['outputs'],**fingerprint([preparation_path])},
                biological_quality='not_evaluated',reference_validation=None)
    try:
        for step in plan['commands']:
            ev.log('Docking workflow: '+step['name'])
            executable=Path(step['argv'][0])
            if not executable.is_file() or not os.access(executable,os.X_OK):raise ValueError('Engine executable unavailable: '+str(executable))
            with (output/(step['name']+'.log')).open('w',encoding='utf-8') as log:
                process=subprocess.run(step['argv'],cwd=output,stdout=log,stderr=subprocess.STDOUT,shell=False)
            if process.returncode or not (output/step['expected']).is_file() or not (output/step['expected']).stat().st_size:
                raise ValueError('External stage failed or expected output missing: '+step['name'])
            result['steps'].append(dict(name=step['name'],status='complete',output_sha256=ev.sha(output/step['expected'])))
            if step['name']=='reference_conversion':
                validation=reference_rmsd(output/'reference-native.sdf',output/'reference-poses.sdf')
                validation['passed']=validation['best_in_place_heavy_atom_rmsd']<=plan['redocking_rmsd_threshold_angstrom']
                result['reference_validation']=validation
                if not validation['passed']:raise ValueError('Reference redocking gate failed; candidate docking not submitted')
            _atomic_json(output/'report.json',result)
        result['status']='complete'
        result['pose_quality']='reference_redocking_gate_only_candidate_poses_unvalidated'
        result['outputs']=fingerprint(p for p in output.iterdir() if p.is_file() and p.name!='report.json')
    except Exception as exc:
        result['status']='failed';result['error']=str(exc);_atomic_json(output/'report.json',result);raise
    _atomic_json(output/'report.json',result)
    return result
