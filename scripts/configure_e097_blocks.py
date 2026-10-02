"""Create workstation-specific sampling and prompt runtime profiles without editing inputs."""
import argparse
import json
import os
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',required=True,help='Completed blocks-20000 directory')
    p.add_argument('--output',required=True,help='New configuration directory')
    p.add_argument('--runtime',default=os.environ.get('AIDD_RUNTIME_PROFILE'),help='Existing chat runtime to copy, preserving search/AF3 settings')
    p.add_argument('--plants-profile',help='Existing reviewed PLANTS receptor/site profile')
    p.add_argument('--spores',default='~/pdb/1/SPORES_64bit',help='Installed SPORES executable')
    p.add_argument('--plants',help='Installed PLANTS executable; otherwise discover PATH/env/home installation')
    a=p.parse_args(); base=Path(a.base).resolve(); out=Path(a.output).resolve(); out.mkdir(parents=True,exist_ok=True)
    cfg=dict(E094=str(base/'final-blocks-e094/backbone'),E095=str(base/'final-blocks-e095-properties'),
             E096=str(base/'final-blocks-e096-joint/blocks'),profiles=str(base/'final-blocks-e096-joint/profiles'),
             descriptors=str(base/'backbone'),validation=str(base/'final-blocks-e094/properties-validation.json'))
    for key,value in cfg.items():
        path=Path(value) if key=='validation' else Path(value)/'report.json'
        if not path.is_file(): raise SystemExit('Missing '+str(path))
    runtime=json.loads(Path(a.runtime).read_text(encoding='utf-8')) if a.runtime else {}
    # Relative paths in copied profiles must remain anchored to the old config.
    if a.runtime:
        old=Path(a.runtime).resolve().parent
        for container,keys in [(runtime,['af3_profile']), (runtime.get('search',{}),['batch','e034','budget_regions']),
                               (runtime.get('pymol',{}),['executable'])]:
            for key in keys:
                if container.get(key) and not Path(container[key]).is_absolute(): container[key]=str((old/container[key]).resolve())
    runtime['block_evaluation']=dict(sampling_profile=str(out/'sampling.json'))
    if a.plants_profile: runtime['block_evaluation']['plants_profile']=str(Path(a.plants_profile).resolve())
    plants=a.plants
    if not plants:
        from aidd_agent.docking_discovery import discover
        commands=discover(os.environ.get('AIDD_DOCKING_PROFILE'))['commands']
        plants=commands.get('PLANTS') or commands.get('plants')
    artifacts=[('sampling.json',cfg)]
    if plants:
        tools=dict(spores_executable=str(Path(a.spores).expanduser().resolve()),spores_mode='complete',
                   plants_executable=str(Path(plants).expanduser().resolve()),padding_angstrom=5,workers=4,
                   timeout_seconds=600,ligand_chemistry_reviewed=False)
        artifacts.append(('receptor-tools.json',tools))
        runtime['block_evaluation']['receptor_tools_profile']=str(out/'receptor-tools.json')
    else: print('PLANTS executable not discovered. Sampling is ready; later supply --plants in a new configuration directory.')
    artifacts.append(('runtime.json',runtime))
    for name,data in artifacts:
        path=out/name
        if path.exists() and json.loads(path.read_text(encoding='utf-8'))!=data:
            raise SystemExit('Different existing config; choose a new --output directory: '+str(path))
        path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
        print(path)


if __name__=='__main__': main()
