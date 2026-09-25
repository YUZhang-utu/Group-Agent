"""Read-only evidence tools and adapters to existing project-scoped executors."""
import csv
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import time
from urllib.parse import urlencode
from urllib.request import Request, build_opener

TOOLS = {
    'structure_workflow': 'Inspect PDB/diversity/consensus/library stages, source branches, real reference proposals and next supported intents. Arguments: optional task_id. Read-only; call before continuing this workflow.',
    'tasks': 'List current-session task IDs, state and report availability. Arguments: {}.',
    'task_report': 'Read an owned task report or one child step report. Arguments: task_id, optional step_id, optional pointer (JSON pointer), offset (default 0), limit (1..100, default 20). Returns child step IDs and discoverable artifacts.',
    'read_artifact': 'Read a discovered JSON or CSV result by artifact_id; optional pointer, offset, limit. No arbitrary paths. Large arrays are paginated.',
    'library_status': 'Inspect configured search.batch catalogs and registry availability without a full-library scan. Arguments: {}.',
    'molecule_lookup': 'Look up an exact molecule ID or exact original source_name in the configured library. Arguments: molecule_id OR source_name; optional limit 1..100. Returns conformer identities and original source provenance, not a new search.',
    'literature_search': 'Search Europe PMC titles and abstracts. Arguments: query (English search expression), optional limit 1..10. Includes publication IDs, URLs and abstract evidence. Not full-text review.',
    'viewer_status': 'Read current bridge state, catalog and latest execution receipt. Arguments: {}.',
    'viewer_read': 'Read a discovered JSON/CSV viewer result by artifact_id, optional pointer, offset, limit.',
    'workflow': 'Execute ONE existing validated Chat decision using the workflow contract. Arguments: decision (the same intent/message/task_id/request and permitted extra fields as the supplied contract). Reuse task IDs. No arbitrary Python or shell. Queued work must not be launched again. PyMOL scene and program operations use the live bridge.',
}


def bounded(value, depth=0):
    if depth>8:return {'truncated': True, 'reason': 'depth; read a narrower JSON pointer'}
    if isinstance(value,dict):
        items=list(value.items());result={k:bounded(v,depth+1) for k,v in items[:60]}
        if len(items)>60:result['_truncated_keys']=len(items)-60
        return result
    if isinstance(value,list):
        return dict(items=[bounded(v,depth+1) for v in value[:20]],total=len(value),truncated=len(value)>20)
    if isinstance(value,str) and len(value)>5000:return value[:5000]+' [truncated]'
    return value


def page(value,args):
    pointer=args.get('pointer','')
    if not isinstance(pointer,str) or len(pointer)>1000 or (pointer and not pointer.startswith('/')):raise ValueError('Use a JSON pointer beginning with /')
    for token in pointer.split('/')[1:]:
        token=token.replace('~1','/').replace('~0','~')
        if isinstance(value,list):
            if not token.isdigit():raise ValueError('Array pointer requires a nonnegative index')
            value=value[int(token)]
        else:value=value[token]
    offset=args.get('offset',0);limit=args.get('limit',20)
    if type(offset) is not int or not 0<=offset<=100000:raise ValueError('offset must be 0..100000')
    if type(limit) is not int or not 1<=limit<=100:raise ValueError('limit must be 1..100')
    if isinstance(value,list):return dict(items=[bounded(v) for v in value[offset:offset+limit]],total=len(value),offset=offset,truncated=offset+limit<len(value))
    if isinstance(value,dict):
        items=list(value.items());return dict(value={k:bounded(v) for k,v in items[offset:offset+limit]},total_keys=len(items),offset=offset,truncated=offset+limit<len(items))
    return dict(value=bounded(value))


class DomainTools:
    def __init__(self,app,sid,request,provider,fetch=None):
        from .prompt_workflow import project_root
        self.app=app;self.sid=sid;self.request=request;self.provider=provider
        ctx=app.context
        self.project=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id']).resolve()
        self.viewer=(app.root/'viewers'/sid).resolve()
        self.artifacts={};self.fetch=fetch or self.fetch_json

    @staticmethod
    def fetch_json(url):
        from .prompt_plan import NoRedirect
        with build_opener(NoRedirect()).open(Request(url,headers={'User-Agent':'AIDD-domain-agent/1.0'}),timeout=20) as stream:
            payload=stream.read(2*1024*1024+1)
        if len(payload)>2*1024*1024:raise ValueError('Literature response exceeds size limit')
        return json.loads(payload)

    def register(self,path,root):
        from .project_context import ensure_within
        path=ensure_within(Path(path),root)
        if not path.is_file() or path.suffix.lower() not in {'.json','.csv'}:return None
        identifier='A'+hashlib.sha256(str(path).encode()).hexdigest()[:16]
        self.artifacts[identifier]=(path,root)
        return dict(artifact_id=identifier,name=path.name,path=str(path),bytes=path.stat().st_size)

    def discover(self,value,root):
        found={}
        def visit(v,depth=0):
            if depth>6 or len(found)>=100:return
            if isinstance(v,dict):
                for child in v.values():visit(child,depth+1)
            elif isinstance(v,list):
                for child in v[:100]:visit(child,depth+1)
            elif isinstance(v,str) and len(v)<2000 and Path(v).is_absolute() and Path(v).suffix.lower() in {'.json','.csv'}:
                try:
                    row=self.register(v,root)
                    if row:found[row['artifact_id']]=row
                except (ValueError,OSError):pass
        visit(value)
        return list(found.values())

    @staticmethod
    def json_file(path):
        if path.stat().st_size>20*1024*1024:raise ValueError('JSON exceeds 20 MiB read limit; use a smaller report')
        return json.loads(path.read_text(encoding='utf-8'))

    def call(self,name,args):
        from .chat_agent import read_json,validate_route
        from .project_context import ensure_within
        if not isinstance(args,dict):raise ValueError('Tool arguments must be an object')
        allowed={
            'structure_workflow':{'task_id'},
            'tasks':set(),'task_report':{'task_id','step_id','pointer','offset','limit'},
            'read_artifact':{'artifact_id','pointer','offset','limit'},'viewer_read':{'artifact_id','pointer','offset','limit'},
            'library_status':set(),'molecule_lookup':{'molecule_id','source_name','limit'},
            'literature_search':{'query','limit'},'viewer_status':set(),'workflow':{'decision'}}
        if name not in allowed or set(args)-allowed[name]:raise ValueError('Unknown tool or arguments')
        if name=='structure_workflow':
            from .structure_workflow_context import inspect
            return inspect(self,args)
        if name=='tasks':
            jobs=self.app.jobs(self.sid)
            return dict(tasks=[{k:j[k] for k in ('id','status','request','report','error')} for j in jobs[-100:]],total=len(jobs),truncated=len(jobs)>100)
        if name=='task_report':
            job=self.app.task(self.sid,args['task_id'])
            if not job.get('report') or not Path(job['report']).is_file():return dict(status=job['status'],message='No report exists yet')
            path=ensure_within(Path(job['report']),self.project);data=self.json_file(path)
            if args.get('step_id'):
                row=data.get('steps',{})[args['step_id']]
                if row.get('result',{}).get('report'):
                    path=ensure_within(Path(row['result']['report']),self.project);data=self.json_file(path)
                else:data=row
            return dict(task_id=job['id'],task_status=job['status'],source=str(path),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),content=page(data,args),
                artifacts=self.discover(data,self.project),step_ids=list(data.get('steps',{})) if isinstance(data,dict) else [])
        if name in {'read_artifact','viewer_read'}:
            path,root=self.artifacts[args['artifact_id']];path=ensure_within(path,root)
            if path.suffix=='.json':content=page(self.json_file(path),args)
            else:
                page([],args)
                if path.stat().st_size>20*1024*1024:raise ValueError('CSV exceeds 20 MiB read limit')
                with path.open(encoding='utf-8-sig',newline='') as handle:
                    import itertools
                    offset=args.get('offset',0);limit=args.get('limit',20)
                    rows=list(itertools.islice(csv.DictReader(handle),offset,offset+limit+1))
                content=dict(items=rows[:limit],offset=offset,has_more=len(rows)>limit)
            return dict(source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),content=content)
        if name in {'library_status','molecule_lookup'}:
            runtime=read_json(self.app.runtime) or {};batch=runtime.get('search',{}).get('batch')
            if not batch:return dict(status='not_configured',message='Configure runtime.search.batch to connect the precomputed molecule library')
            batch=Path(batch).resolve();catalog_path=batch/'artifacts/catalog.json'
            catalog=self.json_file(catalog_path)
            if name=='library_status':
                return dict(status='connected',batch=str(batch),library_id=catalog.get('library_id'),
                    conformers=catalog.get('conformers'),shards=len(catalog.get('shards',[])),
                    registry_available=(batch/'registry.sqlite3').is_file(),chemical_catalog_available=(batch/'chemical/catalog.json').is_file(),
                    catalog_sha256=hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
                    scope='Catalog metadata only, not full-library byte integrity or a scientific search')
            keys=set(args)&{'molecule_id','source_name'}
            if len(keys)!=1:raise ValueError('Supply exactly one molecule_id or exact source_name')
            key=next(iter(keys));value=args[key];limit=args.get('limit',20)
            if not isinstance(value,str) or not 1<=len(value)<=300 or type(limit) is not int or not 1<=limit<=100:raise ValueError('Invalid lookup value/limit')
            column='id' if key=='molecule_id' else 'source_name'
            deadline=time.monotonic()+3
            with sqlite3.connect((batch/'registry.sqlite3').as_uri()+'?mode=ro',uri=True) as db:
                db.row_factory=sqlite3.Row;db.set_progress_handler(lambda:int(time.monotonic()>deadline),10000)
                rows=db.execute(f'SELECT id,source_name,library_id FROM molecule WHERE library_id=? AND {column}=? LIMIT ?',(catalog['library_id'],value,limit+1)).fetchall()
                molecules=[]
                for r in rows[:limit]:
                    conformers=db.execute('SELECT id,conformer_index,source_record_name,source_path,source_record_index,content_sha256,atom_count,bond_count FROM conformer WHERE molecule_id=? ORDER BY conformer_index LIMIT 101',(r['id'],)).fetchall()
                    molecules.append(dict(r)|dict(conformers=[dict(c) for c in conformers[:100]],conformers_truncated=len(conformers)>100))
            return dict(molecules=molecules,truncated=len(rows)>limit,source=str(batch/'registry.sqlite3'))
        if name=='literature_search':
            query=args.get('query');limit=args.get('limit',5)
            if not isinstance(query,str) or not 1<=len(query)<=500 or type(limit) is not int or not 1<=limit<=10:raise ValueError('Invalid literature query/limit')
            url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode(dict(query=query,format='json',resultType='core',pageSize=limit))
            data=self.fetch(url);papers=[]
            for row in data.get('resultList',{}).get('result',[])[:limit]:
                source=row.get('source','');pid=str(row.get('id',''))
                if not re.fullmatch(r'[A-Za-z0-9_-]+',source) or not re.fullmatch(r'[A-Za-z0-9_.-]+',pid):continue
                papers.append(dict(id=source+':'+pid,title=row.get('title'),year=row.get('pubYear'),doi=row.get('doi'),
                    url=f'https://europepmc.org/article/{source}/{pid}',
                    abstract=re.sub('<[^>]+>',' ',row.get('abstractText',''))[:6000]))
            return dict(query=query,url=url,papers=papers,hit_count=data.get('hitCount'),retrieved_at=time.time(),scope='Title/abstract search evidence, not full-text review')
        if name=='viewer_status':
            from .pymol_bridge import bridge_status
            value=bridge_status(self.viewer);saved=read_json(self.viewer/'catalog.json')
            if saved:
                self.app.task(self.sid,saved['task_id'])
                value['catalog']=[{k:r[k] for k in ('id','label','kind')} for r in saved['structures']]
            return dict(state=bounded(value),artifacts=self.discover(value,self.viewer))
        if name=='workflow':
            decision=args['decision'];validate_route(decision)
            before={j['id'] for j in self.app.jobs(self.sid)}
            reply=self.app.execute_decision(self.sid,self.request,self.provider,decision,record=False)
            try:result=json.loads(reply)
            except ValueError:result=dict(message=reply)
            new=[j for j in self.app.jobs(self.sid) if j['id'] not in before]
            if decision['intent']=='resume':
                new=[self.app.task(self.sid,decision['task_id'])]
            if new:return dict(status='queued',tasks=[{k:j[k] for k in ('id','status','request')} for j in new],message=reply)
            # Viewer dispatch is asynchronous. Resolve its receipt once, never resubmit.
            if result.get('status')=='queued' and result.get('operation_id') and result.get('result'):
                from .pymol_agent import wait_result
                result=wait_result(result,seconds=5)
            return dict(result=result,artifacts=self.discover(result,self.viewer if decision['intent'] in {'pymol','pymol_agent'} else self.project))
