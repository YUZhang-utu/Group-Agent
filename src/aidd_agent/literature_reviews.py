"""Project-local, source-checked model review records; never experimental evidence."""
import json
import re
import time

from . import project_knowledge as knowledge
from .literature_citations import label


def initialize(db):
    db.execute('CREATE TABLE IF NOT EXISTS literature_reviews(id TEXT PRIMARY KEY,created REAL,payload TEXT,sha256 TEXT)')


def checked_references(db,directory,refs,verify_source=True):
    checked=[]
    for ref in refs:
        row=db.execute('SELECT d.* FROM documents d JOIN document_state s ON s.id=d.id WHERE d.id=? AND s.archived=0',
                       (ref['document_id'],)).fetchone()
        if not row or row['sha256']!=ref['sha256']:raise ValueError('Saved review source unavailable')
        if verify_source:knowledge.verify_document(db,directory,row)
        else:knowledge.verify_metadata(db,row)
        checked.append(dict(ref,title=row['title'],source=row['source'],bibliography=knowledge.bibliography(db,row['id'])))
    return checked


def save(app,sid,audit,summary,refs):
    if audit.get('status')!='complete' or not audit.get('evidence_ids') or not refs or len(refs)>12:return None
    records=[r for r in audit['steps'] if r.get('step',{}).get('kind')=='tool']
    if any(r['step']['tool'] not in {'project_knowledge','literature_reviews'} or
           r['step'].get('arguments',{}).get('collection','public')!='public' for r in records):return None
    fresh=any(r.get('status')=='complete' and r['step']['tool']=='project_knowledge' and
              any('text' in item for item in r.get('result',{}).get('results',[])) for r in records)
    if not fresh:return None  # Do not recursively save summaries of summaries.
    text=summary.split('\n\nSources:',1)[0]
    if not 1<=len(text)<=24000:return None
    directory=knowledge.folder(app,sid,'public')
    with app.lock,knowledge.connect(directory) as db:
        initialize(db)
        checked_references(db,directory,refs)
        if len(json.dumps(refs,ensure_ascii=False))>4500:return None
        # A literal reference to another document is a dependency even when it
        # appears only in a comparison with an earlier same-chat review.
        known={r['document_id'] for r in refs}
        for identifier in set(re.findall(r'\b[a-f0-9]{32}\b',audit.get('answer','')))-known:
            row=db.execute('SELECT * FROM documents WHERE id=?',(identifier,)).fetchone()
            if not row:return None
            refs.append(dict(document_id=identifier,sha256=row['sha256'],title=row['title'],source=row['source'],pages=[],fresh_pages=[]))
        checked_references(db,directory,refs)
        identifier=audit['id']
        payload=dict(id=identifier,created=time.time(),origin_session=sid,origin_trace=identifier,
            summary=text,references=refs,coverage=audit.get('indexed_reading_coverage',{}),
            verification='Model-generated literature interpretation; not independently verified or experimental evidence')
        inserted=db.execute('INSERT OR IGNORE INTO literature_reviews VALUES(?,?,?,?)',
                   (identifier,payload['created'],json.dumps(payload,ensure_ascii=False),knowledge.digest(payload))).rowcount
    return dict(review_id=identifier,status='saved' if inserted else 'already_saved',scope='Owned project across conversations',verification=payload['verification'])


def handle(app,sid,args):
    if not isinstance(args,dict) or set(args)-{'operation','query','review_id','document_id','offset','limit'}:raise ValueError('Unsupported literature review arguments')
    op=args.get('operation','list');offset=args.get('offset',0);limit=args.get('limit',5)
    if op not in {'list','read','index'} or type(offset) is not int or offset<0 or type(limit) is not int or not 1<=limit<=10:raise ValueError('Use list/read/index, nonnegative offset and limit 1..10')
    if op=='index':return index_previous(app,sid,offset,limit)
    query=args.get('query','')
    if not isinstance(query,str) or len(query)>500:raise ValueError('Use a review query of at most 500 characters')
    for key in ('review_id','document_id'):
        if args.get(key) is not None and (not isinstance(args[key],str) or not re.fullmatch('[a-f0-9]{32}',args[key])):raise ValueError('Invalid review/document ID')
    if op=='read' and not args.get('review_id'):raise ValueError('Supply review_id from the project review index')
    result=dict(collection='public',operation=op,retrieval='saved_model_review_not_original_paper',reviews=[],has_more=False,next_offset=None,
        interpretation='Saved model interpretation, not fresh paper reading or verified experimental findings')
    directory=knowledge.folder(app,sid,'public')
    if not (directory/'index.sqlite').exists():
        if op=='read':raise ValueError('Saved review unavailable')
        return result
    with app.lock,knowledge.connect(directory) as db:
        initialize(db)
        sql='SELECT * FROM literature_reviews'+(' WHERE id=?' if op=='read' else '')+' ORDER BY created DESC,id'
        rows=db.execute(sql,(args['review_id'],) if op=='read' else ())
        matches=[];words=set(knowledge.terms(query))
        for row in rows:
            try:
                payload=json.loads(row['payload'])
                if knowledge.digest(payload)!=row['sha256']:continue
                if args.get('document_id') and args['document_id'] not in {r['document_id'] for r in payload['references']}:continue
                refs=checked_references(db,directory,payload['references'],verify_source=op=='read')
            except (ValueError,KeyError,OSError):continue
            if op=='read':
                text=payload['summary'];end=min(len(text),offset+5000)
                if offset>len(text):raise ValueError('Review offset exceeds saved summary')
                result.update(review_id=payload['id'],summary=text[offset:end],references=refs,
                    created=payload['created'],coverage=payload['coverage'],offset=offset,total_characters=len(text),
                    has_more=end<len(text),next_offset=end if end<len(text) else None,verification=payload['verification'])
                if len(json.dumps(result,ensure_ascii=False))>11000:raise ValueError('Saved review metadata exceeds the model reading budget')
                return result
            score=len(words & set(knowledge.terms(payload['summary']+' '+ ' '.join(label(r)+' '+r['document_id'] for r in refs))))
            if words and not score:continue
            matches.append((score,dict(review_id=payload['id'],created=payload['created'],
                papers=[dict(document_id=r['document_id'],title=label(r),source=r['source']) for r in refs],
                summary_characters=len(payload['summary']),verification=payload['verification']),refs))
        if op=='read':raise ValueError('Saved review unavailable')
        matches.sort(key=lambda item:(item[0],item[1]['created']),reverse=True)
        eligible=[]
        for item in matches:
            try:checked_references(db,directory,item[2])
            except (ValueError,KeyError,OSError):continue
            eligible.append(item)
            if len(eligible)>offset+limit:break
        selected=eligible[offset:offset+limit]
        result.update(reviews=[item[1] for item in selected],offset=offset,
                      has_more=len(eligible)>offset+limit,next_offset=offset+limit if len(eligible)>offset+limit else None)
        while len(json.dumps(result,ensure_ascii=False))>8500 and len(result['reviews'])>1:result['reviews'].pop()
        if len(result['reviews'])<len(selected):result.update(has_more=True,next_offset=offset+len(result['reviews']))
        if len(json.dumps(result,ensure_ascii=False))>11000:raise ValueError('Review index entry exceeds the model budget')
    return result


def brief(app,sid,request):
    found=handle(app,sid,dict(operation='list',query=request[:500],limit=3))
    if not found['reviews']:found=handle(app,sid,dict(operation='list',limit=3))
    return dict(reviews=found['reviews'],instruction='Use literature_reviews read for relevant saved interpretations before rereading papers. '
        'These records are project-scoped and source-checked, not independently verified science. '
        'Read original passages for new claims or verification. Newer summaries do not automatically resolve contradictions.')


def index_previous(app,sid,offset,limit):
    """Adopt eligible old owned-project reviews without rereading a PDF or a model call."""
    from .project_context import ensure_within
    from .literature_citations import present
    from .literature_reading import coverage
    pid=app.context_for(sid)['project_id']
    with app.connect() as db:sessions=[r['id'] for r in db.execute('SELECT id FROM sessions WHERE project_id=?',(pid,))]
    paths=[]
    for session in sessions:
        if re.fullmatch('[a-f0-9]{32}',session):paths.extend((app.root/'agents'/session).glob('*.json'))
    paths.sort(key=lambda p:p.stat().st_mtime,reverse=True)
    saved=0
    for path in paths[offset:offset+limit]:
        try:
            source_session=path.parent.name
            path=ensure_within(path,app.root/'agents')
            if path.parent!=app.root/'agents'/source_session:continue
            if path.stat().st_size>2*1024*1024:continue
            audit=json.loads(path.read_text(encoding='utf-8'))
            if audit.get('id')!=path.stem or audit.get('status')!='complete' or not isinstance(audit.get('answer'),str):continue
            audit['indexed_reading_coverage']=coverage(audit.get('steps',[]))
            text,refs,_=present(audit['answer'],audit['steps'],audit['indexed_reading_coverage'])
            receipt=save(app,path.parent.name,audit,text,refs)
            if receipt and receipt['status']=='saved':saved+=1
        except (ValueError,KeyError,OSError):continue
    return dict(collection='public',operation='index',indexed_reviews=saved,offset=offset,
        has_more=len(paths)>offset+limit,next_offset=offset+limit if len(paths)>offset+limit else None,
        interpretation='Eligible previous public reviews indexed with source checks; no model call or scientific calculation')
