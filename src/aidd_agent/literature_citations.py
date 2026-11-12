"""Paper-level presentation; internal evidence IDs remain in the audit trace."""
import re
from urllib.parse import quote


def label(document):
    meta=document.get('bibliography',{})
    title=meta.get('display_title') or document.get('title') or 'Public paper'
    title=re.sub(r'\.(pdf|txt|md)$','',title,flags=re.I).replace('_',' ')
    title=' '.join(title.split())
    if len(title)>80:title=title[:77]+'...'
    prefix=', '.join(v for v in (meta.get('authors'),meta.get('year')) if v)
    return (prefix+' - ' if prefix else '')+title


def source_link(document):
    title=label(document).translate(str.maketrans({'[':'(',']':')','`':"'",'\n':' ','\r':' '}))
    source=document.get('source','')
    return f'[{title}]({quote(source,safe=":/?#=&%+-._~")})' if source.startswith(('https://','http://')) else title


def references(steps):
    documents={};evidence={}
    for record in steps:
        if record.get('status')!='complete':continue
        result=record.get('result',{});tool=record.get('step',{}).get('tool')
        if result.get('collection')!='public':continue
        if tool=='project_knowledge':rows=result.get('results',[])
        elif tool=='literature_reviews':rows=result.get('references',[])
        else:continue
        for row in rows:
            identifier=row.get('id') or row.get('document_id')
            if not identifier:continue
            doc=documents.setdefault(identifier,dict(document_id=identifier,title=row.get('title','Public paper'),
                source=row.get('source',''),sha256=row.get('sha256'),bibliography=row.get('bibliography',{}),pages=set(),fresh_pages=set()))
            if type(row.get('page')) is int:
                doc['pages'].add(row['page']);doc['fresh_pages'].add(row['page'])
            doc['pages'].update(row.get('pages',[]))
            if record.get('evidence_id'):evidence.setdefault(record['evidence_id'],set()).add(identifier)
    for doc in documents.values():
        doc['pages']=sorted(doc['pages']);doc['fresh_pages']=sorted(doc['fresh_pages'])
    # Distinguish identical titles without using transient tool-call numbers.
    labels={}
    for doc in documents.values():labels.setdefault(label(doc),[]).append(doc)
    for same in labels.values():
        if len(same)>1:
            for doc in same:doc['bibliography']=dict(doc['bibliography'],display_title=(doc['bibliography'].get('display_title') or doc['title'])[:60]+' / '+doc['document_id'][:6])
    return documents,evidence


def present(answer,steps,state,replace_document_ids=True):
    docs,ids=references(steps)
    if not docs:return answer,[],False
    def replace(match):
        for identifier,doc in docs.items():
            if re.fullmatch(re.escape(identifier)+r':page-\d+(?::part-\d+)?',match[1].strip()):return source_link(doc)
        if not re.fullmatch(r'\s*E\d+(?:\s*[,;]\s*E\d+)*(?:\s*,\s*(?:pp?\.?|page)\s*[\d, -]+)?\s*',match[1]):return match[0]
        cited=re.findall(r'\bE\d+\b',match[1]);selected=set()
        if not cited or any(key not in ids for key in cited):return match[0]
        for key in cited:selected.update(ids[key])
        return ' '.join(source_link(docs[key]) for key in sorted(selected))
    text=re.sub(r'\[([^\[\]\n]+)\]',replace,answer)
    if replace_document_ids:
        for identifier,doc in docs.items():
            text=re.sub(re.escape(identifier)+r'(?::page-\d+(?::part-\d+)?)?',lambda _:label(doc),text)
    # Repeated receipts for one paper should not create adjacent identical links.
    for doc in docs.values():
        link=source_link(doc)
        text=re.sub('(?:'+re.escape(link)+r'\s*){2,}',lambda _:link+' ',text)
    details=[]
    for identifier,doc in docs.items():
        pages=[p for p in state.get('pages',[]) if p['document_id']==identifier]
        complete=sorted(p['page'] for p in pages if p['indexed_page_complete'])
        partial=sorted(p['page'] for p in pages if not p['indexed_page_complete'])
        scope=f'This turn: complete page reads {complete}; fragment/partial pages {partial}.' if pages else 'Saved review/source metadata; original paper text was not reread this turn.'
        details.append(source_link(doc)+'\n'+scope+' Recorded source pages: '+str(doc['pages'])+'.')
    text+='\n\nSources: '+'; '.join(source_link(doc) for doc in docs.values())
    text+='\n\n```sources\n'+'\n\n'.join(details)+'\nIndexed text only; section/full-paper completion is not established. Figures and scientific accuracy require review.\n```'
    return text,list(docs.values()),True
