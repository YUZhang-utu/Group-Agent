"""Bounded same-chat comparison context and explicit-page reading receipts."""
import re


def public_records(steps):
    for record in steps:
        step=record.get('step',{})
        if step.get('tool')=='project_knowledge' and record.get('status')=='complete' and record.get('result',{}).get('collection')=='public':
            yield step.get('arguments',{}),record['result']


def coverage(steps):
    pages={}
    for args,result in public_records(steps):
        rows=result.get('results',[])
        for row in rows:
            if 'text' not in row:continue
            key=(row['id'],row['page'])
            page=pages.setdefault(key,dict(document_id=key[0],page=key[1],parts=set(),read_ranges=[],end=None))
            page['parts'].add(row['part'])
        if args.get('operation')!='read' or args.get('page') is None or not rows:continue
        key=(args.get('document_id'),args['page'])
        if key not in pages or any((r['id'],r['page'])!=key for r in rows):continue
        page=pages[key];start=args.get('offset',0);end=start+len(rows)
        page['read_ranges'].append((start,end))
        if not result.get('has_more'):page['end']=end
    summaries=[]
    for page in pages.values():
        reached=0
        for start,end in sorted(page['read_ranges']):
            if start>reached:break
            reached=max(reached,end)
        complete=page['end'] is not None and reached>=page['end']
        summaries.append(dict(document_id=page['document_id'],page=page['page'],passages_seen=len(page['parts']),
            indexed_page_complete=complete,read_from_offset_zero_through=reached,
            mode='complete_page_read' if complete else 'partial_page_read' if page['read_ranges'] else 'retrieved_fragments'))
    return dict(pages=summaries,scope='Current turn only; indexed text coverage, not PDF figures/tables or section completeness')


def next_page_call(steps):
    """Continue only the last selected page; never retry an error or expand a search."""
    records=[r for r in steps if r.get('step',{}).get('tool')=='project_knowledge']
    if not records:return None
    record=records[-1];args=record['step'].get('arguments',{});result=record.get('result',{})
    if record.get('status')!='complete' or result.get('collection')!='public':return None
    if args.get('operation')!='read' or args.get('page') is None or not result.get('results') or not result.get('has_more'):return None
    offset=result.get('next_offset')
    if type(offset) is not int or offset<=args.get('offset',0):return None
    candidate=dict(args,offset=offset)
    if any(r['step'].get('arguments')==candidate for r in records):return None
    return candidate


def previous_summary(history,steps,maximum=18000):
    docs={r['id'] for _,result in public_records(steps) for r in result.get('results',[]) if 'text' in r}
    sources={r['source'] for _,result in public_records(steps) for r in result.get('results',[]) if 'text' in r and r.get('source')}
    for row in reversed(history):
        if row.get('role')!='assistant':continue
        text=row.get('text','')
        if not any(doc in text for doc in docs|sources) or text.lstrip().startswith(('{','[')):continue
        text=text.split('\n\nEvidence sources:',1)[0].split('\n\nAgent trace:',1)[0]
        text=text.split('\n\nSources:',1)[0]
        if len(text)>maximum:
            return dict(status='omitted_over_budget',characters=len(text),reason='Previous answer is not silently truncated')
        return dict(status='available',text=text,
            instruction='Previous assistant answer from this conversation, for comparison only; not source evidence. '
            'Its E-number citations belong to a prior turn. Verify corrections against current retrieved passages; do not reuse old IDs.')
    return dict(status='not_found_in_recent_same_conversation',instruction='No matching previous answer in the bounded recent history; do not invent one')


def methods_feedback(request,steps):
    # This is a reading nudge for an explicit Methods review, not a claim that a
    # heading or experiment has been semantically identified by the software.
    if not re.search(r'\bmethods?\b|\bmethodology\b|\u65b9\u6cd5',request,re.I):return None
    if not re.search(r'\bread\b|\breview\b|\bsummar|\u9605\u8bfb|\u603b\u7ed3',request,re.I):return None
    pages=coverage(steps)['pages']
    fragments=[dict(document_id=p['document_id'],page=p['page']) for p in pages if p['page']!=1 and p['mode']=='retrieved_fragments']
    if not fragments or any(p['page']!=1 and p['indexed_page_complete'] for p in pages):return None
    return ('The user requested Methods reading, but only retrieved fragments are available outside page 1. '
        'Use the search text to select the relevant Methods page(s), then call project_knowledge read with '
        'document_id and page. Selected pages paginate automatically. Do not read every hit blindly. '
        'If the relevant section cannot be located, state that explicitly. Candidate pages: '+str(fragments[:12]))


def coverage_note(state):
    groups={}
    for page in state['pages']:groups.setdefault(page['document_id'],[]).append(page)
    if not groups:return ''
    lines=['Indexed reading coverage (this turn; section/full-paper completion is not established):']
    for doc,pages in list(groups.items())[:8]:
        complete=sorted(p['page'] for p in pages if p['indexed_page_complete'])
        partial=sorted(p['page'] for p in pages if not p['indexed_page_complete'])
        lines.append(f'- {doc}: complete page reads {complete}; fragment/partial pages {partial}.')
    return '\n\n'+'\n'.join(lines)
