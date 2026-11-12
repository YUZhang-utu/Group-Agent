"""Bounded local coordination for explicit read-only saved-review summaries."""
import re


def matches(request):
    text=request.casefold()
    saved=re.search(r'\b(saved|previous|earlier|prior|existing)\b|\u5df2\u4fdd\u5b58|\u4e4b\u524d|\u5df2\u6709',text)
    reviews=re.search(r'\breviews?\b|\bliterature\b|\u6587\u732e|\u5ba1\u9605',text)
    summary=re.search(r'\bsummari[sz]e\b|\bsummary\b|\blearned\b|\brecall\b|\bconclude\b|\u603b\u7ed3|\u56de\u987e|\u7ed3\u8bba',text)
    if not (saved and reviews and summary):return False
    # Mixed work stays with the normal planner. Negated compute instructions do
    # not turn an explicit read-only summary into a compute request.
    affirmative=re.sub(r"\b(?:do not|don't|never)\b[^.\n;]*",'',text)
    affirmative=re.sub(r'(?:\u4e0d\u8981|\u4e0d\u8fd0\u884c|\u7981\u6b62)[^\u3002\n;]*','',affirmative)
    if re.search(r'\b(run|launch|predict|dock|screen|simulate|execute|submit|delete|archive|remove|change)\b|\u8fd0\u884c|\u542f\u52a8|\u9884\u6d4b|\u5bf9\u63a5|\u5220\u9664',affirmative):return False
    if re.search(r'\b(new|latest)\s+(papers|studies|literature)\b',affirmative):return False
    return True


class SavedReviewFlow:
    def __init__(self,request):
        ids=set(re.findall(r'\b[a-f0-9]{32}\b',request))
        self.scope=dict(document_id=next(iter(ids))) if len(ids)==1 else dict(query=request[:500])
        self.queue=[dict(operation='list',limit=3,**self.scope)]
        self.selected=[];self.read={};self.partial=set();self.failed=set()
        direct=re.search(r'\bsaved review(?:\s+id)?\s*[:#]?\s*([a-f0-9]{32})\b',request,re.I)
        if direct:
            self.selected=[direct[1]];self.queue=[dict(operation='read',review_id=direct[1],offset=0)]
        self.allow_index=not re.search(r"\b(?:do not|don't|never)\s+index\b|\u4e0d\u8981\u7d22\u5f15|\u4e0d\u91cd\u65b0\u7d22\u5f15",request,re.I)
        self.index_batches=0;self.index_offset=0;self.index_more=True
        self.list_more=False;self.characters=0;self.reason=None;self.recent_fallback=False

    def next_call(self):
        if self.characters>=20000:
            self.reason='summary_character_budget';self.queue=[]
        return self.queue.pop(0) if self.queue else None

    def accept(self,args,result,status):
        op=args['operation']
        if status!='complete':
            if op=='read':self.failed.add(args['review_id'])
            else:self.reason='lookup_failed';self.queue=[]
            return
        if op=='list':
            found=result.get('reviews',[])[:3]
            if found:
                self.selected=[r['review_id'] for r in found];self.list_more=result.get('has_more',False)
                self.queue.extend(dict(operation='read',review_id=identifier,offset=0) for identifier in self.selected)
            elif 'query' in self.scope:
                self.scope={};self.recent_fallback=True;self.queue.append(dict(operation='list',limit=3))
            elif self.allow_index and self.index_batches<2 and self.index_more:
                self.queue.append(dict(operation='index',offset=self.index_offset,limit=10))
            else:self.reason='no_available_saved_reviews_in_checked_scope'
        elif op=='index':
            self.index_batches+=1;self.index_more=result.get('has_more',False)
            self.index_offset=result.get('next_offset') or 0
            self.queue.append(dict(operation='list',limit=3,**self.scope))
        elif op=='read':
            identifier=args['review_id'];text=result.get('summary','')
            self.characters+=len(text);self.read[identifier]=self.read.get(identifier,0)+len(text)
            if result.get('has_more'):
                self.partial.add(identifier);offset=result.get('next_offset')
                if type(offset) is int and offset>args.get('offset',0) and text and self.read[identifier]<10000:
                    self.queue.append(dict(operation='read',review_id=identifier,offset=offset))
            else:self.partial.discard(identifier)

    def report(self):
        return dict(mode='local_saved_review_coordination_no_json_planning',selected_review_ids=self.selected,
            read_review_ids=sorted(k for k,v in self.read.items() if v),summary_characters_read=self.characters,
            partial_review_ids=sorted(self.partial),unread_review_ids=sorted(set(self.selected)-set(self.read)),
            failed_review_ids=sorted(self.failed),more_review_records_available=self.list_more,
            old_trace_batches_checked=self.index_batches,older_traces_unchecked=self.index_batches>=2 and self.index_more,
            indexing_allowed=self.allow_index,
            recent_project_fallback=self.recent_fallback,stop_reason=self.reason,
            pending_local_operations=len(self.queue),
            interpretation='At most three available project review records; model interpretations, not fresh paper reading. '
            'Explicitly disclose unread/partial records and any fallback to recent project records instead of query matches.')
