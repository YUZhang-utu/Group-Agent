const $=id=>document.getElementById(id);
const fragment=new URLSearchParams(location.hash.slice(1));
if(fragment.get('token')){sessionStorage.setItem('aidd-token',fragment.get('token'));history.replaceState(null,'',location.pathname);}
if(fragment.get('session'))sessionStorage.setItem('aidd-session',fragment.get('session'));
let session=sessionStorage.getItem('aidd-session'),last='',polling=null, busy=false,lastState='';
let attachedFiles=[];
let projectId=null;
function markdown(value){const body=node('div',undefined,'markdown');body.innerHTML=AiddMarkdown.render(value);return body;}
function setBusy(value){busy=value;for(const id of ['send','new','provider','upload-files','project-select','project-create'])if($(id))$(id).disabled=value;for(const b of $('sessions').querySelectorAll('button'))b.disabled=value;}
async function api(path,body){const r=await fetch(path,{method:body?'POST':'GET',headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||''),'Content-Type':'application/json'},body:body?JSON.stringify(body):undefined});const data=await r.json();if(!r.ok)throw Error(data.error||'Request failed');return data;}
function node(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;}
function selectSession(id){session=id;attachedFiles=[];renderAttachments();sessionStorage.setItem('aidd-session',session);last='';lastState='';}
async function newSession(){const r=await api('/api/session',projectId?{project_id:projectId}:{});selectSession(r.session);await refresh();}
async function send(text){
if(busy)return;setBusy(true);$('error').textContent='';
let finished=false,answer='',submitted=false;
const originalAttachments=[...attachedFiles];
live.hidden=false;live.replaceChildren(node('span','MEDCHEM Agent','role'));const status=node('div','Connecting...','stream-status');const preview=markdown('');live.append(status,preview);
try{
  if(!session)await newSession();
  submitted=true;
  const response=await fetch('/api/message/stream',{method:'POST',headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||''),'Content-Type':'application/json'},body:JSON.stringify({session,text,provider:$('provider').value,attachments:originalAttachments.map(file=>file.id)})});
  if(!response.ok){const error=await response.json();throw Error(error.error||'Request failed');}
  if($('prompt').value===text)$('prompt').value='';
  attachedFiles=[];renderAttachments();
  const reader=response.body.getReader(),decoder=new TextDecoder();let buffer='';
  function event(row){
    if(row.type==='answer_delta'){answer+=row.text;preview.innerHTML=AiddMarkdown.render(answer);status.textContent='Writing response...';}
    else if(row.type==='answer_reset'){answer='';preview.replaceChildren();}
    else if(row.type==='project_opened'){projectId=row.project_id;selectSession(row.session);}
    else if(row.type==='done'){finished=true;answer=row.text;preview.innerHTML=AiddMarkdown.render(answer);status.textContent='Response complete';}
    else if(row.type==='error')throw Error(row.text);
    else if(row.text)status.textContent=row.text;
  }
  while(true){const {value,done}=await reader.read();buffer+=decoder.decode(value||new Uint8Array(),{stream:!done});let end;while((end=buffer.indexOf('\n'))>=0){const line=buffer.slice(0,end);buffer=buffer.slice(end+1);if(line.trim())event(JSON.parse(line));}if(done)break;}
  if(buffer.trim())event(JSON.parse(buffer));
  if(!finished)throw Error('Connection ended before the final response. Check saved tasks before resubmitting.');
}catch(e){$('error').textContent=e.message+(submitted?' Existing tasks and uploaded files remain on this host.':'');}
finally{setBusy(false);await refresh(true);live.hidden=true;}
}
function render(state){modeNotice.hidden=state.interface_mode!=='isolated_test';$('task-count').textContent=state.tasks.length;$('task-count').hidden=!state.tasks.length;$('compute').textContent=state.compute_enabled?'Compute enabled':'Preparation only';$('sessions').replaceChildren();for(const s of state.sessions){const b=node('button',s.title,s.id===session?'selected':'');b.disabled=busy;b.onclick=()=>{session=s.id;attachedFiles=[];renderAttachments();sessionStorage.setItem('aidd-session',session);last='';refresh();};$('sessions').append(b);}
const key=JSON.stringify(state.messages);if(key!==last){$('messages').replaceChildren();for(const m of state.messages){const box=node('div',undefined,'message '+m.role);box.append(node('span',m.role==='user'?'You':'MEDCHEM Agent', 'role'),m.role==='assistant'?markdown(m.text):node('span',m.text));$('messages').append(box);}if(!state.messages.length)$('messages').append(node('p','Ask a research question.','empty'));$('messages').scrollTop=$('messages').scrollHeight;last=key;}
$('coverage').replaceChildren();for(const w of state.workflows){const n=node('div',undefined,'coverage-item');n.append(node('strong',w.name),node('div',w.status),node('div',w.scope));$('coverage').append(n);}
$('tasks').replaceChildren();
if(!state.tasks.length)$('tasks').append(node('p','Your tasks will appear here.','empty'));
for(const t of [...state.tasks].reverse())$('tasks').append(taskCard(t));
}
function taskCard(t){
  const card=node('article',undefined,'task '+t.status),request=t.request||'';
  card.append(node('div',t.status+' / '+t.provider,'task-status'),node('p',request.slice(0,180)+(request.length>180?'...':'')),node('div','Task '+t.id,'path'));
  if(request.length>180){const full=node('details');full.append(node('summary','Full request'),node('pre',request));card.append(full);}
  if(t.status==='queued')card.append(node('p','Waiting for the local worker.'));
  if(t.status==='planning')card.append(node('p','Preparing and validating the workflow plan. Scientific execution has not started.'));
  for(const [id,step] of Object.entries(t.details?.steps||{})){
    card.append(node('p',id+' / '+step.action+': '+step.status));
    for(const key of ['input_json','model','manifest','confidence_analysis','names','ranking','original_mol2_directory','posed_mol2_directory','export_manifest'])if(step.result?.[key]){
      const path=step.result[key];card.append(node('div',key+': '+path,'path'));
      const button=node('button','Copy '+key+' path');button.onclick=()=>navigator.clipboard.writeText(path).catch(e=>$('error').textContent=e.message);card.append(button);
    }
  }
  if(t.error)card.append(node('p',t.error,'task-error'));
  for(const k of ['plan','report','log'])if(t[k]){card.append(node('div',k+': '+t[k],'path'));const b=node('button','Copy '+k+' path');b.onclick=()=>navigator.clipboard.writeText(t[k]).catch(e=>$('error').textContent=e.message);card.append(b);}
  if(t.details){const d=node('details');d.append(node('summary','Full execution report'),node('pre',JSON.stringify(t.details,null,2)));card.append(d);}
  if(t.progress){const d=node('details');d.append(node('summary','Execution log'),node('pre',t.progress));card.append(d);}
  const active=['queued','planning','running'].includes(t.status);
  if(active||['failed','blocked','interrupted','cancelled'].includes(t.status)){
    const action=active?'Cancel':'Resume',b=node('button',action+' task');b.onclick=()=>send('/'+action.toLowerCase()+' '+t.id);card.append(b);
  }
  const results=node('button','Show results');results.onclick=()=>send('/results '+t.id);card.append(results);
  return card;
}
function refresh(force=false){
if(polling)return force?polling.then(()=>refresh(true)):polling;
polling=(async()=>{try{
  const state=await api('/api/state'+(session?'?session='+encodeURIComponent(session):''));
  const key=session+JSON.stringify(state);
  if(key!==lastState){render(state);renderProjects(state);renderDocking(state.docking_results);renderWorkflowResults(state.workflow_results);renderViewer(state.viewer);renderAgent(state.agent_runs);renderChains(state.structure_chains);lastState=key;}
  $('connection').textContent='Connected';
}catch(e){$('connection').textContent='Disconnected';$('error').textContent=e.message;}finally{polling=null;}})();
return polling;
}
$('new').onclick=()=>newSession().catch(e=>$('error').textContent=e.message);$('compose').onsubmit=e=>{e.preventDefault();if($('prompt').value.trim())send($('prompt').value);};$('prompt').onkeydown=e=>{if(e.key==='Enter'&&!e.shiftKey){e.preventDefault();$('compose').requestSubmit();}};
const live=node('div',undefined,'message assistant live-answer');live.hidden=true;live.setAttribute('aria-live','polite');$('messages').after(live);
const modeNotice=node('p','INTERFACE TEST: synthetic responses, no live model or scientific execution. Open the configured Agent server to run real tasks.','fixture-notice');modeNotice.hidden=true;$('messages').before(modeNotice);
const uploadBar=node('div',undefined,'upload-bar'),uploadButton=node('button','Attach files');uploadButton.id='upload-files';uploadButton.type='button';
const fileInput=node('input');fileInput.type='file';fileInput.multiple=true;fileInput.hidden=true;fileInput.accept='.pdb,.cif,.mmcif,.sdf,.mol2,.mol,.smi,.smiles,.csv,.tsv,.json,.fasta,.fa,.txt';
const attachmentList=node('div',undefined,'attachment-list');attachmentList.setAttribute('aria-live','polite');
uploadBar.append(uploadButton,node('span','Drop structures, sequences or tables here. Files stay on the Agent host.','upload-hint'),fileInput);$('prompt').before(uploadBar,attachmentList);
function renderAttachments(){attachmentList.replaceChildren();for(const file of attachedFiles){const chip=node('div',undefined,'attachment-chip');chip.append(node('span',file.name+' ('+Math.ceil(file.bytes/1024)+' KB)'));const remove=node('button','Remove');remove.type='button';remove.disabled=busy;remove.onclick=()=>{attachedFiles=attachedFiles.filter(row=>row.id!==file.id);renderAttachments();};chip.append(remove);attachmentList.append(chip);}}
async function uploadFiles(files){
  if(busy)return;
  if(attachedFiles.length+files.length>8){$('error').textContent='Attach at most eight files per message.';return;}
  setBusy(true);$('error').textContent='';
  try{if(!session)await newSession();for(const file of files){if(!file.size||file.size>50*1024*1024)throw Error('Use nonempty files up to 50 MiB. Keep full libraries on the host and supply their paths.');uploadButton.textContent='Uploading '+file.name+'...';const response=await fetch('/api/uploads?session='+encodeURIComponent(session)+'&name='+encodeURIComponent(file.name),{method:'POST',headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||''),'Content-Type':'application/octet-stream'},body:file});const record=await response.json();if(!response.ok)throw Error(record.error||'Upload failed');attachedFiles.push(record);renderAttachments();}}
  catch(e){$('error').textContent=e.message;}
  finally{uploadButton.textContent='Attach files';fileInput.value='';setBusy(false);renderAttachments();}
}
uploadButton.onclick=()=>fileInput.click();fileInput.onchange=()=>uploadFiles([...fileInput.files]);
for(const event of ['dragenter','dragover'])$('compose').addEventListener(event,e=>{if(e.dataTransfer.types.includes('Files')){e.preventDefault();$('compose').classList.add('drop-active');}});
$('compose').addEventListener('dragleave',e=>{if(!$('compose').contains(e.relatedTarget))$('compose').classList.remove('drop-active');});
$('compose').addEventListener('drop',e=>{e.preventDefault();$('compose').classList.remove('drop-active');uploadFiles([...e.dataTransfer.files]);});
document.addEventListener('click',async e=>{const button=e.target.closest('[data-copy-code]');if(!button)return;try{await navigator.clipboard.writeText(button.closest('.code-block').querySelector('code').textContent);button.textContent='Copied';setTimeout(()=>button.textContent='Copy code',1500);}catch(error){$('error').textContent='Clipboard unavailable. Select and copy the code.';}});
const projectBox=node('section',undefined,'project-box');projectBox.setAttribute('aria-label','Research project');
const projectLabel=node('label','TARGET PROJECT');projectLabel.htmlFor='project-select';
const projectSelect=node('select');projectSelect.id='project-select';
const projectCreate=node('button','+ New project');projectCreate.id='project-create';projectCreate.type='button';
const memoryPanel=node('details',undefined,'project-memory');memoryPanel.append(node('summary','Project knowledge'));
const memoryBody=node('div');memoryPanel.append(memoryBody);projectBox.append(projectLabel,projectSelect,projectCreate,memoryPanel);$('new').before(projectBox);
const knowledgeButton=node('button','Knowledge libraries');knowledgeButton.type='button';projectBox.append(knowledgeButton);
const knowledgeDialog=node('dialog',undefined,'project-dialog direction-dialog');document.body.append(knowledgeDialog);
knowledgeButton.onclick=async()=>{
  try{
    if(!session)await newSession();const knowledgeSession=session;knowledgeDialog.replaceChildren();
    const knowledgeCall=arguments_=>api('/api/knowledge',{session:knowledgeSession,arguments:arguments_});
    const close=node('button','Close');close.onclick=()=>knowledgeDialog.close();
    const state=node('p'),notice=node('p'),output=node('div',undefined,'knowledge-results');notice.setAttribute('aria-live','polite');output.setAttribute('aria-live','polite');
    let generation=0,currentArgs=null;
    async function status(){await knowledgeCall({operation:'status'});state.textContent='Public documents can be used by the cloud model. Internal documents stay in this local interface. Cloud use for the rest of this project remains available.';}
    const collection=node('select');collection.setAttribute('aria-label','Document classification');collection.append(new Option('Internal / group-only','internal'),new Option('Public literature','public'));
    const source=node('input');source.type='url';source.placeholder='Public source URL (required for public documents)';source.setAttribute('aria-label','Public source URL');
    const confirmation=node('input');confirmation.type='checkbox';const confirmationLabel=node('label','I confirm this document and its title/source contain only public information.');confirmationLabel.prepend(confirmation);
    const file=node('input');file.type='file';file.accept='.pdf,.txt,.md,.csv,.tsv,.json';file.setAttribute('aria-label','Knowledge document');
    const ingest=node('button','Import into selected library');ingest.onclick=async()=>{
      ingest.disabled=true;
      for(const control of [collection,source,confirmation,file])control.disabled=true;
      try{
        if(!file.files.length)throw Error('Choose a document.');
        const chosen=file.files[0];if(chosen.size>20*1024*1024)throw Error('Use a document no larger than 20 MiB.');
        if(collection.value==='public'&&(!confirmation.checked||!source.value))throw Error('Supply a public source URL and confirm the public classification.');
        const query=new URLSearchParams({session:knowledgeSession,name:chosen.name,collection:collection.value,source:source.value,public_confirmed:String(confirmation.checked)});
        const r=await fetch('/api/knowledge-upload?'+query,{method:'POST',headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||''),'Content-Type':'application/octet-stream'},body:chosen});
        const result=await r.json();if(!r.ok)throw Error(result.error||'Import failed');
        notice.textContent=(result.status==='already_indexed'?'Existing document verified and reindexed: ':'Document indexed: ')+result.document_id+(result.archived?' (remains archived)':'');
        if(result.previous_public_copy_archived)notice.textContent+=' Its identical public copy was archived. Earlier chat/provider disclosures are not erased.';
        file.value='';confirmation.checked=false;await load({operation:'list',collection:collection.value,limit:10,offset:0,include_archived:archived.checked});
      }catch(e){notice.textContent=e.message;}finally{ingest.disabled=false;for(const control of [collection,source,confirmation,file])control.disabled=false;}
    };
    const query=node('input');query.placeholder='Search words in the selected library';query.setAttribute('aria-label','Local knowledge query');
    const archived=node('input');archived.type='checkbox';const archivedLabel=node('label','Include archived documents (local view only)');archivedLabel.prepend(archived);
    async function load(args){
      const call=++generation;currentArgs={...args};output.replaceChildren(node('p','Reading the selected local library...'));
      try{
        const result=await knowledgeCall(args);if(call!==generation)return;output.replaceChildren();
        if(!result.results.length)output.append(node('p','No matching passages or documents. Try another query or page.'));
        if(args.operation==='read'){
          const page=node('input');page.type='number';page.min='1';page.placeholder='All pages';page.value=args.page||'';page.setAttribute('aria-label','Document page');
          const go=node('button','Read page');go.onclick=()=>load({...args,offset:0,page:page.value?Number(page.value):undefined});output.append(node('p','Page number (leave blank to read across the document)'),page,go);
        }
        for(const row of result.results){
          const card=node('article',undefined,'knowledge-card');card.append(node('h4',((row.bibliography||{}).display_title||row.title)+(row.archived?' (archived)':'')),node('div',row.citation||row.id,'path'));
          if(row.source)card.append(node('p','Source: '+row.source));
          if(row.page)card.append(node('p','Page '+row.page+' · Passage '+(row.part+1)));
          else card.append(node('p','Last indexed page: '+row.last_indexed_page));
          if(row.text)card.append(node('p',row.text,'knowledge-passage'));
          const read=node('button',row.page?'Read this page':'Read document');read.onclick=()=>load({operation:'read',collection:args.collection,document_id:row.id,page:row.page,limit:5,offset:0,include_archived:archived.checked});
          const change=node('button',row.archived?'Restore':'Archive');change.onclick=async()=>{change.disabled=true;try{const receipt=await knowledgeCall({operation:row.archived?'restore':'archive',collection:args.collection,document_id:row.id});notice.textContent=receipt.status+'. '+receipt.history;await load({...currentArgs,offset:0});}catch(e){notice.textContent=e.message;change.disabled=false;}};
          card.append(read,change);
          const edit=node('details');edit.append(node('summary','Citation details'));
          const fields={};for(const [key,caption] of [['display_title','Short paper title'],['authors','Author label (optional)'],['year','Publication year (optional)']]){
            const input=node('input');input.value=(row.bibliography||{})[key]||'';input.placeholder=caption;input.setAttribute('aria-label',caption);fields[key]=input;edit.append(input);
          }
          const save=node('button','Save citation details');save.onclick=async()=>{save.disabled=true;try{
            await knowledgeCall({operation:'bibliography',collection:args.collection,document_id:row.id,...Object.fromEntries(Object.entries(fields).map(([k,v])=>[k,v.value]))});
            notice.textContent='Citation details saved. Original document content is unchanged.';await load(currentArgs);
          }catch(e){notice.textContent=e.message;}finally{save.disabled=false;}};
          edit.append(save);card.append(edit);output.append(card);
        }
        const nav=node('div');if(result.offset>0){const previous=node('button','Previous');previous.onclick=()=>load({...args,offset:Math.max(0,result.offset-result.limit)});nav.append(previous);}
        if(result.has_more){const next=node('button','Next');next.onclick=()=>load({...args,offset:result.next_offset});nav.append(next);}
        output.append(nav);
      }catch(e){if(call===generation)output.replaceChildren(node('p',e.message,'task-error'));}
    }
    const search=node('button','Search locally');search.onclick=()=>load({operation:'search',collection:collection.value,query:query.value,limit:5,offset:0,include_archived:archived.checked});
    const list=node('button','List documents');list.onclick=()=>load({operation:'list',collection:collection.value,limit:10,offset:0,include_archived:archived.checked});
    const savedReviews=node('button','Saved public reviews');savedReviews.disabled=true;
    const indexReviews=node('button','Index previous public reviews');indexReviews.disabled=true;
    const reviewCall=arguments_=>api('/api/literature-reviews',{session:knowledgeSession,arguments:arguments_});
    async function reviewList(offset=0){
      const call=++generation;output.replaceChildren(node('p','Reading saved project reviews...'));
      try{
        const result=await reviewCall({operation:'list',offset,limit:5});if(call!==generation)return;output.replaceChildren();
        if(!result.reviews.length)output.append(node('p','No available saved public reviews. New paper reviews are recorded automatically; older ones can be indexed below.'));
        for(const review of result.reviews){
          const card=node('article',undefined,'knowledge-card');card.append(node('h4',review.papers.map(p=>p.title).join(' / ')),node('p',review.verification));
          const read=node('button','Read saved review');read.onclick=()=>readReview(review.review_id,0);card.append(read);output.append(card);
        }
        if(offset){const back=node('button','Previous reviews');back.onclick=()=>reviewList(Math.max(0,offset-5));output.append(back);}
        if(result.has_more){const next=node('button','Next reviews');next.onclick=()=>reviewList(result.next_offset);output.append(next);}
      }catch(e){if(call===generation)output.replaceChildren(node('p',e.message,'task-error'));}
    }
    async function readReview(reviewId,offset){
      const call=++generation;
      try{
        const result=await reviewCall({operation:'read',review_id:reviewId,offset});if(call!==generation)return;
        output.replaceChildren(node('p',result.verification),markdown(result.summary));
        const back=node('button','Back to saved reviews');back.onclick=()=>reviewList();output.append(back);
        if(result.has_more){const next=node('button','Continue saved review');next.onclick=()=>readReview(reviewId,result.next_offset);output.append(next);}
      }catch(e){if(call===generation)output.replaceChildren(node('p',e.message,'task-error'));}
    }
    savedReviews.onclick=()=>reviewList();let indexOffset=0;
    indexReviews.onclick=async()=>{indexReviews.disabled=true;const call=++generation;try{
      const result=await reviewCall({operation:'index',offset:indexOffset,limit:10});if(call!==generation)return;indexOffset=result.next_offset||0;
      notice.textContent=result.indexed_reviews+' public review(s) indexed. '+(result.has_more?'Click again to inspect the next batch.':'Index pass finished.');
      await reviewList();
    }catch(e){notice.textContent=e.message;}finally{indexReviews.disabled=collection.value!=='public';}};
    collection.onchange=()=>{generation++;query.value='';confirmation.checked=false;currentArgs=null;output.replaceChildren();notice.textContent='';savedReviews.disabled=indexReviews.disabled=collection.value!=='public';};
    query.onkeydown=e=>{if(e.key==='Enter'){e.preventDefault();search.click();}};
    knowledgeDialog.append(close,node('h2','Public and internal knowledge'),state,node('p','Internal search results are displayed here only and are not added to chat or project memory. Do not paste internal content into ordinary cloud chat. This classification applies to documents imported here.'),collection,source,confirmationLabel,file,ingest,notice,node('h3','Local retrieval'),archivedLabel,query,search,list,savedReviews,indexReviews,output);
    await status();knowledgeDialog.showModal();
  }catch(e){$('error').textContent=e.message;}
};
const directionButton=node('button','Project direction & learning');directionButton.type='button';projectBox.append(directionButton);
const directionDialog=node('dialog',undefined,'project-dialog direction-dialog'),directionBody=node('div');directionDialog.append(directionBody);document.body.append(directionDialog);
let guidanceState=null,guidanceSession=null;
async function guidanceCall(arguments_){return api('/api/project-guidance',{session:guidanceSession,arguments:arguments_});}
function guidanceAction(label,operation,after){const b=node('button',label);b.type='button';b.onclick=async()=>{b.disabled=true;try{await guidanceCall(operation);if(after)after();await loadGuidance();await refresh(true);}catch(e){directionBody.prepend(node('p',e.message,'task-error'));}finally{b.disabled=false;}};return b;}
async function loadGuidance(){
  guidanceState=await guidanceCall({operation:'read'});const m=guidanceState.mission;directionBody.replaceChildren();
  const close=node('button','Close');close.type='button';close.onclick=()=>directionDialog.close();directionBody.append(close,node('h2','Project direction'),node('p','Revision '+m.revision+' · Target '+m.target),node('p','Original goal: '+m.baseline.objective));
  const form=node('form');const fields={};
  for(const [key,label,value] of [['objective','Current goal',m.objective],['target','Target',m.target],['excluded_work','Outside this project',m.excluded_work],['constraints','Constraints (one per line)',m.constraints.join('\n')]]){const l=node('label',label),input=node('textarea');input.value=value;input.required=['objective','target'].includes(key);fields[key]=input;l.append(input);form.append(l);}
  const reqs=node('div');form.append(node('h3','Deliverables and acceptance'),reqs);const requirements=[];
  for(const r of m.requirements){const box=node('fieldset');box.append(node('legend',r.id));const text=node('textarea'),acceptance=node('textarea');text.value=r.text;acceptance.value=r.acceptance;for(const [label,input] of [['Required outcome',text],['Definition of done',acceptance]]){const l=node('label',label);l.append(input);box.append(l);}requirements.push({r,text,acceptance});reqs.append(box);}
  const safeguard=node('details');safeguard.append(node('summary','Execution safeguards'),node('pre',JSON.stringify(m.rules,null,2)),node('p','Ask in chat to propose changes to exact sequences, seeds, protein copies, library coverage or permitted operations. Review the concrete proposal here before applying it.'));form.append(safeguard);
  const save=node('button','Save my project changes');save.type='submit';form.append(save);form.onsubmit=async e=>{e.preventDefault();save.disabled=true;try{const mission={objective:fields.objective.value,target:fields.target.value,excluded_work:fields.excluded_work.value,constraints:fields.constraints.value.split('\n').filter(s=>s.trim()),rules:m.rules,requirements:requirements.map(({r,text,acceptance})=>({...r,text:text.value,acceptance:acceptance.value}))};await guidanceCall({operation:'apply',base_revision:m.revision,mission});await loadGuidance();}catch(error){form.prepend(node('p',error.message,'task-error'));}finally{save.disabled=false;}};directionBody.append(form);
  directionBody.append(node('h3','Proposed changes'));
  if(!guidanceState.proposals.length)directionBody.append(node('p','No pending changes. The agent cannot change this goal by itself.'));
  for(const p of guidanceState.proposals){const d=node('details');d.append(node('summary','Review proposal '+p.id),node('pre',JSON.stringify(p.mission,null,2)),guidanceAction('Apply this proposal',{operation:'adopt',proposal_id:p.id}));directionBody.append(d);}
  directionBody.append(node('h3','What this project has learned'));
  const acceptBox=node('details');acceptBox.append(node('summary','Record my requirement acceptance'));
  const requirementSelect=node('select');requirementSelect.setAttribute('aria-label','Requirement');for(const r of m.requirements)requirementSelect.append(new Option(r.id+' · '+r.text,r.id));
  const evidenceSelect=node('select');evidenceSelect.setAttribute('aria-label','Evidence');for(const e of guidanceState.learning.filter(e=>e.kind==='evidence'))evidenceSelect.append(new Option(e.id+' · task '+e.task_id,e.id));
  const acceptanceReason=node('textarea');acceptanceReason.placeholder='Why this evidence meets the requirement';acceptanceReason.setAttribute('aria-label','Acceptance rationale');
  const acceptButton=node('button','Record my acceptance');acceptButton.type='button';acceptButton.onclick=async()=>{acceptButton.disabled=true;try{await guidanceCall({operation:'accept_requirement',base_revision:m.revision,requirement_id:requirementSelect.value,evidence_ids:[evidenceSelect.value],reason:acceptanceReason.value});await loadGuidance();}catch(e){acceptBox.append(node('p',e.message,'task-error'));}finally{acceptButton.disabled=false;}};
  acceptBox.append(requirementSelect,evidenceSelect,acceptanceReason,node('p','Acceptance applies to this mission revision and does not establish experimental validity.'),acceptButton);directionBody.append(acceptBox);
  for(const entry of [...guidanceState.learning].reverse()){if(entry.kind==='task_event')continue;const d=node('details');d.append(node('summary',entry.kind+' · '+(entry.text||entry.task_id||entry.id)+(entry.resolved?' · archived':'')),node('pre',JSON.stringify(entry,null,2)));if(!entry.resolved&&['constraint','preference','correction','observation','hypothesis','next_step'].includes(entry.kind))d.append(guidanceAction('Archive this feedback',{operation:'resolve_feedback',entry_id:entry.id}));directionBody.append(d);}
  directionBody.append(node('p','Dossier: '+guidanceState.dossier.path,'path'),node('p','Dossier status: '+guidanceState.dossier.status));
  const cards=node('details');cards.append(node('summary','Task alignment and reuse records'),node('pre',JSON.stringify(guidanceState.task_cards.map(({mission_snapshot,output_hashes,...c})=>c),null,2)));directionBody.append(cards);
}
directionButton.onclick=async()=>{try{if(!session)await newSession();guidanceSession=session;await loadGuidance();directionDialog.showModal();}catch(e){$('error').textContent=e.message;}};
const projectDialog=node('dialog',undefined,'project-dialog');const projectForm=node('form');const projectError=node('p',undefined,'task-error');projectError.setAttribute('role','alert');
projectForm.append(node('h2','New target project'),node('p','Conversations, tasks and saved knowledge belong to this target.'));
for(const [key,label] of [['name','Project name'],['target','Target / accession'],['organism','Organism (optional)'],['objective','Research objective']]){const l=node('label',label);l.htmlFor='project-'+key;const input=node(key==='objective'?'textarea':'input');input.id='project-'+key;input.name=key;input.required=key!=='organism';input.maxLength=2000;projectForm.append(l,input);}
const createSubmit=node('button','Create project');createSubmit.type='submit';const cancelProject=node('button','Cancel');cancelProject.type='button';cancelProject.onclick=()=>projectDialog.close();projectForm.append(projectError,createSubmit,cancelProject);projectDialog.append(projectForm);document.body.append(projectDialog);
projectCreate.onclick=()=>{projectError.textContent='';projectDialog.showModal();};
async function openProject(arguments_){setBusy(true);try{const result=await api('/api/project',{session,arguments:arguments_});projectId=result.project_id;selectSession(result.session);await refresh(true);return true;}catch(e){$('error').textContent=e.message;return false;}finally{setBusy(false);}}
projectSelect.onchange=()=>openProject({operation:'open',project_id:projectSelect.value});
projectForm.onsubmit=async e=>{e.preventDefault();projectError.textContent='';createSubmit.disabled=true;const args={operation:'create'};for(const key of ['name','target','organism','objective'])if($('project-'+key).value.trim())args[key]=$('project-'+key).value.trim();try{if(await openProject(args))projectDialog.close();else projectError.textContent=$('error').textContent;}finally{createSubmit.disabled=false;}};
async function controlResearchFlow(id,operation){setBusy(true);try{await api('/api/research-workflow',{session,arguments:{operation,workflow_id:id}});await refresh(true);}catch(e){$('error').textContent=e.message;}finally{setBusy(false);}}
function renderProjects(state){
  projectId=state.project_id;projectSelect.replaceChildren();
  for(const p of state.projects||[]){const option=node('option',p.name+(p.target?' / '+p.target:''));option.value=p.project_id;option.selected=p.project_id===projectId;projectSelect.append(option);}
  memoryBody.replaceChildren();const memory=state.project_memory;if(!memory)return;
  memoryBody.append(node('p',memory.project.objective||'No objective recorded.'),node('p',memory.total_tasks+' project tasks'));
  for(const note of memory.notes||[])memoryBody.append(node('p',note.kind+': '+note.text+(note.task_id?' (task '+note.task_id+')':'')));
  for(const set of memory.candidate_sets||[])memoryBody.append(node('p',set.name+' / '+set.count+' molecules'),node('div',set.set_id,'path'));
  for(const flow of memory.workflows||[]){const entry=node('div',undefined,'project-flow');entry.append(node('strong',flow.status+': '+flow.objective),node('div',flow.id,'path'));if(flow.error)entry.append(node('p',flow.error,'task-error'));for(const step of flow.nodes)entry.append(node('p',step.id+': '+step.status+(step.task_id?' / '+step.task_id:'')));if(['active','blocked','paused'].includes(flow.status)){const action=flow.status==='active'?'pause':'resume';const b=node('button',action==='pause'?'Pause continuation':'Resume continuation');b.disabled=busy;b.onclick=()=>controlResearchFlow(flow.id,action);entry.append(b);}memoryBody.append(entry);}
}
refresh();setInterval(refresh,2000);

$('attach').onclick=async()=>{try{if(!session)await newSession();const kind=$('attachment-kind').value;if(['kras_export','af3_batch','af3_prediction'].includes(kind))await api('/api/workflow-results',{session,arguments:{operation:'attach',kind,report:$('existing').value.trim()}});else if(kind==='equiscore')await api('/api/block-campaign',{session,arguments:{operation:'import',report:$('existing').value.trim()}});else if(kind==='docking')await api('/api/block-results',{session,arguments:{operation:'attach',report:$('existing').value.trim()}});else await api('/api/import',{session,plan:$('existing').value.trim()});await refresh();}catch(e){$('error').textContent=e.message;}};

function renderViewer(viewer){
if(!viewer||(!viewer.connected&&!viewer.latest))return;
const card=node('article',undefined,'task');card.append(node('strong','Desktop PyMOL'),node('p',viewer.connected?'Connected':'Not connected'));
if(viewer.latest){const result=viewer.latest;card.append(node('pre',JSON.stringify(result,null,2)));
for(const [kind,path] of Object.entries(result.artifacts||{})){const button=node('button','Download '+kind);button.onclick=async()=>{
try{const url='/api/viewer/artifact?session='+encodeURIComponent(session)+'&id='+encodeURIComponent(result.id)+'&artifact='+encodeURIComponent(kind);
const response=await fetch(url,{headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||'')}});if(!response.ok)throw Error('Artifact unavailable');
const blob=URL.createObjectURL(await response.blob());const a=node('a');a.href=blob;a.download=path.split(/[\\/]/).pop();a.click();setTimeout(()=>URL.revokeObjectURL(blob),1000);
}catch(e){$('error').textContent=e.message;}};card.append(button);}}
$('tasks').prepend(card);
}

function renderAgent(runs){
for(const run of [...(runs||[])].reverse()){
const card=node('article',undefined,'task');
card.append(node('strong','Research request: '+(run.status==='pending'?'submitted; see task status':run.status)),node('p',(run.request||'').slice(0,180)));
for(const step of run.steps||[])card.append(node('div',step.tool+' / '+step.status+' ? '+(step.error||step.purpose||'')));
const button=node('button','Download tool trace');button.onclick=async()=>{
try{const data=await api('/api/agent/trace?session='+encodeURIComponent(session)+'&id='+encodeURIComponent(run.id));
const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'}));const a=node('a');a.href=url;a.download=run.id+'.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}catch(e){$('error').textContent=e.message;}};card.append(button);$('tasks').prepend(card);
}}

function renderChains(chains){
for(const chain of [...(chains||[])].reverse()){
const card=node('article',undefined,'task');
card.append(node('strong','Structure workflow: '+chain.status),node('div','Chain '+chain.id,'path'),node('p','Goal: '+chain.goal),node('div','Current task: '+chain.current_task,'path'));
if(chain.error)card.append(node('p',chain.error));
for(const step of chain.history||[])card.append(node('div',step.intent+': '+step.source_task+' -> '+step.task_id));
card.append(node('p','Chain controls affect future stages; current tasks are controlled separately.'));
$('tasks').prepend(card);
}}

// Examples fill the composer; only Send submits a request.
for(const button of document.querySelectorAll('[data-prompt]'))button.onclick=()=>{$('prompt').value=button.dataset.prompt;$('prompt').focus();};
$('attachment-kind').append(new Option('Completed EquiScore analysis','equiscore'));
$('attachment-kind').onchange=()=>{const kind=$('attachment-kind').value;const external=['kras_export','af3_batch','af3_prediction'].includes(kind);$('existing-label').textContent=external?'Report or manifest JSON path on this host':kind==='equiscore'?'EquiScore analysis report.json path':kind==='docking'?'PLANTS report.json path':'Existing plan.json path';$('existing').placeholder=external?'/configured/import/root/report.json':kind==='equiscore'?'/project/runs/.../equiscore-analysis/report.json':kind==='docking'?'/project/runs/.../plants-final/report.json':'/project/runs/PROMPT-.../plan.json';};
const coverage=$('coverage-menu');let coveragePinned=false;
for(const [value,label] of [['kras_export','KRAS pool export report'],['af3_batch','AF3 ligand batch report'],['af3_prediction','AF3 prediction manifest']])$('attachment-kind').append(new Option(label,value));
const workflowResults=node('div');workflowResults.id='workflow-results';$('docking-results').before(workflowResults);
const sourceChecks=new Map();
let workflowResultsKey=null;
function renderWorkflowResults(rows){
const key=JSON.stringify(rows||[]);if(key===workflowResultsKey)return;workflowResultsKey=key;
workflowResults.replaceChildren();
for(const row of rows||[]){
const card=node('article',undefined,'task');
card.append(node('strong',row.kind+' / '+row.summary.status),node('div',row.source,'path'),node('p','Saved report evidence. Model and MOL2 file contents have not been verified.'));
if(row.summary.targets)for(const [name,data] of Object.entries(row.summary.targets))card.append(node('p',name+': '+data.processed.toLocaleString()+' processed; '+data.exported.toLocaleString()+' exported; '+data.quarantined.toLocaleString()+' quarantined; '+data.batch_count+' MOL2 files.'));
if(row.summary.mode)card.append(node('p',row.summary.recorded_jobs+' recorded jobs; '+row.summary.generated_models+' reported models; '+row.summary.failed_jobs+' failed jobs.'));
const details=node('details');details.append(node('summary','Counts, confidence and reported file paths'),node('pre',JSON.stringify(row.summary,null,2)));card.append(details);
const check=node('button','Check source report');const feedback=node('p',sourceChecks.get(row.id)||'');check.onclick=async()=>{check.disabled=true;try{const result=await api('/api/workflow-results',{session,arguments:{operation:'check',attachment_id:row.id}});const message=result.source_status+' (checked '+new Date(result.checked_at*1000).toLocaleTimeString()+')';sourceChecks.set(row.id,message);feedback.textContent=message;}catch(e){feedback.textContent=e.message;}finally{check.disabled=false;}};
card.append(check,feedback);workflowResults.append(card);
}
}
coverage.addEventListener('pointerenter',event=>{if(event.pointerType==='mouse')coverage.open=true;});
coverage.addEventListener('pointerleave',()=>{if(!coveragePinned&&!coverage.contains(document.activeElement))coverage.open=false;});
coverage.addEventListener('focusin',()=>{coverage.open=true;});
coverage.addEventListener('focusout',event=>{if(!coveragePinned&&!coverage.contains(event.relatedTarget))coverage.open=false;});
coverage.querySelector('summary').addEventListener('click',event=>{event.preventDefault();coveragePinned=!coveragePinned;coverage.open=coveragePinned;});
coverage.addEventListener('keydown',event=>{if(event.key==='Escape'){coveragePinned=false;coverage.open=false;coverage.querySelector('summary').focus();coverage.open=false;}});
document.addEventListener('pointerdown',event=>{if(!coverage.contains(event.target)){coveragePinned=false;coverage.open=false;}});
function renderDocking(rows){
$('docking-results').replaceChildren();
for(const row of rows||[]){const card=node('article',undefined,'task');card.append(node('strong','Attached docking result'),node('div',row.id,'path'),node('p',row.summary.status+' / '+row.summary.scored.toLocaleString()+' scores / '+row.summary.failed_jobs+' failed jobs'),node('p',row.verification),node('div',row.report,'path'));const button=node('button','Analyze saved scores');button.onclick=async()=>{button.disabled=true;try{await api('/api/block-results',{session,arguments:{operation:'analyze',attachment_id:row.id}});await refresh();}catch(e){$('error').textContent=e.message;}finally{button.disabled=false;}};card.append(button);$('docking-results').append(card);}
}

// Hover previews results; clicking pins the drawer for sustained interaction.
const resultsMenu=$('results-menu');let resultsPinned=false,resultsCloseTimer;
function closeResults(restoreFocus=false){clearTimeout(resultsCloseTimer);resultsPinned=false;if(restoreFocus)resultsMenu.querySelector('summary').focus();resultsMenu.open=false;}
resultsMenu.addEventListener('pointerenter',event=>{clearTimeout(resultsCloseTimer);if(event.pointerType==='mouse')resultsMenu.open=true;});
resultsMenu.addEventListener('pointerleave',()=>{clearTimeout(resultsCloseTimer);resultsCloseTimer=setTimeout(()=>{if(!resultsPinned&&!resultsMenu.contains(document.activeElement))resultsMenu.open=false;},200);});
resultsMenu.addEventListener('focusin',()=>{clearTimeout(resultsCloseTimer);resultsMenu.open=true;});
resultsMenu.addEventListener('focusout',event=>{if(!resultsPinned&&!resultsMenu.contains(event.relatedTarget))resultsMenu.open=false;});
resultsMenu.querySelector('summary').addEventListener('click',event=>{event.preventDefault();clearTimeout(resultsCloseTimer);resultsPinned=!resultsPinned;resultsMenu.open=resultsPinned;});
resultsMenu.addEventListener('keydown',event=>{if(event.key==='Escape'){event.preventDefault();closeResults(true);}});
$('close-results').onclick=()=>closeResults(true);
document.addEventListener('pointerdown',event=>{if(!resultsMenu.contains(event.target))closeResults();});

// Explicit controls use the same owned tool contract as research chat.
const campaignPanel=node('details',undefined,'attach-menu');
campaignPanel.append(node('summary','Block workflows (development deferred)'));
const campaignFields={};
function campaignField(key,label,options){
  const control=node(options?'select':'input');control.id='campaign-'+key;
  if(options)for(const [value,text] of options)control.append(new Option(text,value));
  const caption=node('label',label);caption.htmlFor=control.id;
  campaignPanel.append(caption,control);campaignFields[key]=control;return control;
}
campaignField('task_id','Imported analysis task',[]);
campaignField('query_task_id','Confirmed MDM2 query task',[]);
campaignField('scheme','Partition',[['E094','E094'],['E095','E095'],['E096','E096']]);
campaignField('receptor','Receptor',[['1RV1_B','1RV1_B'],['7NA2_A','7NA2_A'],['5J7F_A','5J7F_A']]);
campaignField('method','Block selection',[['union','Union of both scores'],['intersection','Shared leading blocks'],['equiscore','EquiScore'],['chemplp','ChemPLP']]);
campaignField('blocks','Leading blocks per score',[['5','Top 5'],['10','Top 10']]);
const conformerInput=campaignField('conformers','Conformers to export (multiple per molecule allowed)');
conformerInput.type='number';conformerInput.min='1';conformerInput.max='1000000';conformerInput.value='100000';
const campaignOutput=node('pre');campaignOutput.setAttribute('aria-live','polite');
async function campaignCall(arguments_){if(!session)await newSession();return api('/api/block-campaign',{session,arguments:arguments_});}
function campaignButton(label,callback){const button=node('button',label);button.type='button';button.onclick=async()=>{button.disabled=true;try{await callback();}catch(e){campaignOutput.textContent=e.message;}finally{button.disabled=false;}};campaignPanel.append(button);}
campaignButton('Load saved analyses and queries',async()=>{
  const data=await campaignCall({operation:'list'});
  campaignFields.task_id.replaceChildren(new Option('Select an imported analysis',''));
  campaignFields.query_task_id.replaceChildren(new Option('Select the confirmed query',''));
  for(const task of data.tasks.filter(t=>t.status==='complete'&&t.request.startsWith('Import saved EquiScore')))
    campaignFields.task_id.append(new Option(task.id,task.id));
  for(const query of data.query_tasks)campaignFields.query_task_id.append(new Option(query.task_id+' / '+query.kind,query.task_id));
  campaignOutput.textContent=JSON.stringify(data,null,2);
});
function campaignSelection(operation){return {operation,task_id:campaignFields.task_id.value,scheme:campaignFields.scheme.value,receptor:campaignFields.receptor.value,method:campaignFields.method.value,blocks:Number(campaignFields.blocks.value)};}
campaignButton('Preview leading blocks',async()=>{campaignOutput.textContent=JSON.stringify(await campaignCall(campaignSelection('preview')),null,2);});
campaignButton('Search and dock conformers',async()=>{
  const args={...campaignSelection('start'),query_task_id:campaignFields.query_task_id.value,conformers:Number(conformerInput.value),dock:true};
  campaignOutput.textContent=JSON.stringify(await campaignCall(args),null,2);await refresh();
});
campaignPanel.append(node('p','Blocks: mean of Top-10 distinct molecules. Search budget: conformers across the selected blocks. Docking uses the selected receptor.','attach-note'),campaignOutput);
$('docking-results').before(campaignPanel);
