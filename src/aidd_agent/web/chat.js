const $=id=>document.getElementById(id);
const fragment=new URLSearchParams(location.hash.slice(1));
if(fragment.get('token')){sessionStorage.setItem('aidd-token',fragment.get('token'));history.replaceState(null,'',location.pathname);}
if(fragment.get('session'))sessionStorage.setItem('aidd-session',fragment.get('session'));
let session=sessionStorage.getItem('aidd-session'),last='',polling=null, busy=false,lastState='';
let attachedFiles=[];
function markdown(value){const body=node('div',undefined,'markdown');body.innerHTML=AiddMarkdown.render(value);return body;}
function setBusy(value){busy=value;for(const id of ['send','new','provider','upload-files'])if($(id))$(id).disabled=value;for(const b of $('sessions').querySelectorAll('button'))b.disabled=value;}
async function api(path,body){const r=await fetch(path,{method:body?'POST':'GET',headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||''),'Content-Type':'application/json'},body:body?JSON.stringify(body):undefined});const data=await r.json();if(!r.ok)throw Error(data.error||'Request failed');return data;}
function node(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;}
async function newSession(){const r=await api('/api/session',{});session=r.session;attachedFiles=[];renderAttachments();sessionStorage.setItem('aidd-session',session);last='';await refresh();}
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
  if(key!==lastState){render(state);renderDocking(state.docking_results);renderWorkflowResults(state.workflow_results);renderViewer(state.viewer);renderAgent(state.agent_runs);renderChains(state.structure_chains);lastState=key;}
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
