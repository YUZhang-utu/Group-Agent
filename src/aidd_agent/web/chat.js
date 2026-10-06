const $=id=>document.getElementById(id);
const fragment=new URLSearchParams(location.hash.slice(1));
if(fragment.get('token')){sessionStorage.setItem('aidd-token',fragment.get('token'));history.replaceState(null,'',location.pathname);}
if(fragment.get('session'))sessionStorage.setItem('aidd-session',fragment.get('session'));
let session=sessionStorage.getItem('aidd-session'),last='',polling=false;
async function api(path,body){const r=await fetch(path,{method:body?'POST':'GET',headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||''),'Content-Type':'application/json'},body:body?JSON.stringify(body):undefined});const data=await r.json();if(!r.ok)throw Error(data.error||'Request failed');return data;}
function node(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;}
async function newSession(){const r=await api('/api/session',{});session=r.session;sessionStorage.setItem('aidd-session',session);last='';await refresh();}
async function send(text){if($('send').disabled)return;$('error').textContent='';$('send').disabled=true;try{if(!session)await newSession();await api('/api/message',{session,text,provider:$('provider').value});if($('prompt').value===text)$('prompt').value='';await refresh();}catch(e){$('error').textContent=e.message;}finally{$('send').disabled=false;}}
function render(state){$('task-count').textContent=state.tasks.length;$('task-count').hidden=!state.tasks.length;$('compute').textContent=state.compute_enabled?'Compute enabled':'Preparation only';$('sessions').replaceChildren();for(const s of state.sessions){const b=node('button',s.title,s.id===session?'selected':'');b.onclick=()=>{session=s.id;sessionStorage.setItem('aidd-session',session);last='';refresh();};$('sessions').append(b);}
const key=JSON.stringify(state.messages);if(key!==last){$('messages').replaceChildren();for(const m of state.messages){const box=node('div',undefined,'message '+m.role);box.append(node('span',m.role==='user'?'You':'MEDCHEM Agent', 'role'),node('span',m.text));$('messages').append(box);}if(!state.messages.length)$('messages').append(node('p','Ask a research question.','empty'));$('messages').scrollTop=$('messages').scrollHeight;last=key;}
$('coverage').replaceChildren();for(const w of state.workflows){const n=node('div',undefined,'coverage-item');n.append(node('strong',w.name),node('div',w.status),node('div',w.scope));$('coverage').append(n);}
$('tasks').replaceChildren();if(!state.tasks.length)$('tasks').append(node('p','Your tasks will appear here.','empty'));for(const t of [...state.tasks].reverse()){const card=node('article',undefined,'task '+t.status);card.append(node('div',t.status+' / '+t.provider,'task-status'),node('p',t.request),node('div','Task '+t.id,'path'));for(const k of ['plan','report','log'])if(t[k]){card.append(node('div',k+': '+t[k],'path'));const b=node('button','Copy '+k+' path');b.onclick=()=>navigator.clipboard.writeText(t[k]).catch(e=>$('error').textContent=e.message);card.append(b);}if(t.error)card.append(node('p',t.error));if(t.details){const d=node('details');d.append(node('summary','Steps and result artifacts'),node('pre',JSON.stringify(t.details,null,2)));card.append(d);}if(t.progress){const d=node('details');d.append(node('summary','Execution log'),node('pre',t.progress));card.append(d);}const action=['queued','planning','running'].includes(t.status)?'Cancel':'Resume';const b=node('button',action+' task');b.onclick=()=>send('/'+action.toLowerCase()+' '+t.id);card.append(b);$('tasks').append(card);}}
async function refresh(){if(polling)return;polling=true;try{const state=await api('/api/state'+(session?'?session='+encodeURIComponent(session):''));render(state);renderDocking(state.docking_results);renderViewer(state.viewer);renderAgent(state.agent_runs);renderChains(state.structure_chains);$('connection').textContent='Connected';}catch(e){$('connection').textContent='Disconnected';$('error').textContent=e.message;}finally{polling=false;}}
$('new').onclick=()=>newSession().catch(e=>$('error').textContent=e.message);$('compose').onsubmit=e=>{e.preventDefault();if($('prompt').value.trim())send($('prompt').value);};$('prompt').onkeydown=e=>{if(e.key==='Enter'&&!e.shiftKey){e.preventDefault();$('compose').requestSubmit();}};
refresh();setInterval(refresh,2000);

$('attach').onclick=async()=>{try{if(!session)await newSession();if($('attachment-kind').value==='equiscore')await api('/api/block-campaign',{session,arguments:{operation:'import',report:$('existing').value.trim()}});else if($('attachment-kind').value==='docking')await api('/api/block-results',{session,arguments:{operation:'attach',report:$('existing').value.trim()}});else await api('/api/import',{session,plan:$('existing').value.trim()});await refresh();}catch(e){$('error').textContent=e.message;}};

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
card.append(node('strong','Research agent: '+run.status),node('p',run.request));
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
$('attachment-kind').onchange=()=>{const kind=$('attachment-kind').value;$('existing-label').textContent=kind==='equiscore'?'EquiScore analysis report.json path':kind==='docking'?'PLANTS report.json path':'Existing plan.json path';$('existing').placeholder=kind==='equiscore'?'/project/runs/.../equiscore-analysis/report.json':kind==='docking'?'/project/runs/.../plants-final/report.json':'/project/runs/PROMPT-.../plan.json';};
const coverage=$('coverage-menu');let coveragePinned=false;
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
campaignPanel.append(node('summary','Leading blocks to 3D search'));
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
