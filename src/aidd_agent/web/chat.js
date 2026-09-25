const $=id=>document.getElementById(id);
const fragment=new URLSearchParams(location.hash.slice(1));
if(fragment.get('token')){sessionStorage.setItem('aidd-token',fragment.get('token'));history.replaceState(null,'',location.pathname);}
if(fragment.get('session'))sessionStorage.setItem('aidd-session',fragment.get('session'));
let session=sessionStorage.getItem('aidd-session'),last='',polling=false;
async function api(path,body){const r=await fetch(path,{method:body?'POST':'GET',headers:{'Authorization':'Bearer '+(sessionStorage.getItem('aidd-token')||''),'Content-Type':'application/json'},body:body?JSON.stringify(body):undefined});const data=await r.json();if(!r.ok)throw Error(data.error||'Request failed');return data;}
function node(tag,text,cls){const n=document.createElement(tag);if(text!==undefined)n.textContent=text;if(cls)n.className=cls;return n;}
async function newSession(){const r=await api('/api/session',{});session=r.session;sessionStorage.setItem('aidd-session',session);last='';await refresh();}
async function send(text){if($('send').disabled)return;$('error').textContent='';$('send').disabled=true;try{if(!session)await newSession();await api('/api/message',{session,text,provider:$('provider').value});if($('prompt').value===text)$('prompt').value='';await refresh();}catch(e){$('error').textContent=e.message;}finally{$('send').disabled=false;}}
function render(state){$('compute').textContent=state.compute_enabled?'Compute enabled':'Preparation only';$('sessions').replaceChildren();for(const s of state.sessions){const b=node('button',s.title,s.id===session?'selected':'');b.onclick=()=>{session=s.id;sessionStorage.setItem('aidd-session',session);last='';refresh();};$('sessions').append(b);}
const key=JSON.stringify(state.messages);if(key!==last){$('messages').replaceChildren();for(const m of state.messages){const box=node('div',undefined,'message '+m.role);box.append(node('span',m.role==='user'?'You':'AIDD', 'role'),node('span',m.text));$('messages').append(box);}if(!state.messages.length)$('messages').append(node('p','Describe a task, ask what is supported, or request the status of your latest run.','empty'));$('messages').scrollTop=$('messages').scrollHeight;last=key;}
$('coverage').replaceChildren();for(const w of state.workflows){const n=node('div',undefined,'coverage-item');n.append(node('strong',w.name),node('div',w.status),node('div',w.scope));$('coverage').append(n);}
$('tasks').replaceChildren();if(!state.tasks.length)$('tasks').append(node('p','No tasks yet. Your plans and results will appear here.','empty'));for(const t of [...state.tasks].reverse()){const card=node('article',undefined,'task '+t.status);card.append(node('div',t.status+' / '+t.provider,'task-status'),node('p',t.request),node('div','Task '+t.id,'path'));for(const k of ['plan','report','log'])if(t[k]){card.append(node('div',k+': '+t[k],'path'));const b=node('button','Copy '+k+' path');b.onclick=()=>navigator.clipboard.writeText(t[k]).catch(e=>$('error').textContent=e.message);card.append(b);}if(t.error)card.append(node('p',t.error));if(t.details){const d=node('details');d.append(node('summary','Steps and result artifacts'),node('pre',JSON.stringify(t.details,null,2)));card.append(d);}if(t.progress){const d=node('details');d.append(node('summary','Execution log'),node('pre',t.progress));card.append(d);}const action=['queued','planning','running'].includes(t.status)?'Cancel':'Resume';const b=node('button',action+' task');b.onclick=()=>send('/'+action.toLowerCase()+' '+t.id);card.append(b);$('tasks').append(card);}}
async function refresh(){if(polling)return;polling=true;try{const state=await api('/api/state'+(session?'?session='+encodeURIComponent(session):''));render(state);renderViewer(state.viewer);renderAgent(state.agent_runs);$('connection').textContent='Connected';}catch(e){$('connection').textContent='Disconnected';$('error').textContent=e.message;}finally{polling=false;}}
$('new').onclick=()=>newSession().catch(e=>$('error').textContent=e.message);$('compose').onsubmit=e=>{e.preventDefault();if($('prompt').value.trim())send($('prompt').value);};$('prompt').onkeydown=e=>{if(e.key==='Enter'&&!e.shiftKey){e.preventDefault();$('compose').requestSubmit();}};
refresh();setInterval(refresh,2000);

$('attach').onclick=async()=>{try{if(!session)await newSession();await api('/api/import',{session,plan:$('existing').value});await refresh();}catch(e){$('error').textContent=e.message;}};

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
