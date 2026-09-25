"""Loopback-only chat workbench with persistent sessions and background jobs."""
from __future__ import annotations
import argparse
import hmac
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import secrets
from urllib.parse import urlsplit, parse_qs

from .chat_agent import ChatAgent
from .prompt_workflow import file_lock


def make_server(agent, port=8765, token=None):
    secret = token or secrets.token_urlsafe(32)
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args): pass
        def send(self, status, value, content_type="application/json"):
            data = json.dumps(value, ensure_ascii=False).encode() if content_type == "application/json" else value
            self.send_response(status)
            self.send_header("Content-Type", content_type+"; charset=utf-8")
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "no-referrer")
            self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' blob:; frame-ancestors 'none'; base-uri 'none'")
            self.end_headers(); self.wfile.write(data)
        def authorized(self):
            return hmac.compare_digest(self.headers.get("Authorization", ""), "Bearer "+secret)
        def do_GET(self):
            path = urlsplit(self.path)
            assets = {"/":("chat.html","text/html"),"/chat.js":("chat.js","text/javascript"),"/chat.css":("chat.css","text/css")}
            if path.path in assets:
                name, mime = assets[path.path]
                self.send(200,(Path(__file__).parent / "web" / name).read_bytes(),mime); return
            if not self.authorized(): self.send(401,{"error":"Open the local access link printed by the server."}); return
            if path.path=='/api/viewer/artifact':
                try:
                    from .project_context import ensure_within
                    import re
                    query=parse_qs(path.query);sid=query['session'][0];agent.session(sid)
                    identifier=query['id'][0];key=query['artifact'][0]
                    if not re.fullmatch(r'\d+-[a-f0-9]{12}',identifier):raise ValueError('Invalid viewer operation ID')
                    root=agent.root/'viewers'/sid
                    receipt=json.loads((root/(identifier+'.result.json')).read_text(encoding='utf-8'))
                    file=ensure_within(Path(receipt['artifacts'][key]),root)
                    mime={'.png':'image/png','.pse':'application/octet-stream','.json':'application/json','.csv':'text/csv','.py':'text/plain'}[file.suffix]
                    if mime=='application/json':self.send(200,json.loads(file.read_text(encoding='utf-8')))
                    else:self.send(200,file.read_bytes(),mime)
                except (ValueError,KeyError,OSError):self.send(400,{'error':'Viewer artifact unavailable'})
                return
            if path.path != "/api/state": self.send(404,{"error":"Not found"}); return
            try: self.send(200,agent.snapshot(parse_qs(path.query).get("session",[None])[0]))
            except ValueError as exc: self.send(400,{"error":str(exc)})
        def do_POST(self):
            if not self.authorized(): self.send(401,{"error":"Invalid access token"}); return
            origin = self.headers.get("Origin")
            if origin and origin not in {f"http://127.0.0.1:{self.server.server_port}",f"http://localhost:{self.server.server_port}"}:
                self.send(403,{"error":"Cross-origin requests are not allowed"}); return
            try:
                size = int(self.headers.get("Content-Length","0"))
                if not 0 < size <= 40000: raise ValueError("Invalid request size")
                body = json.loads(self.rfile.read(size))
                if not isinstance(body,dict): raise ValueError("Expected JSON object")
                if self.path == "/api/session": result={"session":agent.new_session()}
                elif self.path == "/api/import": result={"task":agent.attach(body["session"],body["plan"])}
                elif self.path == "/api/message":
                    result={"message":agent.ask(body["session"],body["text"],body.get("provider","deepseek"))}
                else: self.send(404,{"error":"Not found"}); return
                self.send(200,result)
            except (ValueError,KeyError,RuntimeError) as exc: self.send(400,{"error":str(exc)})
            except Exception: self.send(500,{"error":"Request failed. Check model configuration, credentials and server environment; no task success is implied."})
    server = ThreadingHTTPServer(("127.0.0.1",port),Handler)
    server.daemon_threads=True
    server.access_token=secret
    return server


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--storage-root",type=Path,default=Path("/mnt/local/hand/yuzhang/aidd/prompt-workspace"))
    p.add_argument("--runtime",type=Path,default=os.environ.get("AIDD_RUNTIME_PROFILE"))
    p.add_argument("--llm-config-dir",type=Path,default=os.environ.get("AIDD_LLM_CONFIG_DIR"))
    p.add_argument("--port",type=int,default=8765)
    p.add_argument("--allow-compute",action="store_true")
    args=p.parse_args()
    if not 0 < args.port < 65536: p.error("Invalid port")
    lock=args.storage_root.resolve()/"chat/server.lock"
    lock.parent.mkdir(parents=True,exist_ok=True)
    with file_lock(lock):
        agent=ChatAgent(args.storage_root,args.runtime,args.llm_config_dir,args.allow_compute)
        server=None
        try:
            server=make_server(agent,args.port)
            print(f"AIDD chat: http://127.0.0.1:{args.port}/#token={server.access_token}",flush=True)
            print("Keep this process running. Closing the browser does not stop tasks.",flush=True)
            server.serve_forever(poll_interval=.5)
        except KeyboardInterrupt: pass
        finally:
            if server: server.server_close()
            agent.close()
    return 0


if __name__ == "__main__": raise SystemExit(main())
