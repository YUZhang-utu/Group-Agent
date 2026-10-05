"""Start an isolated, credential-free UI fixture; never use the real Chat storage."""
import json
from pathlib import Path
import tempfile

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.chat_web import make_server

root = Path(tempfile.mkdtemp(prefix='medchem-ui-'))
app = ChatAgent(root, start=False)
server = make_server(app, 0, token='isolated-ui-fixture')
receipt = Path('data/e098-ui/server.json')
receipt.parent.mkdir(parents=True, exist_ok=True)
receipt.write_text(json.dumps(dict(url=f'http://127.0.0.1:{server.server_port}')), encoding='utf-8')
print('Isolated UI fixture ready', flush=True)
try:
    server.serve_forever()
finally:
    server.server_close()
    app.close()
