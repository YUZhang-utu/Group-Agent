"""Open an isolated Chat UI fixture without credentials or scientific execution."""
import argparse
import json
from pathlib import Path
import tempfile
import time

from aidd_agent import domain_agent
from aidd_agent.chat_agent import ChatAgent
from aidd_agent.chat_web import make_server


def fixture_planner(prompt, profile, **kwargs):
    context = json.loads(prompt)
    if not context['evidence']:
        return dict(kind='tool', tool='tasks', arguments={}, purpose='Inspecting this isolated conversation task list.'), {}
    answer = ('## Interface test completed\n\n**Synthetic model response.** The real local task-list tool was called. '
              'No scientific calculation was started.\n\n'
              '| Check | Result |\n| --- | --- |\n| Tool connection | Local task list returned |\n'
              '| Streaming | Answer arrives in chunks |\n\n'
              '- Files remain in the temporary Agent Project.\n- Reload to check conversation persistence.\n\n'
              '```text\nThis code block can be copied.\n```\n\n')
    if 'Attached local file references' in context['request']:
        answer += 'Attachment references reached the planner; raw file bodies were not inserted.\n'
    answer += 'This fixture does not validate a live provider or GPU runtime.'
    value = dict(kind='answer', answer=answer, evidence_ids=[context['evidence'][0]['id']])
    raw = json.dumps(value)
    for offset in range(0, len(raw), 16):
        if kwargs.get('on_content'):
            kwargs['on_content'](raw[offset:offset+16])
        time.sleep(.12)
    return value, {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=0)
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix='aidd-chat-feature-qa-'))
    app = ChatAgent(root, start=False)
    domain_agent.chat_plan = fixture_planner
    server = make_server(app, args.port, token='isolated-feature-fixture')
    print('ISOLATED FIXTURE: no live provider or scientific execution', flush=True)
    print('Storage: '+str(root), flush=True)
    print(f'http://127.0.0.1:{server.server_port}/#token=isolated-feature-fixture', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        app.close()


if __name__ == '__main__':
    main()
