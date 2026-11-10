"""Bounded provider SSE decoding and answer-only structured JSON previews."""
import json

LIMIT = 2 * 1024 * 1024


def answer_prefix(content):
    """Return only a top-level answer string, never tool arguments or reasoning."""
    decoder = json.JSONDecoder()
    text = content.lstrip()
    if not text.startswith('{'):
        return ''
    pos, values = 1, {}
    while pos < len(text):
        while pos < len(text) and text[pos] in ' \r\n\t,':
            pos += 1
        try:
            key, end = decoder.raw_decode(text, pos)
        except ValueError:
            break
        if not isinstance(key, str):
            break
        pos = end
        while pos < len(text) and text[pos].isspace():
            pos += 1
        if pos >= len(text) or text[pos] != ':':
            break
        pos += 1
        while pos < len(text) and text[pos].isspace():
            pos += 1
        try:
            value, pos = decoder.raw_decode(text, pos)
            values[key] = value
        except ValueError:
            if key == 'answer' and text[pos:pos+1] == '"':
                # Decode only complete escapes; retain unfinished escapes for the next delta.
                result, i = [], pos + 1
                while i < len(text):
                    char = text[i]
                    if char == '"' or ord(char) < 32:
                        break
                    if char == '\\':
                        size = 6 if text[i+1:i+2] == 'u' else 2
                        if i + size > len(text):
                            break
                        try:
                            char = json.loads('"' + text[i:i+size] + '"')
                        except ValueError:
                            break
                        if len(char) == 1 and 0xD800 <= ord(char) <= 0xDBFF:
                            if i + 12 > len(text):
                                break
                            try:
                                char = json.loads('"' + text[i:i+12] + '"')
                            except ValueError:
                                break
                            size = 12
                        i += size
                    else:
                        i += 1
                    result.append(char)
                values[key] = ''.join(result)
            break
    answer = values.get('answer', '')
    return answer if values.get('kind') == 'answer' and isinstance(answer, str) else ''


def read_sse(response, on_content):
    """Collect Chat Completions SSE while forwarding actual content deltas."""
    total, content, data = 0, [], []
    result, finish = {}, None

    def consume(lines):
        nonlocal finish, result
        payload = '\n'.join(lines)
        if not payload or payload == '[DONE]':
            return
        event = json.loads(payload)
        if 'error' in event:
            raise ValueError('Provider reported a streaming error')
        if event.get('id'):
            result['id'] = event['id']
        if event.get('usage'):
            result['usage'] = event['usage']
        for choice in event.get('choices', []):
            if choice.get('index', 0) != 0:
                continue
            if choice.get('finish_reason'):
                finish = choice['finish_reason']
            delta = choice.get('delta', {}).get('content')
            if delta:
                if not isinstance(delta, str):
                    raise ValueError('Invalid provider content delta')
                content.append(delta)
                on_content(delta)

    while True:
        raw = response.readline(LIMIT + 1)
        if not raw:
            break
        total += len(raw)
        if total > LIMIT:
            raise ValueError('LLM response too large')
        line = raw.decode('utf-8').rstrip('\r\n')
        if not line:
            consume(data)
            data = []
        elif line.startswith('data:'):
            data.append(line[5:].lstrip(' '))
    if data:
        consume(data)
    result['choices'] = [dict(finish_reason=finish, message=dict(content=''.join(content)))]
    return result
