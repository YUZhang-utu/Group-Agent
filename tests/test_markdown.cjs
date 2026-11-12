const test = require('node:test');
const assert = require('node:assert/strict');
const {render} = require('../src/aidd_agent/web/markdown.js');
test('headings, emphasis, lists and tables are readable', () => {
  const html = render('# Results\n\n**Complete** and `CL`\n\n- one\n- two\n\n| Molecule | Score |\n| --- | --- |\n| A | 1 |');
  for (const token of ['<h2>Results</h2>','<strong>Complete</strong>','<code>CL</code>','<ul>','<th>Molecule</th>','<td>1</td>']) assert.ok(html.includes(token));
});
test('raw HTML, script URLs and attributes are not executable', () => {
  const html = render('<img src=x onerror=alert(1)>\n\n[bad](javascript:alert) [quote](https://example.test/\"onmouseover=\"x)');
  assert.ok(!html.includes('<img'));
  assert.ok(!html.includes('href="javascript:'));
  assert.ok(html.includes('&lt;img'));
  assert.ok(!html.includes('href="https://example.test/"onmouseover'));
});
test('code keeps its exact text escaped and includes copy control', () => {
  const html = render('```python\nprint("<safe>")\n```');
  assert.ok(html.includes('data-copy-code'));
  assert.ok(html.includes('print(&quot;&lt;safe&gt;&quot;)'));
});
test('unfinished fences remain safe during streaming', () => {
  assert.ok(render('```html\n<script>x').includes('&lt;script&gt;x'));
});
test('source details collapse safely with clickable paper links', () => {
  const html=render('Result.\n\n```sources\n[Paper](https://example.org/paper)\nPage 13\n<img src=x onerror=alert(1)>\n[bad](javascript:alert)\n```');
  assert.ok(html.includes('<details class="literature-sources">'));
  assert.ok(html.includes('<summary>Sources and reading details</summary>'));
  assert.ok(html.includes('href="https://example.org/paper"'));
  assert.ok(!html.includes('<img')&&!html.includes('href="javascript:'));
  assert.ok(!html.includes(' open'));
});
