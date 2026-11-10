/* A deliberately small Markdown renderer: raw HTML and remote images stay inert. */
(function (root) {
  'use strict';
  const escape = value => String(value).replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  function link(url) {
    if (/^(https?:\/\/|mailto:)/i.test(url) && !/[\x00-\x20\x7f]/.test(url)) return escape(url);
    return null;
  }
  function inline(text) {
    const pattern = /(`+)([^`]+)\1|\[([^\]\n]+)\]\(([^\s)]+)\)|\*\*([^*\n]+)\*\*|__([^_\n]+)__|\*([^*\n]+)\*|~~([^~\n]+)~~/g;
    let output = '', offset = 0, match;
    while ((match = pattern.exec(text))) {
      output += escape(text.slice(offset, match.index));
      if (match[1]) output += '<code>'+escape(match[2])+'</code>';
      else if (match[3]) {
        const href = link(match[4]);
        output += href ? '<a href="'+href+'" target="_blank" rel="noopener noreferrer">'+escape(match[3])+'</a>' : escape(match[0]);
      } else if (match[5] || match[6]) output += '<strong>'+escape(match[5] || match[6])+'</strong>';
      else if (match[7]) output += '<em>'+escape(match[7])+'</em>';
      else output += '<del>'+escape(match[8])+'</del>';
      offset = pattern.lastIndex;
    }
    return output + escape(text.slice(offset));
  }
  const cells = line => line.trim().replace(/^\|/, '').replace(/\|$/, '').split(/(?<!\\)\|/).map(cell => cell.trim().replace(/\\\|/g, '|'));
  const separator = line => line.includes('|') && cells(line).every(cell => /^:?-{3,}:?$/.test(cell));
  function render(source) {
    const lines = String(source).replace(/\r\n?/g, '\n').split('\n');
    let output = '', i = 0;
    while (i < lines.length) {
      const line = lines[i];
      if (!line.trim()) { i++; continue; }
      const fence = line.match(/^\s*(`{3,}|~{3,})([\w+-]*)\s*$/);
      if (fence) {
        const code = []; i++;
        while (i < lines.length && !(lines[i].trim().length >= fence[1].length && [...lines[i].trim()].every(c => c === fence[1][0]))) code.push(lines[i++]);
        if (i < lines.length) i++;
        output += '<div class="code-block"><div class="code-toolbar"><span>'+escape(fence[2] || 'Code')+'</span><button type="button" data-copy-code>Copy code</button></div><pre><code>'+escape(code.join('\n'))+'</code></pre></div>';
        continue;
      }
      const heading = line.match(/^(#{1,6})\s+(.+)$/);
      if (heading) { const level = Math.min(heading[1].length+1, 6); output += '<h'+level+'>'+inline(heading[2])+'</h'+level+'>'; i++; continue; }
      if (/^\s*(?:---+|\*\*\*+|___+)\s*$/.test(line)) { output += '<hr>'; i++; continue; }
      if (i+1 < lines.length && line.includes('|') && separator(lines[i+1])) {
        const headers = cells(line); i += 2;
        output += '<div class="table-scroll"><table><thead><tr>'+headers.map(c => '<th>'+inline(c)+'</th>').join('')+'</tr></thead><tbody>';
        while (i < lines.length && lines[i].trim() && lines[i].includes('|')) {
          const row = cells(lines[i++]);
          output += '<tr>'+headers.map((_, j) => '<td>'+inline(row[j] || '')+'</td>').join('')+'</tr>';
        }
        output += '</tbody></table></div>'; continue;
      }
      const item = line.match(/^\s*(?:([-+*])|([0-9]+)[.)])\s+(.+)$/);
      if (item) {
        const ordered = !!item[2], tag = ordered ? 'ol' : 'ul';
        output += '<'+tag+(ordered ? ' start="'+Number(item[2])+'"' : '')+'>';
        while (i < lines.length) {
          const next = lines[i].match(/^\s*(?:([-+*])|([0-9]+)[.)])\s+(.+)$/);
          if (!next || !!next[2] !== ordered) break;
          output += '<li>'+inline(next[3])+'</li>'; i++;
        }
        output += '</'+tag+'>'; continue;
      }
      if (/^>\s?/.test(line)) {
        const quote = [];
        while (i < lines.length && /^>\s?/.test(lines[i])) quote.push(lines[i++].replace(/^>\s?/, ''));
        output += '<blockquote>'+quote.map(inline).join('<br>')+'</blockquote>'; continue;
      }
      const paragraph = [line]; i++;
      while (i < lines.length && lines[i].trim() && !/^(#{1,6}\s|\s*`{3,}|\s*~{3,}|>\s?|\s*[-+*]\s|\s*\d+[.)]\s)/.test(lines[i]) && !(i+1 < lines.length && separator(lines[i+1]))) paragraph.push(lines[i++]);
      output += '<p>'+paragraph.map(inline).join('<br>')+'</p>';
    }
    return output;
  }
  root.AiddMarkdown = {render};
  if (typeof module !== 'undefined') module.exports = {render};
})(typeof window !== 'undefined' ? window : globalThis);
