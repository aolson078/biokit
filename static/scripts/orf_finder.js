document.addEventListener('DOMContentLoaded', () => {
  const runBtn = document.getElementById('run-orf');
  if (!runBtn) return;

  const statusEl = document.getElementById('orf-status');
  const resultsEl = document.getElementById('orf-results');

  runBtn.addEventListener('click', async () => {
    const sequence = document.getElementById('orf-sequence').value;
    const minLength = Number(document.getElementById('orf-min-length').value || 90);
    statusEl.textContent = 'Running analysis...';
    resultsEl.innerHTML = '';

    try {
      const response = await fetch('/api/orf-finder', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence, min_length: minLength })
      });
      const data = await response.json();
      if (!response.ok) throw new Error(data.error || 'Request failed');

      statusEl.textContent = `Found ${data.count} ORFs (showing up to 100).`;
      if (!data.orfs.length) {
        resultsEl.innerHTML = '<p class="text-muted">No ORFs found for this threshold.</p>';
        return;
      }

      const rows = data.orfs.map((orf) => `
        <tr>
          <td>${orf.start}</td><td>${orf.end}</td><td>${orf.strand}</td><td>${orf.frame}</td>
          <td>${orf.length_nt}</td><td>${orf.length_aa}</td><td><code>${orf.protein}</code></td>
        </tr>
      `).join('');

      resultsEl.innerHTML = `
        <table class="table table-striped table-sm">
          <thead><tr><th>Start</th><th>End</th><th>Strand</th><th>Frame</th><th>Length (nt)</th><th>Length (aa)</th><th>Protein</th></tr></thead>
          <tbody>${rows}</tbody>
        </table>`;
    } catch (error) {
      statusEl.textContent = `Error: ${error.message}`;
    }
  });
});
