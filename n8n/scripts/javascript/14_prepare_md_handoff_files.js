/** n8n Code node: Prepare MD Handoff Files. */
// n8n supplies globals such as $json, $items and $binary at runtime.

async function run() {
  const bundle = $('Build MD Handoff Bundle').first().json;
  const folder = $('Create MD Handoff Folder').first().json;
  if (!folder.id) throw new Error('Google Drive MD_Handoff folder was not created');
  if (!bundle.md_files_b64 || !Object.keys(bundle.md_files_b64).length) throw new Error('MD handoff bundle is empty');
  const mime = (name) => name.endsWith('.json') ? 'application/json'
    : name.endsWith('.html') ? 'text/html'
    : name.endsWith('.svg') ? 'image/svg+xml'
    : name.endsWith('.md') || name.endsWith('.log') || name.endsWith('.pdb') || name.endsWith('.pdbqt') || name.endsWith('.cif') ? 'text/plain'
    : 'application/octet-stream';
  return Object.entries(bundle.md_files_b64).map(([fileName, fileB64]) => ({
    json: { fileName, folderId: folder.id, fileB64, mimeType: mime(fileName) }
  }));
}

module.exports = { run };
