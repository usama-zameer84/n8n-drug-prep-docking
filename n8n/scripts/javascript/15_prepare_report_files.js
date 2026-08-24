/** n8n Code node: Prepare Report Files. */
// n8n supplies globals such as $json, $items and $binary at runtime.

async function run() {
  const summary = $('Build Scientific Report Package').first().json;
  const folder = $('Create Run Folder (Drive)').first().json;
  if (!folder.id) throw new Error('Google Drive run folder was not created');
  if (!summary.reports_b64 || !Object.keys(summary.reports_b64).length) throw new Error('Report package is empty');
  return Object.entries(summary.reports_b64).map(([fileName, fileB64]) => ({
    json: {
      fileName,
      folderId: folder.id,
      fileB64,
      mimeType: fileName.endsWith('.html') ? 'text/html' : fileName.endsWith('.sdf') ? 'chemical/x-mdl-sdfile' : fileName.endsWith('.pdb') ? 'chemical/x-pdb' : 'application/json',
    }
  }));
}

module.exports = { run };
