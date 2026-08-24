/** n8n Code node: Attach Drive Metadata. */
// n8n supplies globals such as $json, $items and $binary at runtime.

async function run() {
  const source = $('SMILES Drop (Drive Trigger)').first().json;
  const extracted = $input.first().json;
  const text = typeof extracted.data === 'string'
    ? extracted.data
    : (typeof extracted.text === 'string' ? extracted.text : String(extracted.data ?? ''));
  if (!text.trim()) throw new Error('The Drive file is empty or is not extractable as plain text');
  return [{ json: {
    data: text,
    source_file_id: source.id,
    source_file_name: source.name || source.fileName || 'drive-input.txt',
    source_mime_type: source.mimeType || null,
    source_modified_time: source.modifiedTime || null,
    source_md5: source.md5Checksum || null,
  } }];
}

module.exports = { run };
