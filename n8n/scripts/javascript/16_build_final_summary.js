/** n8n Code node: Build Final Summary. */
// n8n supplies globals such as $json, $items and $binary at runtime.

async function run() {
  const science = $('Build Scientific Report Package').first().json;
  const folder = $('Create Run Folder (Drive)').first().json;
  const mdFolder = $('Create MD Handoff Folder').first().json;
  const reportFiles = $('Upload Reports (Drive)').all().map(({ json }) => ({
    id: json.id,
    name: json.name,
    mimeType: json.mimeType,
    url: json.webViewLink || `https://drive.google.com/file/d/${json.id}/view`,
  }));
  const mdFiles = $('Upload MD Handoff Files').all().map(({ json }) => ({
    id: json.id,
    name: json.name,
    mimeType: json.mimeType,
    url: json.webViewLink || `https://drive.google.com/file/d/${json.id}/view`,
  }));
  return [{ json: {
    ...science.report_summary,
    run_dir: science.run_dir,
    input_dir: science.input_dir,
    output_dir: science.output_dir,
    reports_dir: science.reports_dir,
    drive: {
      folderId: folder.id,
      folderName: folder.name,
      folderUrl: `https://drive.google.com/drive/folders/${folder.id}`,
      reportFiles,
    },
    md_handoff: {
      status: 'STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY',
      folderId: mdFolder.id,
      folderName: mdFolder.name,
      folderUrl: `https://drive.google.com/drive/folders/${mdFolder.id}`,
      files: mdFiles,
      chemistryAuthority: 'best_pose_ligand.sdf',
      warning: 'No force-field topology, solvent, ions, equilibration, trajectory, or MD result is included.',
    },
  } }];
}

module.exports = { run };
