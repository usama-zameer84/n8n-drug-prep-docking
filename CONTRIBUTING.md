# Contributing

Changes are welcome when they keep the workflow reproducible and the scientific limits explicit.

## Before opening a pull request

1. Edit the workflow JSON or the relevant Code-node source.
2. Run `python3 n8n/tools/extract_code_nodes.py` after any Code-node change.
3. Run `python3 -m unittest discover -s tests -v`.
4. Document changes to input fields, ranking logic, report fields, or MD handoff files.
5. Do not commit credentials, tokens, Drive folder IDs, patient data, proprietary compounds, or unpublished structures.

For changes to scoring or receptor selection, include a small public test case and explain the expected scientific effect. A passing workflow run is required, but it does not replace review of the receptor, box, pose, and chemistry.
