"""Verify the downloaded workshop files and report the hosted environment.
No scientific experiment or package installation is performed.
"""
from pathlib import Path
import hashlib, importlib, json, platform, sys
root=Path(__file__).resolve().parent
manifest=json.loads((root/'runtime-manifest.json').read_text())
failures=[]
for item in manifest['files']:
    p=root/item['path']
    valid=p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
    print(('PASS' if valid else 'FAIL')+' file '+item['path'])
    if not valid: failures.append(item['path'])
print('Python:',platform.python_version())
for name in ['numpy','pandas','matplotlib','scipy','sklearn']:
    try:
        module=importlib.import_module(name)
        print('Available:',name,getattr(module,'__version__','unknown'))
    except ImportError:
        failures.append(name)
        print('Missing dependency:',name)
if failures:
    print('Setup incomplete:',', '.join(failures));sys.exit(1)
print('Setup verified. No research result has been computed.')
print('Keep classroom.py and workshop.py in this folder. Begin the input audit in the lesson.')
