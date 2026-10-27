#!/usr/bin/env python3
"""Fetch publisher supplement URLs explicitly supplied by the user; validate content type and hash."""
import argparse,hashlib,json,urllib.request
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--url',required=True);p.add_argument('--output',required=True);p.add_argument('--expect',choices=['xlsx','pdf','zip'],required=True);a=p.parse_args()
req=urllib.request.Request(a.url,headers={'User-Agent':'Mozilla/5.0 research supplementary reproduction'})
with urllib.request.urlopen(req,timeout=90) as response:blob=response.read(150_000_000)
magic={'xlsx':b'PK\x03\x04','zip':b'PK\x03\x04','pdf':b'%PDF-'}
if not blob.startswith(magic[a.expect]):raise RuntimeError('Wrong file signature: publisher may have returned an HTML/login page; file NOT accepted')
path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(blob)
print(json.dumps({'file':str(path),'size':len(blob),'sha256':hashlib.sha256(blob).hexdigest(),'source_url':a.url},indent=2))
