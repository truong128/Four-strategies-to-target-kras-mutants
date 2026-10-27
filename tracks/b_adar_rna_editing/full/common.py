"""Shared filesystem, execution, and evidence/traceability utilities."""
import csv, hashlib, json, os, re, shutil, subprocess, time
from pathlib import Path

class GateError(RuntimeError):
    pass

def rows(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def csvout(path, data, fields=None):
    p=Path(path);p.parent.mkdir(parents=True, exist_ok=True)
    data=list(data)
    cols=fields or (list(data[0]) if data else [])
    if not cols:raise ValueError('CSV columns required')
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(data)
    return str(p)

def jsonout(path, obj):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,sort_keys=True,default=str)+'\n')
    return str(p)

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def which(tool):
    exe=shutil.which(tool)
    if not exe:raise GateError(f'Missing executable: {tool}. Install tool or use --plan-only.')
    return exe

def run(command, log, execute=False, cwd=None):
    log=Path(log);log.parent.mkdir(parents=True,exist_ok=True)
    jsonout(str(log)+'.command.json',{'argv':list(map(str,command)),'cwd':str(cwd) if cwd else None,'execute':execute})
    if not execute:return 'PLAN_ONLY'
    with log.open('w') as fp:
        try:
            result=subprocess.run([str(x) for x in command],cwd=cwd,stdout=fp,stderr=subprocess.STDOUT,check=False)
        except FileNotFoundError as exc:
            raise GateError(f'Executable not found: {command[0]}') from exc
    if result.returncode!=0:raise GateError(f'Command failed ({result.returncode}), see {log}')
    return 'EXECUTED'

def path_required(value,label):
    if not value:raise GateError(f'{label} is required')
    p=Path(value).expanduser().resolve()
    if not p.is_file():raise GateError(f'{label} missing: {p}')
    return p

def bounded_id(raw):
    x=re.sub(r'[^A-Za-z0-9._-]','_',str(raw))
    if not x or x.startswith('.'):raise GateError('Unsafe ID')
    return x[:130]
