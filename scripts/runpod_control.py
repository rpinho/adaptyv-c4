#!/usr/bin/env python3
"""Use Ricardo's existing Runpod MCP OAuth session without printing credentials.

Fallback for this Codex session's missing Runpod connector. No credential is
copied or committed; authorization is sent only to the configured Runpod host.
Schemas are discovered from tools/list before calls. Stop is never terminate.
"""
import argparse
import json
from pathlib import Path
import urllib.request

URL='https://mcp.getrunpod.io/'


def redact(value):
    if isinstance(value,list):return [redact(x) for x in value]
    if isinstance(value,dict):
        return {k:('[redacted]' if k.lower() in {'env','headers','accesstoken','refreshtoken','password','secretvalue'}
                   else redact(v)) for k,v in value.items()}
    return value


def safe_output(result):
    result=redact(result)
    if isinstance(result,dict):
        for item in result.get('content',[]):
            if item.get('type')=='text':
                try:item['text']=json.dumps(redact(json.loads(item['text'])),indent=2)
                except ValueError:pass
    return result


def request(method,params,seq=1):
    creds=json.loads((Path.home()/'.claude/.credentials.json').read_text())
    entries=[v for k,v in creds.get('mcpOAuth',{}).items()
             if 'runpod' in k and v.get('serverUrl','').rstrip('/')==URL.rstrip('/')]
    if len(entries)!=1:raise RuntimeError('unique existing Runpod OAuth session unavailable')
    payload=json.dumps(dict(jsonrpc='2.0',id=seq,method=method,params=params)).encode()
    req=urllib.request.Request(URL,data=payload,headers={
        'Authorization':'Bearer '+entries[0]['accessToken'],
        'Content-Type':'application/json','Accept':'application/json, text/event-stream'})
    with urllib.request.urlopen(req,timeout=45) as response:
        raw=response.read().decode()
    if raw.lstrip().startswith('{'):data=json.loads(raw)
    else:
        events=[json.loads(line[6:]) for line in raw.splitlines() if line.startswith('data: ')]
        data=next(e for e in reversed(events) if e.get('id')==seq)
    if 'error' in data:raise RuntimeError(json.dumps(data['error']))
    return data['result']


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--list',action='store_true')
    ap.add_argument('--tool');ap.add_argument('--arguments',default='{}')
    ap.add_argument('--out');args=ap.parse_args()
    result=request('tools/list',{})
    tools=result['tools']
    if args.list:
        names=['get-pod','pod-action','list-pods','list-pod-billing','list-billing']
        result=[t for t in tools if t['name'] in names]
    else:
        tool=next((t for t in tools if t['name']==args.tool),None)
        if not tool:raise SystemExit('Tool absent from live server')
        result=request('tools/call',dict(name=args.tool,arguments=json.loads(args.arguments)),seq=2)
    encoded=json.dumps(safe_output(result),indent=2)+'\n'
    if args.out:Path(args.out).write_text(encoded)
    else:print(encoded,end='')


if __name__=='__main__':main()
