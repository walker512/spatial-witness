"""Bounded, exact-threshold point-sample comparison. No registration or surface claim."""
import csv, hashlib, io, json, os, re, stat
from pathlib import Path
INPUT_LIMIT=2*1024*1024
OUTPUT_LIMIT=16*1024*1024
class Invalid(ValueError):pass

def read(path,limit=INPUT_LIMIT):
    try:
        if not stat.S_ISREG(Path(path).lstat().st_mode):raise Invalid('Regular non-symlink file required.')
        fd=os.open(path,os.O_RDONLY|getattr(os,'O_NONBLOCK',0)|getattr(os,'O_NOFOLLOW',0)|getattr(os,'O_BINARY',0))
        try:
            info=os.fstat(fd)
            if not stat.S_ISREG(info.st_mode) or info.st_size>limit:raise Invalid('File type or size limit exceeded.')
            with os.fdopen(fd,'rb') as f:fd=None;raw=f.read(limit+1)
        finally:
            if fd is not None:os.close(fd)
        if len(raw)>limit:raise Invalid('File size limit exceeded.')
        raw.decode('utf8');return raw
    except (OSError,UnicodeError):raise Invalid('Cannot safely read bounded UTF-8 file.') from None

def canonical(v):return (json.dumps(v,sort_keys=True,ensure_ascii=True,indent=2,allow_nan=False)+'\n').encode('ascii')
def digest(raw):return hashlib.sha256(raw).hexdigest()
def fields(d,names):
    if type(d) is not dict or set(d)!=set(names):raise Invalid('Unexpected or missing fields.')
def integer(x,lo,hi):
    if type(x) is not int or not lo<=x<=hi:raise Invalid('Integer outside allowed bounds.')
def parse(raw):
    def pairs(items):
        out={}
        for k,v in items:
            if k in out:raise Invalid('Duplicate JSON keys.')
            k.encode('utf8');out[k]=v
        return out
    def bad(x):raise Invalid('Non-finite JSON.')
    try:
        if len(raw)>INPUT_LIMIT:raise Invalid('Input too large.')
        d=json.loads(raw.decode('utf8'),object_pairs_hook=pairs,parse_constant=bad)
        todo=[(d,0)]
        while todo:
            v,n=todo.pop()
            if isinstance(v,(dict,list)):
                if n>=16:raise Invalid('JSON too deep.')
                todo.extend((x,n+1) for x in (v.values() if isinstance(v,dict) else v))
            elif isinstance(v,str):v.encode('utf8')
            elif isinstance(v,float):raise Invalid('Only integer contract values supported.')
        return d
    except (ValueError,UnicodeError,RecursionError):raise Invalid('Invalid bounded JSON.') from None

def contract(raw):
    d=parse(raw);fields(d,['schema','frame','reference_unit','reconstruction_unit','threshold_um','minimum_precision_bp','minimum_recall_bp','regions'])
    integer(d['schema'],1,1)
    if type(d['frame']) is not str or not re.fullmatch('[A-Za-z][A-Za-z0-9_.-]{0,47}',d['frame']):raise Invalid('Explicit common frame identifier required.')
    for k in ('reference_unit','reconstruction_unit'):
        if type(d[k]) is not str or d[k] not in ('m','mm'):raise Invalid('Unit must be m or mm.')
    integer(d['threshold_um'],1,10**12)
    for k in ('minimum_precision_bp','minimum_recall_bp'):integer(d[k],0,10000)
    if type(d['regions']) is not list or len(d['regions'])>16:raise Invalid('At most 16 regions.')
    ids=set()
    for r in d['regions']:
        fields(r,['id','min_um','max_um','minimum_recall_bp'])
        if type(r['id']) is not str or not re.fullmatch('[A-Za-z][A-Za-z0-9_.-]{0,47}',r['id']) or r['id'] in ids:raise Invalid('Invalid or duplicate region ID.')
        ids.add(r['id']);integer(r['minimum_recall_bp'],0,10000)
        for k in ('min_um','max_um'):
            if type(r[k]) is not list or len(r[k])!=3:raise Invalid('Region requires three coordinates.')
            for x in r[k]:integer(x,-10**12,10**12)
        if any(a>b for a,b in zip(r['min_um'],r['max_um'])):raise Invalid('Inverted region bounds.')
    return d

def coordinate(text,unit):
    precision=6 if unit=='m' else 3
    if not re.fullmatch(r'[+-]?\d{1,13}(?:\.\d{1,'+str(precision)+r'})?',text):raise Invalid('Expected bounded fixed decimal coordinate.')
    negative=text.startswith('-');pieces=text.lstrip('+-').split('.')
    value=int(pieces[0])*10**precision+(int(pieces[1].ljust(precision,'0')) if len(pieces)>1 else 0)
    if negative:value=-value
    if abs(value)>10**12:raise Invalid('Coordinate exceeds one million metres.')
    return value

def cloud(raw,unit):
    if unit not in ('m','mm'):raise Invalid('Unsupported unit.')
    try:
        if len(raw)>INPUT_LIMIT:raise Invalid('Cloud input too large.')
        rows=csv.reader(io.StringIO(raw.decode('utf8'),newline=''),strict=True)
        if next(rows,None)!=['x','y','z']:raise Invalid('CSV header must be exactly x,y,z.')
        result=[]
        for row in rows:
            if len(row)!=3:raise Invalid('Every CSV row must have three coordinates.')
            result.append(tuple(coordinate(x,unit) for x in row))
            if len(result)>1000:raise Invalid('At most 1000 points per cloud.')
        return result
    except (csv.Error,UnicodeError):raise Invalid('Invalid point CSV.') from None

def summary(points):
    return {'count':len(points),'duplicate_rows':len(points)-len(set(points)), 'bounds_um':None if not points else {'min':[min(p[i] for p in points) for i in range(3)],'max':[max(p[i] for p in points) for i in range(3)]}}

def nearest(left,right):
    out=[]
    for p in left:
        best=None;index=None
        for i,q in enumerate(right):
            sq=sum((a-b)**2 for a,b in zip(p,q))
            if best is None or sq<best:best,index=sq,i
        out.append({'nearest_index':index,'distance_squared_um2':best})
    return out

def metric(rows,threshold):
    n=len(rows);hit=sum(r['distance_squared_um2'] is not None and r['distance_squared_um2']<threshold**2 for r in rows)
    return {'matched':hit,'total':n,'ratio':None if not n else hit/n}
def passes(m,bp):return m['total']>0 and m['matched']*10000>=bp*m['total']

def evaluate(craw,rraw,praw):
    c=contract(craw);r=cloud(rraw,c['reference_unit']);p=cloud(praw,c['reconstruction_unit'])
    if not r:raise Invalid('Reference must contain points; no geometric quality claim without a reference.')
    forward=nearest(p,r);backward=nearest(r,p)
    precision=metric(forward,c['threshold_um']);recall=metric(backward,c['threshold_um'])
    regions=[]
    for reg in c['regions']:
        indices=[i for i,v in enumerate(r) if all(a<=x<=b for x,a,b in zip(v,reg['min_um'],reg['max_um']))]
        m=metric([backward[i] for i in indices],c['threshold_um'])
        regions.append({'id':reg['id'],'reference_indices':indices,'recall':m,'minimum_recall_bp':reg['minimum_recall_bp'],'status':'unknown_no_reference_samples' if not indices else ('pass' if passes(m,reg['minimum_recall_bp']) else 'fail')})
    ok=passes(precision,c['minimum_precision_bp']) and passes(recall,c['minimum_recall_bp']) and all(z['status']=='pass' for z in regions)
    # Exact rational F-score converted only for presentation; gates use integer hit counts.
    den=precision['matched']*recall['total']+recall['matched']*precision['total']
    fscore=None if not p else (0.0 if not den else 2*precision['matched']*recall['matched']/den)
    return {'schema':1,'input_sha256':{'contract':digest(craw),'reference':digest(rraw),'reconstruction':digest(praw)},'contract':c,'reference':summary(r),'reconstruction':summary(p),'precision':precision,'recall':recall,'fscore':fscore,'regions':regions,'status':'pass' if ok else 'fail_or_unknown','reconstruction_to_reference':forward,'reference_to_reconstruction':backward}

def outputs(report):
    f=io.StringIO(newline='');w=csv.writer(f,lineterminator='\n');w.writerow(['direction','point_index','nearest_index','distance_squared_um2','within_strict_threshold'])
    t=report['contract']['threshold_um']**2
    for key in ('reconstruction_to_reference','reference_to_reconstruction'):
        for i,r in enumerate(report[key]):
            sq=r['distance_squared_um2'];w.writerow([key,i,r['nearest_index'],sq,sq is not None and sq<t])
    lines=['# Spatial Witness','', 'Status: '+report['status'], '', 'Point-sample comparison only. Alignment, units and reference truth are user declarations.', 'Strict distance < threshold; gates use exact integer arithmetic. Duplicates retain sample weight.', '', f"Common declared frame: {report['contract']['frame']}",f"Threshold: {report['contract']['threshold_um']} micrometres",'']
    for k in ('precision','recall'):
        m=report[k];lines.append(f"{k}: {m['matched']}/{m['total']}"+(' (undefined)' if not m['total'] else ''))
    lines+=['','Regions (reference samples inside inclusive boxes; matched against the entire reconstruction):']
    for r in report['regions']:lines.append(f"- {r['id']}: {r['status']}, {r['recall']['matched']}/{r['recall']['total']}")
    lines+=['','This bundle contains the full source coordinates. It is not anonymous, a surface certificate, or a benchmark submission.','']
    return {'report.json':canonical(report),'distances.csv':f.getvalue().encode('utf8'),'review.md':'\n'.join(lines).encode('utf8')}

def build(config,reference,reconstruction,out):
    raws=[read(x) for x in (config,reference,reconstruction)];report=evaluate(*raws)
    data=dict(zip(('contract.json','reference.csv','reconstruction.csv'),raws));data.update(outputs(report))
    if any(len(raw)>OUTPUT_LIMIT for raw in data.values()):raise Invalid('Generated report too large.')
    Path(out).mkdir(parents=False,exist_ok=False)
    for name,raw in data.items():(Path(out)/name).write_bytes(raw)
    return report

def verify(path):
    p=Path(path);report=evaluate(read(p/'contract.json'),read(p/'reference.csv'),read(p/'reconstruction.csv'))
    for name,raw in outputs(report).items():
        if read(p/name,OUTPUT_LIMIT)!=raw:raise Invalid('Bundle output does not match source recomputation.')
    return report
