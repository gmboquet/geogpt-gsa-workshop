"""Participant evidence inspection and non-overwriting experiment records.

The numerical methods remain in workshop.py. No language-model text is generated.
"""
from pathlib import Path
import hashlib, io, json, re, time, threading
import pandas as pd
import numpy as np
import workshop
ROOT = Path(__file__).resolve().parent
RUNS = ROOT / 'participant-runs'
_RUN_LOCK = threading.RLock()

def evidence(case):
    if case == 'fossils':
        raw=(ROOT/'data/raw/fossils/M0027A-forams.tsv').read_text()
        df=pd.read_csv(io.StringIO(raw.split('*/',1)[1].strip()),sep='\t')
        cols=[c for c in df if c in ['Sample label','Depth sed [m]','Preserv','Foram plankt [#]','P. sicana','C. dissimilis']]
        return df[cols].copy()
    if case == 'resolution':
        _,_,_,_,_,logs,_=workshop.load_seismic()
        return pd.DataFrame({'depth_m':logs[:,0],'density_kg_m3':logs[:,1],
            'sonic_us_m':logs[:,2],'velocity_m_s':1e6/logs[:,2],
            'impedance_kg_m2_s':logs[:,1]*1e6/logs[:,2]})
    if case == 'spatial':
        return pd.read_csv(ROOT/'data/prepared/gorkha-grid.csv')
    raise ValueError('Choose fossils, resolution or spatial')

def run(case, label, question, prediction, **settings):
    """Serialize shared plotting/output state; record success after all artifacts exist."""
    with _RUN_LOCK:
        return _run(case, label, question, prediction, **settings)


def _run(case, label, question, prediction, **settings):
    if not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,59}', label):
        raise ValueError('Use a short label containing letters, digits, hyphens or underscores.')
    if not question.strip() or not prediction.strip():
        raise ValueError('Write your question and prediction before running the experiment.')
    fn = {'fossils': workshop.fossils, 'resolution': workshop.resolution, 'spatial': workshop.spatial}[case]
    json.dumps(settings, allow_nan=False)
    folder = RUNS / case / label
    folder.mkdir(parents=True, exist_ok=False)
    record = dict(case=case, label=label, question=question,
                  prediction_before_execution=prediction, settings=settings, status='running',
                  started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                  method_sha256=hashlib.sha256((ROOT / 'workshop.py').read_bytes()).hexdigest(),
                  runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (folder / 'experiment-record.json').write_text(json.dumps(record, indent=2, allow_nan=False))
    old = workshop.OUT
    workshop.OUT = folder
    start = time.perf_counter()
    try:
        result = fn(**settings)
        if case == 'spatial':
            table = pd.DataFrame(result['metrics'])
            table['skill_vs_constant'] = 1 - table.brier / table.baseline_brier
            table.to_csv(folder / 'skill-comparison.csv', index=False)
            pred = pd.read_csv(folder / '03-holdout-predictions.csv')
            rows = []
            for col in ['pred_terrain', 'pred_terrain_shaking']:
                bins = pd.cut(pred[col], bins=np.linspace(0, 1, 6), include_lowest=True)
                for interval, g in pred.groupby(bins, observed=True):
                    rows.append(dict(model=col, probability_bin=str(interval), n=len(g),
                                     mean_prediction=g[col].mean(), observed_fraction=g.landslide_present.mean()))
            pd.DataFrame(rows).to_csv(folder / 'reliability-bins.csv', index=False)
        else:
            table = pd.DataFrame([{k: v for k, v in result.items()
                                  if not isinstance(v, (dict, list)) and k not in ['conclusion', 'limitation', 'seconds']}])
        # Refuse non-finite result metadata rather than advertise a plausible completed run.
        json.dumps(result, allow_nan=False)
        record.update(status='completed', result=result)
    except Exception as exc:
        record.update(status='failed', error=str(exc))
        raise
    finally:
        workshop.OUT = old
        record['wall_seconds'] = time.perf_counter() - start
        (folder / 'experiment-record.json').write_text(json.dumps(record, indent=2, allow_nan=False))
    print('Saved experiment:', folder.relative_to(ROOT))
    return table

def compare(case, labels=None):
    rows=[]
    selected=None if labels is None else set(labels)
    found=set()
    for p in sorted((RUNS/case).glob('*/experiment-record.json')):
        r=json.loads(p.read_text())
        if selected is not None and r['label'] not in selected:continue
        found.add(r['label'])
        if r.get('status')!='completed':continue
        result=r['result']
        if case=='spatial':
            for m in result['metrics']:
                if m['split']=='geographic':
                    rows.append({'label':r['label'],**r['settings'],**m,'skill_vs_constant':1-m['brier']/m['baseline_brier']})
        else:
            rows.append({'label':r['label'],**{k:v for k,v in result.items() if k not in ['conclusion','limitation','seconds']}})
    if selected is not None and found != selected:
        raise ValueError('Unknown run labels: '+', '.join(sorted(selected-found)))
    return pd.DataFrame(rows)

def export_memo(case, question, evidence_used, decision, alternative, limitation, next_test, label='research-memo'):
    if case not in ['fossils','resolution','spatial'] or not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_-]{0,59}',label):
        raise ValueError('Choose a valid case and a unique short memo label.')
    fields=locals().copy();fields.pop('case');fields.pop('label')
    if not all(str(v).strip() for v in fields.values()):
        raise ValueError('Complete every memo field, including the surviving alternative and next test.')
    folder=RUNS/case;folder.mkdir(parents=True,exist_ok=True)
    path=folder/(label+'.md')
    text=('# Research decision — '+case+'\n\n'+'\n\n'.join('## '+k.replace('_',' ').title()+'\n\n'+str(v) for k,v in fields.items()))
    with path.open('x') as f:f.write(text)
    print('Saved',path.relative_to(ROOT));return path
