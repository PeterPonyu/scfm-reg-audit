#!/usr/bin/env python3
"""Post hoc extension of the genuine-brain component streams; immutable R1 inputs."""
import argparse
import copy
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
from scipy.stats import rankdata
from screg_validation.batched_null import batched_components
from screg_validation.common import digest, load_json, sha256, write_json
from screg_validation.statistics import adjust, partial_spearman


def inputs(root):
    panel = load_json(root / 'panels.json')['brain_selected_panel']
    genes = panel['genes']
    tf = np.array([genes.index(g) for g in panel['tf_genes']])
    a = dict(np.load(root / 'graphs.npz', allow_pickle=False))
    ii, jj = np.repeat(tf, len(genes)), np.tile(np.arange(len(genes)), len(tf))
    ii, jj = ii[ii != jj], jj[ii != jj]
    structural = np.column_stack([a[k][jj] for k in ['peak_count', 'gene_length', 'detection', 'gc']])
    controls = np.column_stack([structural, rankdata(a['coexpression'][ii, jj])])
    names = ['geneformer_embed', 'random_init_geneformer', 'scgpt_encoder', 'uce_encoder']
    scores = np.stack([np.load(root / (name + '.npz'), allow_pickle=False)['graph'] for name in names])
    return a, tf, ii, jj, structural, controls, names, scores


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--pilot', action='store_true')
    args = parser.parse_args()
    root, out = args.root, args.out
    out.mkdir(parents=True, exist_ok=True)
    plan = load_json(root / 'posthoc-plan.json')
    for name, expected in plan['input_sha256'].items():
        assert sha256(root / name) == expected, name
    a, tf, ii, jj, structural, controls, names, scores = inputs(root)
    old = load_json(root / 'results.json')
    assert len(old['rows']) == 72 and sum(r['measured'] for r in old['rows']) == 16
    started = time.monotonic()
    verification = []
    runtime = {'host': platform.node(), 'python': platform.python_version(), 'numpy': np.__version__,
               'pid': os.getpid(), 'backend': 'cuda', 'float64': True, 'tf32': False}
    for degree in ['own', 'none']:
        for null in ['restricted', 'unrestricted']:
            key = degree + '-' + null
            prefix = load_json(root / ('prefix-' + key + '.json'))
            assert prefix['next_draw'] == 999
            assert prefix['state_sha256'] == digest({k:v for k,v in prefix.items() if k != 'state_sha256'})
            strata = np.zeros(len(a['proxy']), int)
            if null == 'restricted':
                ob = np.digitize(a['peak_count'], np.unique(np.quantile(a['peak_count'], [.2,.4,.6,.8])))
                gb = np.digitize(a['gc'], np.unique(np.quantile(a['gc'], [1/3,2/3])))
                strata = 15 * np.isin(np.arange(len(a['proxy'])), tf) + 3 * ob + gb
            kw = dict(strata=strata, seed=20261001, degree_mode=degree)
            if args.pilot:
                checks = {}
                for backend in ['cpu', 'cuda']:
                    rows = batched_components(scores, a['proxy'], tf, controls, permutations=8, backend=backend, **kw)
                    error = max(float(np.max(np.abs(np.array([r['nulls'][k] for r in rows]).T - np.array(prefix['nulls'][k])[:8]))) for k in ['M','D'])
                    assert error < 1e-10, (key, backend, error)
                    oldrows = {r['readout']:r for r in old['rows'] if r['measured'] and r['degree_mode']==degree and r['null_mode']==null}
                    obs_error = max(abs(r['observed'] - oldrows[n]['observed']) for n,r in zip(names,rows))
                    assert obs_error < 1e-10
                    # Verify all saved tail counts against the published component p-values.
                    for i, name in enumerate(names):
                        for component in ['M','D']:
                            count = int((np.abs(np.array(prefix['nulls'][component])[:,i]) >= abs(oldrows[name]['observed'])).sum())
                            assert (count+1)/1000 == oldrows[name]['p'+component]
                    checks[backend] = {'prefix_first8_max_abs':error, 'observed_max_abs':obs_error}
                c = structural if degree == 'none' else np.column_stack([structural,(a['proxy']>0).sum(1)[ii],(a['proxy']>0).sum(0)[jj]])
                c = np.column_stack([c, rankdata(a['coexpression'][ii,jj])])
                independent_error = max(abs(partial_spearman(s[ii,jj],a['proxy'][ii,jj],c)-oldrows[n]['observed']) for n,s in zip(names,scores))
                assert independent_error < 1e-10
                verification.append({'group':key,'checks':checks,'direct_projection_max_abs':independent_error})
                print('PILOT',key,verification[-1],flush=True)
                continue
            assert (out/'pilot.json').exists(), 'Numerical portability pilot must pass first'
            pilot=load_json(out/'pilot.json')
            assert pilot['passed'] and pilot['plan_sha256']==digest(plan)
            context = {'plan_sha256':digest(plan),'group':key,'readouts':names,'prefix_sha256':sha256(root/('prefix-'+key+'.json'))}
            target=out/key
            target.mkdir(exist_ok=True)
            ck=target/'draw-progress.json'
            if not ck.exists():
                state=copy.deepcopy(prefix)
                state['context_sha256']=digest(context)
                state.pop('state_sha256')
                write_json(ck,{**state,'state_sha256':digest(state)})
            last=[0.]
            def check():
                elapsed=time.monotonic()-started
                if elapsed > 28800: raise RuntimeError('Eight-hour fixed compute budget reached')
                if elapsed-last[0]>60:
                    saved=load_json(ck)['next_draw']
                    print('PROGRESS',key,saved,'of 9999; elapsed',round(elapsed,1),flush=True)
                    last[0]=elapsed
            if (target/'rows.json').exists():
                print('COMPLETE previously',key,flush=True)
                continue
            rows=batched_components(scores,a['proxy'],tf,controls,permutations=9999,backend='cuda',
                                    check=check,checkpoint=ck,checkpoint_context=context,**kw)
            saved=load_json(ck)
            assert saved['next_draw']==9999
            for k in ['M','D']: assert saved['nulls'][k][:999]==prefix['nulls'][k]
            for name,row in zip(names,rows):
                row.pop('nulls')
                row.update(readout=name,null_mode=null)
            write_json(target/'rows.json',rows)
            print('FINISHED',key,'min IUT',min(r['pIUT'] for r in rows),flush=True)
    if args.pilot:
        write_json(out/'pilot.json',{'passed':True,'plan_sha256':digest(plan),'runtime':runtime,'groups':verification,'seconds':time.monotonic()-started})
        return
    rows=copy.deepcopy(old['rows'])
    for row in rows:
        if row['measured']:
            newer=next(r for r in load_json(out/(row['degree_mode']+'-'+row['null_mode'])/'rows.json') if r['readout']==row['readout'])
            assert abs(row['observed']-newer['observed']) < 1e-10
            row['registered_999']={k:row[k] for k in ['pM','pD','pIUT','q_IUT_BY','p_IUT_Holm','observed']}
            row.update(newer,correction_input_p=newer['pIUT'])
    for method,field in [('BY','q_IUT_BY'),('Holm','p_IUT_Holm')]:
        for row,value in zip(rows,adjust([r['correction_input_p'] for r in rows],method)):row[field]=float(value)
        measured=[r for r in rows if r['measured']]
        for row,value in zip(measured,adjust([r['pIUT'] for r in measured],method)):row[field+'_measured16']=float(value)
    report={'status':'complete','post_hoc':True,'plan_sha256':digest(plan),'runtime':runtime,
            'permutations':9999,'registered_hypotheses':72,'measured_hypotheses':16,'rows':rows,
            'scope':'Resolution sensitivity; original 999-draw analysis retained. Measured-only family is exploratory.',
            'seconds':time.monotonic()-started}
    write_json(out/'results.json',report)
    write_json(out/'completion.json',{'plan_sha256':digest(plan),'results_sha256':sha256(out/'results.json'),
               'outputs':{str(p.relative_to(out)):sha256(p) for p in out.rglob('*.json') if p.name!='completion.json'},
               'prefix_preserved':True,'all_four_groups_complete':True})
    print('ALL COMPLETE', report['seconds'],flush=True)


if __name__=='__main__': main()
