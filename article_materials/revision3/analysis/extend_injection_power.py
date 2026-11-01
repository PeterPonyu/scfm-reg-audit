#!/usr/bin/env python3
"""Increase fixed-panel injection replicates using explicit isolated inputs.

Float64 CUDA evaluates FWL cross-products; SciPy ranks and NumPy SVD preserve
average ties and rank-deficient designs. This is a synthetic sensitivity study,
not donor uncertainty or a validation of randomization exchangeability.
"""
import argparse
import hashlib
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import torch
from scipy.stats import rankdata, beta

A_GRID = [0., .001, .002, .003, .004, .005, .006, .008, .010, .015]
NEW_REPS = 80
N_PERM = 999
NEW_SEED = [20260927, 5613, 4]


def zscore(v):
    v = np.asarray(v, dtype=np.float64)
    return (v-v.mean())/(v.std()+1e-9)


def basis(x):
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    return u[:, s > np.finfo(float).eps * max(x.shape) * s[0]]


def residual(y, q):
    return y - q @ (q.T @ y)


def health():
    raw = subprocess.check_output(['nvidia-smi', '--query-gpu=temperature.gpu,enforced.power.limit,clocks_event_reasons.hw_thermal_slowdown,clocks_event_reasons.sw_thermal_slowdown', '--format=csv,noheader,nounits'], text=True).strip()
    vals=[x.strip() for x in raw.split(',')]
    if float(vals[0]) >= 82 or float(vals[1]) < 90 or any(x != 'Not Active' for x in vals[2:]):
        raise RuntimeError('GPU readiness changed; stop without enqueueing more work: '+raw)
    return raw


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def make_rows(g, old, out, validation_only):
    ii,jj=g['ii'],g['jj']; G=g['G']
    od=(G>0).sum(1).astype(np.float32); ind=(G>0).sum(0).astype(np.float32)
    fixed=np.column_stack([np.ones(len(ii)),rankdata(g['co']),*[zscore(g[k][jj]) for k in ('pc','gl','det','gc')]])
    degrees=np.column_stack([zscore(od[ii]),zscore(ind[jj])])
    C1=np.column_stack([fixed,degrees]); y=rankdata(G[ii,jj])
    # Reproduce the published planted-signal definition, including its original solve.
    za=zscore(y-C1@np.linalg.solve(C1.T@C1,C1.T@y))
    seeds=list(map(int,old['noise_seed'])); doses=list(old['alpha']); reps=list(old['rep'])
    if not validation_only:
        seq=np.random.SeedSequence(NEW_SEED+[0 if g['tissue']=='brain' else 1])
        added=[int(s.generate_state(1,dtype=np.uint64)[0]) for s in seq.spawn(len(A_GRID)*NEW_REPS)]
        assert not set(seeds)&set(added)
        seeds+=added
        doses += [a for a in A_GRID for _ in range(NEW_REPS)]
        reps += [r+20 for _ in A_GRID for r in range(NEW_REPS)]
    x=np.lib.format.open_memmap(out/'rank_rows.npy',mode='w+',dtype='float64',shape=(len(seeds),len(ii)))
    q=basis(C1); yr=residual(y,q);yr-=yr.mean();yn=np.linalg.norm(yr)
    obs=[]
    for j,(s,a) in enumerate(zip(seeds,doses)):
        v=rankdata(a*za+(1-a)*zscore(np.random.default_rng(s).standard_normal(len(za))))
        x[j]=v
        # Exact reference OLS observable (same as original pcorr implementation).
        xr=residual(v,q); xr-=xr.mean()
        obs.append(float(xr@yr)/(np.linalg.norm(xr)*yn))
        if j%100==0: log(g['tissue'],'ranked rows',j,'/',len(seeds))
    x.flush()
    err=float(np.max(np.abs(np.array(obs[:200])-old['observed'])))
    assert err < 1e-9, ('observed equivalence',err)
    log(g['tissue'],'observed original 200 max error',err)
    return x,fixed,od,ind,np.array(obs),np.array(doses),np.array(reps),np.array(seeds,dtype=np.uint64),err


def null_vectors(g,fixed,od,ind,nt,seed):
    ii,jj=g['ii'],g['jj'];G=g['G'];ng=G.shape[0]
    xx=fixed if nt=='mantel' else np.column_stack([fixed,zscore(od[ii])])
    q=basis(xx)
    seq=np.random.SeedSequence(seed);gp=G.copy()
    tf=np.unique(g['tf']);nonself={t:np.r_[0:t,t+1:ng] for t in tf}
    for k in range(N_PERM):
        repseed=int(seq.spawn(1)[0].generate_state(1,dtype=np.uint64)[0]);rng=np.random.default_rng(repseed)
        if nt=='mantel':
            p=rng.permutation(ng); y=rankdata(G[p[ii],p[jj]])
            ch=np.column_stack([zscore(od[p][ii]),zscore(ind[p][jj])])
        else:
            for t in tf:
                idx=nonself[t];v=gp[t,idx].copy();rng.shuffle(v);gp[t,idx]=v
            y=rankdata(gp[ii,jj]);newind=(gp>0).sum(0).astype(np.float64)
            ch=zscore(newind[jj])[:,None]
        qc=basis(residual(ch,q));yr=residual(residual(y,q),qc);yr-=yr.mean()
        yield np.column_stack([yr/np.linalg.norm(yr),qc])


def run_null(g,x,fixed,od,ind,nt,seed,out):
    health();xx=fixed if nt=='mantel' else np.column_stack([fixed,zscore(od[g['ii']])])
    q=torch.as_tensor(basis(xx).T,device='cuda',dtype=torch.float64)
    # Row-major: each row is one fixed synthetic graph rank vector.
    xr=torch.as_tensor(np.asarray(x),device='cuda',dtype=torch.float64)
    xr.sub_((xr@q.T)@q);xr.sub_(xr.mean(dim=1,keepdim=True));xr.div_(torch.linalg.vector_norm(xr,dim=1,keepdim=True))
    result=np.empty((N_PERM,len(x)),dtype=np.float64)
    batch=[];start=0
    for k,v in enumerate(null_vectors(g,fixed,od,ind,nt,seed)):
        batch.append(v)
        if len(batch)==16 or k==N_PERM-1:
            mat=torch.as_tensor(np.column_stack(batch),device='cuda',dtype=torch.float64)
            cross=(xr@mat).cpu().numpy();col=0
            for j,b in enumerate(batch):
                c=cross[:,col:col+b.shape[1]]
                den=np.sqrt(np.maximum(0,1-np.sum(c[:,1:]**2,axis=1)))
                assert np.all(den>0)
                result[start+j]=c[:,0]/den;col+=b.shape[1]
            start=k+1;batch.clear();del mat
        if (k+1)%100==0:
            log(g['tissue'],nt,k+1,'/',N_PERM,'GPU',health())
    del xr,q;torch.cuda.empty_cache()
    np.save(out/(nt+'_nulls.npy'),result)
    return result


def exact_ci(k,n,alpha=.05):
    return [0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),
            1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))]


def wilson(k,n):
    z=1.959963984540054;p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d
    r=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [max(0.,float(c-r)),min(1.,float(c+r))]


def summarize(obs,doses,pm,pd):
    curves=[]
    for a in A_GRID:
        s=doses==a;row={'alpha_injected':a,'n':int(s.sum()),'mean_rho':float(obs[s].mean()),'sd_rho':float(obs[s].std())}
        for tag,threshold in pm['levels'].items():
            k=int(((pm['p'][s]<=threshold)&(pd[s]<=threshold)).sum())
            row[tag]={'count':k,'proportion':k/int(s.sum()),'wilson95':wilson(k,int(s.sum())),
                      'simultaneous95_per_curve':exact_ci(k,int(s.sum()),.05/len(A_GRID))}
        curves.append(row)
    return curves


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--validation-only',action='store_true');args=ap.parse_args()
    torch.set_num_threads(4);torch.backends.cuda.matmul.allow_tf32=False
    args.out.mkdir(parents=True,exist_ok=True)
    plan={'design':'80 new independent noise replicates per dose added to the archived 20; ten fixed doses; 999 shared randomizations per null/tissue',
          'new_seed_contract':NEW_SEED,'new_reps_per_dose':NEW_REPS,'target_total':100,'alpha_grid':A_GRID,
          'success':'Reproduce all 400 archived observed statistics and 799200 archived null statistics to absolute tolerance 1e-9; exact original p-values; no adaptive sample size.',
          'scope':'Conditional on fixed graph, panel, confounds and shared null streams; neither biological power nor donor variability.',
          'health_at_start':health()}
    (args.out/'design.json').write_text(json.dumps(plan,indent=2)+'\n')
    seedroot=np.random.SeedSequence([20260924,2]);s=[int(x.generate_state(1,dtype=np.uint64)[0]) for x in seedroot.spawn(6)]
    report={'design':plan,'tissues':{},'validation':{},'torch':torch.__version__,'device':torch.cuda.get_device_name(0)}
    for t,ss in [('brain',s[:3]),('pbmc',s[3:])]:
        td=args.out/t;td.mkdir(exist_ok=True)
        g=dict(np.load(args.inputs/(t+'.npz')));g['tissue']=t
        old={nt:dict(np.load(args.inputs/f'power_{t}_{nt}.npz')) for nt in ['mantel','degree']}
        x,fixed,od,ind,obs,doses,reps,seeds,err=make_rows(g,old['mantel'],td,args.validation_only)
        arrays={};valid={'observed_max_abs_error':err,'nulls':{}}
        for nt,seed in [('mantel',ss[1]),('degree',ss[2])]:
            null=run_null(g,x,fixed,od,ind,nt,seed,td)
            maxerr=float(np.max(np.abs(null[:,:200]-old[nt]['nulls'])))
            p=(np.sum(np.abs(null)>=np.abs(obs)[None,:],axis=0)+1)/(N_PERM+1)
            op=(np.sum(np.abs(old[nt]['nulls'])>=np.abs(old[nt]['observed'])[None,:],axis=0)+1)/(N_PERM+1)
            assert maxerr<1e-9,(t,nt,maxerr)
            assert np.array_equal(p[:200],op),(t,nt,'p-values differ')
            valid['nulls'][nt]={'max_abs_error':maxerr,'original_p_values_exact':True,'comparisons':int(old[nt]['nulls'].size)}
            arrays[nt]=p
            log(t,nt,'equivalence PASSED',maxerr)
        np.savez_compressed(td/'results.npz',observed=obs,alpha=doses,rep=reps,noise_seed=seeds,pM=arrays['mantel'],pD=arrays['degree'])
        report['validation'][t]=valid
        levels={'a05':.05,'a05m':.05/(8 if t=='brain' else 5)}
        report['tissues'][t]={'curves':summarize(obs,doses,{'levels':levels,'p':arrays['mantel']},arrays['degree'])}
        report['tissues'][t]['original20']=summarize(obs[:200],doses[:200],{'levels':levels,'p':arrays['mantel'][:200]},arrays['degree'][:200])
        report['tissues'][t]['new80']=summarize(obs[200:],doses[200:],{'levels':levels,'p':arrays['mantel'][200:]},arrays['degree'][200:]) if not args.validation_only else []
        (args.out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
        del x
        (td/'rank_rows.npy').unlink()
    report['status']='PASS';report['health_at_end']=health()
    report['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (args.out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    log('DONE',args.out/'summary.json')

if __name__=='__main__':
    main()
