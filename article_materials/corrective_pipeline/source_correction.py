#!/usr/bin/env python3
"""Separately registered source correction; never overwrite historical evidence."""
from __future__ import annotations
import argparse
import hashlib
import time
from contextlib import contextmanager
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr
from screg_validation.common import ContractError, authorize, cache_receipt, digest, load_json, sha256, write_json
from screg_validation.h5ad import index, read_rows, counts_normalization as _counts_normalization, column, parse_peak
from screg_validation.graphs import read_coordinates, read_meme, locus_links, build_motif_hits, gene_coexpression
from screg_validation.models import require_resources, tokenizable_genes, run_readout, READOUT_RESOURCES
from screg_validation.hardware import require_gpu
from screg_validation.batched_null import batched_components
from screg_validation.statistics import adjust, partial_spearman
from screg_validation.replication import coverage_status, source_signature
from provenance_fingerprint import unique_symbols

HERE = Path(__file__).resolve().parent
# This pinned adapter counts with-replacement token appearances, not cells.
UCE_ADAPTER_SHA256 = '9df70d6f1e161be91cb2579d6a1bdc99e36256317376fefe6c5695fee780c439'
UCE_TOKENS_PER_CELL = 1024


def counts_normalization(values, scale):
    """Detach CSR storage before the frozen normalizer sorts sparse indices.

    SciPy dtype conversion may share indices while copying data. Sorting that
    converted matrix could otherwise permute the source indices without source
    data. A canonical private copy preserves exact source gene/count mapping.
    """
    private = sparse.csr_matrix(values, copy=True)
    private.sum_duplicates()
    private.sort_indices()
    private.eliminate_zeros()
    return _counts_normalization(private, scale)


@contextmanager
def budget_stage(out, task, maximum):
    """Persist cumulative runtime for authorized resumes; stale running state blocks."""
    path=out.parent/(out.name+'-'+task+'-budget.json')
    state=load_json(path) if path.exists() else {'elapsed_seconds':0.,'invocations':0,'running':False}
    if state.get('running') or state['elapsed_seconds']>=maximum:
        raise ContractError('cumulative stage budget exhausted or prior invocation was not finalized')
    start=time.monotonic();prior=float(state['elapsed_seconds'])
    write_json(path,{**state,'running':True,'invocations':state['invocations']+1,'maximum':maximum})
    def check():
        if prior+time.monotonic()-start>=maximum:
            raise ContractError('cumulative corrective stage runtime exceeded registered budget')
    try:
        yield check
    finally:
        write_json(path,{'elapsed_seconds':prior+time.monotonic()-start,
            'invocations':state['invocations']+1,'running':False,'maximum':maximum})


def valid_coexpression(log):
    log=sparse.csr_matrix(log,dtype=float)
    graph=gene_coexpression(log)
    variance=np.asarray(log.multiply(log).mean(0)).ravel()-np.asarray(log.mean(0)).ravel()**2
    missing=np.flatnonzero(variance<=0)
    graph[missing,:]=np.nan;graph[:,missing]=np.nan
    return graph,variance


def identity(plan, resources):
    if plan['resources_sha256'] != digest(resources):
        raise ContractError('corrective input registry differs from frozen plan')
    require_resources(resources, list(resources))
    return {'plan': digest(plan), 'resources': digest(resources), 'code': digest({
        'core': source_signature(), 'source_correction.py': sha256(__file__),
        'provenance_fingerprint.py': sha256(HERE/'provenance_fingerprint.py')})}


def finish(root, name, ctx, paths):
    write_json(root/name, {'context_sha256':digest(ctx), 'outputs':{p:sha256(root/p) for p in paths}})


def receipt(root, name, ctx, outputs=None):
    record=cache_receipt(root/name,ctx)
    if not record or (outputs is not None and set(record['outputs']) != set(outputs)):
        raise ContractError('complete corrective receipt required: '+name)
    return record


def official_metadata(resources, name):
    p=Path(resources[name]['path'])
    record=load_json(str(p)+'.download.json')
    if record['sha256'] != resources[name]['sha256'] or record['bytes'] != resources[name]['bytes']:
        raise ContractError('official metadata download record differs')
    frame=pd.read_csv(p,dtype=str,keep_default_na=False)
    if 'Barcode' not in frame or not len(frame) or frame.Barcode.duplicated().any():
        raise ContractError('official metadata exact unique barcodes required')
    return frame


def exact_metadata(frame,cells):
    if len(cells)!=len(set(cells)) or not set(frame.Barcode)<=set(cells):
        raise ContractError('metadata/cell identity mismatch')
    lookup={c:i for i,c in enumerate(cells)}
    return {c:lookup[c] for c in frame.Barcode}


def balanced_selection(rna,atac,rule,seed):
    valid=lambda frame:frame.loc[~frame['Cell.Type'].str.strip().str.casefold().isin({'','na','nan','none','unknown','null'})]
    rna,atac=valid(rna),valid(atac)
    rs,ats=rule['sample_key_rna'],rule['sample_key_atac']
    samples=sorted(set(rna[rs])&set(atac[ats]))
    if len(samples)<3 or any(not s.strip() for s in samples):
        raise ContractError('shared sample blocks absent')
    types=set(rna[rule['cell_type_key']])&set(atac[rule['cell_type_key']])
    for label in list(types):
        if not label or any(not (((frame[sk]==sample)&(frame['Cell.Type']==label)).any())
            for frame,sk in [(rna,rs),(atac,ats)] for sample in samples):
            types.remove(label)
    types=sorted(types)
    if not types:raise ContractError('no common sample/type composition')
    quota={label:min(rule['max_cells_per_sample']//len(types),min(int(((frame[sk]==sample)&(frame['Cell.Type']==label)).sum())
        for frame,sk in [(rna,rs),(atac,ats)] for sample in samples)) for label in types}
    if min(quota.values())<1:raise ContractError('empty sample/type quota')
    result={}
    for assay,frame,sk in [('rna',rna,rs),('atac',atac,ats)]:
        rows=[]
        for sample in samples:
            for label in types:
                pool=frame.loc[(frame[sk]==sample)&(frame['Cell.Type']==label)]
                chosen=sorted(pool.to_dict('records'),key=lambda row:hashlib.sha256(f'{seed}:{assay}:{sample}:{label}:{row["Barcode"]}'.encode()).hexdigest())[:quota[label]]
                rows.extend(chosen)
        result[assay]=rows
    return result,{'samples':samples,'types':types,'per_type_quota':quota,
        'cells_per_sample':sum(quota.values()),'sample_count':len(samples),
        'pairing':'same biological sample identifiers; separate nuclei, no cell pairing or verified donor crosswalk'}


def read_10x_rows(path,rows):
    with h5py.File(path,'r') as f:
        m=f['matrix'];shape=tuple(int(v) for v in m['shape'][:])
        cells=m['barcodes'].asstr()[:].tolist()
        names=m['features/name'].asstr()[:].tolist();ids=m['features/id'].asstr()[:].tolist()
        if len(set(ids))!=len(ids) or len(cells)!=len(set(cells)) or shape!=(len(names),len(cells)):
            raise ContractError('10x feature or cell identity mismatch')
        if set(m['features/feature_type'].asstr()[:])!={'Gene Expression'}:
            raise ContractError('RNA feature type differs')
        genome=sorted(set(m['features/genome'].asstr()[:]))
        if genome!=['GRCh38.p12.premrna2']:raise ContractError('official brain RNA reference differs')
        if len(rows)!=len(set(rows)) or min(rows)<0 or max(rows)>=len(cells):
            raise ContractError('10x selection invalid')
        data,ind,ptr=[],[],[0]
        for row in rows:
            lo,hi=map(int,m['indptr'][row:row+2]);data.append(m['data'][lo:hi]);ind.append(m['indices'][lo:hi]);ptr.append(ptr[-1]+hi-lo)
        raw=sparse.csr_matrix((np.concatenate(data),np.concatenate(ind),ptr),shape=(len(rows),len(names)))
    counts_normalization(raw,'counts')
    return raw,unique_symbols(names),{'features':len(names),'cells':len(cells),'feature_ids_sha256':digest(ids),'genome':genome,
        'feature_names_sha256':digest(names),'feature_naming':'exact official name order with standard uniqueness replay, no suffix stripping'}


def panel_spec(genes,tfgenes):
    return {'genes':genes,'tf_genes':tfgenes,'panel_sha256':digest({'genes':genes,'tf_genes':tfgenes})}


def subset_genes(values,allgenes,genes):
    lookup={g:i for i,g in enumerate(allgenes)}
    mapping=np.array([lookup.get(g,-1) for g in genes])
    rows=np.flatnonzero(mapping>=0)
    # Preserve declared panel positions and explicit missing features. Zeros in
    # counts do not become observed model graph scores; inference records NaN.
    if len(rows)==len(genes):return values[:,mapping].tocsr()
    coo=values[:,mapping[rows]].tocoo()
    return sparse.csr_matrix((coo.data,(coo.row,rows[coo.col])),shape=(values.shape[0],len(genes)))


def prepare(plan,resources,out,ctx,check=None):
    if (out/'prepared.receipt.json').exists():
        return require_prepared(out,ctx,plan)[0]
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ContractError('uncertified partial preparation requires a new namespace')
    rna=official_metadata(resources,'brain_rna_metadata');atac=official_metadata(resources,'brain_atac_metadata')
    bp=Path(resources['brain_rna']['path'])
    down=load_json(str(bp)+'.download.json')
    if down['sha256']!=resources['brain_rna']['sha256'] or down['bytes']!=bp.stat().st_size:
        raise ContractError('official RNA download identity differs')
    with h5py.File(bp) as f:bc=f['matrix/barcodes'].asstr()[:].tolist()
    with h5py.File(resources['brain_atac']['path']) as f:
        ac=index(f,'obs');peaks=index(f,'var')
        genome=set(column(f,'var','genome'))
        if not genome or any('GRCh38' not in v and 'hg38' not in v for v in genome):raise ContractError('brain ATAC reference absent')
    rlookup=exact_metadata(rna,bc);alookup=exact_metadata(atac,ac)
    selection,composition=balanced_selection(rna,atac,plan['brain'],plan['seed'])
    rrows=[rlookup[r['Barcode']] for r in selection['rna']];arows=[alookup[r['Barcode']] for r in selection['atac']]
    raw,allgenes,source=read_10x_rows(bp,rrows);cp,log=counts_normalization(raw,'counts')
    detection=np.asarray((raw>0).mean(0)).ravel();det=dict(zip(allgenes,detection))
    coordinates=read_coordinates(resources['gene_coordinates']['path'],plan['panel']['promoter_bp'])
    candidates=sorted(set(coordinates)&tokenizable_genes(resources))
    opportunity=locus_links(peaks,candidates,coordinates).getnnz(0)
    detected=[g for g,count in zip(candidates,opportunity) if count>0 and det.get(g,0)>=plan['panel']['detection_floor']]
    tf_symbols={g for m in read_meme(resources['meme']['path']).values() for g in m['symbols']}
    tf=sorted(set(detected)&tf_symbols)
    if len(tf)>plan['panel']['max_genes']:raise ContractError('TF count exceeds frozen cap')
    non=sorted(set(detected)-set(tf),key=lambda g:(-det[g],g))
    newgenes=sorted(tf+non[:plan['panel']['max_genes']-len(tf)])
    if len(newgenes)<plan['panel']['min_genes'] or len(tf)<plan['panel']['min_tfs']:
        raise ContractError('source-corrected panel support below registered minimum')
    historical=load_json(resources['original_panel']['path']);hg=historical['genes'];ht=[g for g in hg if historical['tf_flags'][g]]
    panels={'historical_panel_corrected_source':panel_spec(hg,ht),'brain_selected_panel':panel_spec(newgenes,tf)}
    selection_report={'brain':composition,'rna_cells':selection['rna'],'atac_cells':selection['atac'],
        'brain_rna_metadata_missing_cells':len(bc)-len(rna),'brain_atac_metadata_missing_cells':len(ac)-len(atac),
        'new_panel_selection_scale':'detection over the registered metadata-balanced genuine brain RNA cell pool; no model outcomes',
        'historical_overlap':len(set(newgenes)&set(hg)),'withdrawn_readouts':plan['models']['excluded_historical_readouts'],
        'source':source,'brain_rna_sha256':resources['brain_rna']['sha256'],'brain_atac_sha256':resources['brain_atac']['sha256']}
    write_json(out/'panels.json',panels);write_json(out/'selection.json',selection_report)
    pbmeta=pd.read_csv(resources['pbmc_metadata']['path'],dtype=str,keep_default_na=False)
    if pbmeta.Barcode.duplicated().any():raise ContractError('PBMC metadata has duplicate barcode')
    with h5py.File(resources['pbmc_rna']['path']) as f:pbc=index(f,'obs');pgenes=index(f,'var')
    with h5py.File(resources['pbmc_atac']['path']) as f:pac=index(f,'obs');ppeaks=index(f,'var')
    if pbc!=pac:raise ContractError('PBMC paired source exact row order mismatch')
    pmeta=pbmeta.loc[pbmeta.Barcode.isin(pbc)&(~pbmeta['Cell.Type'].isin(['','NA','nan','NaN']))]
    pchosen=sorted(pmeta.to_dict('records'),key=lambda r:hashlib.sha256(f'{plan["seed"]}:pbmc:{r["Barcode"]}'.encode()).hexdigest())[:plan['pbmc']['max_cells']]
    if not pchosen:raise ContractError('no typed paired PBMC cells')
    plook={c:i for i,c in enumerate(pbc)};prows=[plook[r['Barcode']] for r in pchosen]
    praw=read_rows(resources['pbmc_rna']['path'],'X',prows);pcp,plog=counts_normalization(praw,'counts')
    selection_report['pbmc_cells']=pchosen;write_json(out/'selection.json',selection_report)
    outputs=['panels.json','selection.json'];jobs=[]
    for arm,aps in plan['arms'].items():
        panel=panels[arm];genes=panel['genes'];tfgenes=panel['tf_genes']
        for tissue in aps['tissues']:
            job=arm+'-'+tissue;d=out/job;d.mkdir();jobs.append(job)
            if tissue=='brain':
                counts,norm,lv,universe=raw,cp,log,allgenes
                apath,asels,pks,labels=resources['brain_atac']['path'],arows,peaks,[r['Cell.Type'] for r in selection['atac']]
            else:
                counts,norm,lv,universe=praw,pcp,plog,pgenes
                apath,asels,pks,labels=resources['pbmc_atac']['path'],prows,ppeaks,[r['Cell.Type'] for r in pchosen]
            subraw,subcp,sublog=[subset_genes(v,universe,genes) for v in (counts,norm,lv)]
            links=locus_links(pks,genes,coordinates);keep=np.flatnonzero(links.getnnz(1));links=links[keep].tocsr()
            kept=[pks[i] for i in keep];windows=[]
            for p in kept:
                c,lo,hi=parse_peak(p);mid=(lo+hi)//2;w=plan['proxy']['peak_width_bp']//2;windows.append(f'{c}:{max(0,mid-w)}-{mid+w}')
            hits,gc,length=build_motif_hits(windows,tfgenes,resources['genome_fasta']['path'],resources['meme']['path'],pvalue=plan['proxy']['motif_pvalue'],check=check)
            ax=read_rows(apath,'X',asels,keep)
            if not np.isfinite(ax.data).all() or np.any(ax.data<0) or not np.allclose(ax.data,np.round(ax.data),atol=1e-6,rtol=0):raise ContractError('selected ATAC values are not raw integer counts')
            types=sorted(set(labels));lab=np.array(labels);access=np.mean([np.log1p(np.asarray(ax[lab==t].mean(0)).ravel()) for t in types],axis=0)
            block=(hits.T@links.multiply(access[:,None])).toarray();graph=np.zeros((len(genes),len(genes)));tfrows=[genes.index(g) for g in tfgenes];graph[tfrows]=block;np.fill_diagonal(graph,0)
            ng=np.maximum(1,links.getnnz(0));gcv=np.asarray(links.T@gc).ravel()/ng
            de=np.array([float((counts[:,universe.index(g)]>0).mean()) if g in universe else 0 for g in genes])
            sparse.save_npz(d/'raw.npz',subraw);sparse.save_npz(d/'cp10k.npz',subcp)
            coexpression,variance=valid_coexpression(sublog)
            np.savez_compressed(d/'graphs.npz',proxy=graph,coexpression=coexpression,rna_variance=variance,peak_count=links.getnnz(0),gene_length=np.array([coordinates[g][3] for g in genes]),gc=gcv,detection=de)
            write_json(d/'graph.json',{'arm':arm,'tissue':tissue,'cells':len(asels),'types':types,'rna_cells':counts.shape[0],
                'panel_sha256':panel['panel_sha256'],'peak_ids_sha256':digest(kept),'window_ids_sha256':digest(windows),'linked_peaks':len(kept),
                'motif_parameters':plan['proxy'],'ordered_tf_sha256':digest(tfgenes),'source_feature_position':[universe.index(g) if g in universe else None for g in genes],'missing_rna_features':[g for g in genes if g not in universe],'undefined_coexpression_genes':[g for g,v in zip(genes,variance) if v<=0],
                'scope':'metadata-restricted source-corrected motif-accessibility proxy; no causal reference or same-cell brain pairing'})
            outputs += [job+'/'+n for n in ['raw.npz','cp10k.npz','graphs.npz','graph.json']]
            print('prepared',job,'genes',len(genes),'TFs',len(tfgenes),'peaks',len(kept),flush=True)
    report={'status':'prepared','jobs':jobs,'panels_sha256':sha256(out/'panels.json'),'selection_sha256':sha256(out/'selection.json'),
        'source_context':ctx,'claims':'Separate corrective analysis; historical biological scope withdrawn, historical outputs untouched'}
    write_json(out/'prepared.json',report);outputs.append('prepared.json');finish(out,'prepared.receipt.json',ctx,outputs)
    require_prepared(out,ctx,plan)
    return report


def require_prepared(out,ctx,plan):
    report=load_json(out/'prepared.json');panels=load_json(out/'panels.json')
    jobs=[arm+'-'+t for arm,rule in plan['arms'].items() for t in rule['tissues']]
    expected=['prepared.json','panels.json','selection.json']+[job+'/'+f for job in jobs for f in ['raw.npz','cp10k.npz','graphs.npz','graph.json']]
    receipt(out,'prepared.receipt.json',ctx,expected)
    if report.get('status')!='prepared' or report.get('jobs')!=jobs or set(panels)!=set(plan['arms']) or report['source_context']!=ctx or report['panels_sha256']!=sha256(out/'panels.json') or report['selection_sha256']!=sha256(out/'selection.json'):
        raise ContractError('corrective preparation identity drift')
    for arm,panel in panels.items():
        genes,tf=panel['genes'],panel['tf_genes']
        if not genes or len(genes)!=len(set(genes)) or len(tf)!=len(set(tf)) or not set(tf)<=set(genes) or panel['panel_sha256']!=digest({'genes':genes,'tf_genes':tf}):
            raise ContractError('corrective panel identity drift')
        n=len(genes)
        for tissue in plan['arms'][arm]['tissues']:
            d=out/(arm+'-'+tissue);a=np.load(d/'graphs.npz',allow_pickle=False)
            if a['proxy'].shape!=(n,n) or not np.isfinite(a['proxy']).all() or a['coexpression'].shape!=(n,n) or np.isinf(a['coexpression']).any():
                raise ContractError('corrective prepared graph shape or finiteness differs')
            if sparse.load_npz(d/'raw.npz').shape!=sparse.load_npz(d/'cp10k.npz').shape or sparse.load_npz(d/'raw.npz').shape[1]!=n:
                raise ContractError('corrective count dimensions differ from panel')
    return report,panels


def pilot(plan,resources,out,ctx,check=None):
    prepared,panels=require_prepared(out,ctx,plan);d=out/'pilot';d.mkdir(exist_ok=True)
    if (d/'pilot.receipt.json').exists():require_pilot(plan,out,ctx);return load_json(d/'pilot.json')
    outputs=[];rows=[]
    for arm in plan['arms']:
        job=arm+'-brain';panel=panels[arm]
        cp=sparse.load_npz(out/job/'cp10k.npz')[:plan['models']['batch']];raw=sparse.load_npz(out/job/'raw.npz')[:plan['models']['batch']]
        for readout in plan['models']['readouts']:
            td=d/arm/readout;td.mkdir(parents=True,exist_ok=True);tc={**ctx,'readout':readout,'job':job,'panel_sha256':panel['panel_sha256']}
            if (td/'result.receipt.json').exists():receipt(td,'result.receipt.json',tc,['cuda.npz','cpu.npz','pilot.json']);r=load_json(td/'pilot.json')
            else:
                kwargs={'batch':plan['models']['batch'],'seed':plan['seed'],'gpu_policy':plan['gpu'],'check':check}
                g,t=run_readout(readout,panel['genes'],cp,raw,resources,device='cuda',**kwargs)
                c,ct=run_readout(readout,panel['genes'],cp,raw,resources,device='cpu',**kwargs)
                finite=np.isfinite(g)&~np.eye(len(g),dtype=bool);pc=plan['pilot']
                if check:check()
                if not np.array_equal(np.isfinite(g),np.isfinite(c)) or finite.sum()<pc['min_comparable_pairs'] or not np.allclose(g[finite],c[finite],atol=pc['cpu_cuda_atol'],rtol=pc['cpu_cuda_rtol']):
                    raise ContractError('source corrective CPU/CUDA pilot mismatch')
                r={'readout':readout,'job':job,'status':'pilot_success','cuda':t,'cpu':ct,'pairs':int(finite.sum()),'max_difference':float(np.abs(g[finite]-c[finite]).max()),'panel_sha256':panel['panel_sha256'],'scope':'engineering pilot, not biological result'}
                np.savez_compressed(td/'cuda.npz',graph=g);np.savez_compressed(td/'cpu.npz',graph=c);write_json(td/'pilot.json',r);finish(td,'result.receipt.json',tc,['cuda.npz','cpu.npz','pilot.json'])
            rows.append(r);outputs += [arm+'/'+readout+'/'+f for f in ['cuda.npz','cpu.npz','pilot.json','result.receipt.json']]
    report={'status':'pilot_complete','readouts':rows,'source_context':ctx};write_json(d/'pilot.json',report);outputs.append('pilot.json');finish(d,'pilot.receipt.json',ctx,outputs);require_pilot(plan,out,ctx,resources);return report


def validate_telemetry(row, readout, cells, device, *, resources=None, batch=2):
    counts=np.asarray(row.get('coverage_per_gene',[]))
    if row.get('readout')!=readout or row.get('cells')!=cells or row.get('device')!=device or row.get('context')!='frozen_panel' or row.get('batches')!=int(np.ceil(cells/batch)):
        raise ContractError('corrective telemetry execution identity differs')
    maximum = cells
    if readout == 'uce_encoder':
        adapter = HERE/'screg_validation/vendor/fm_readout_uce.py'
        if sha256(adapter) != UCE_ADAPTER_SHA256:
            raise ContractError('UCE token-count semantics require the frozen adapter')
        maximum = cells * UCE_TOKENS_PER_CELL
    if counts.ndim!=1 or not len(counts) or not np.isfinite(counts).all() or np.any(counts<0) or np.any(counts>maximum) or not np.equal(counts,np.round(counts)).all():
        raise ContractError('corrective telemetry coverage counts invalid')
    if resources is not None and row.get('resources_sha256')!=digest({n:resources[n]['sha256'] for n in READOUT_RESOURCES[readout]}):
        raise ContractError('corrective telemetry checkpoint identity differs')
    if device=='cuda' and any(not isinstance(row.get(k),dict) or row[k].get('ready') is not True or row[k].get('blockers') for k in ['gpu_before','gpu_after']):
        raise ContractError('corrective GPU readiness evidence blocked')


def require_pilot(plan,out,ctx,resources=None):
    _,panels=require_prepared(out,ctx,plan)
    d=out/'pilot';expected=['pilot.json']+[arm+'/'+r+'/'+f for arm in plan['arms'] for r in plan['models']['readouts'] for f in ['cuda.npz','cpu.npz','pilot.json','result.receipt.json']]
    receipt(d,'pilot.receipt.json',ctx,expected)
    summary=load_json(d/'pilot.json')
    expected_pairs={(arm+'-brain',r) for arm in plan['arms'] for r in plan['models']['readouts']}
    summary_rows=summary.get('readouts',[])
    if summary.get('status')!='pilot_complete' or summary.get('source_context')!=ctx or len(summary_rows)!=len(expected_pairs) or {(row.get('job'),row.get('readout')) for row in summary_rows}!=expected_pairs:raise ContractError('global corrective pilot identity incomplete')
    for arm in plan['arms']:
        for r in plan['models']['readouts']:
            td=d/arm/r;tc={**ctx,'readout':r,'job':arm+'-brain','panel_sha256':panels[arm]['panel_sha256']};receipt(td,'result.receipt.json',tc,['cuda.npz','cpu.npz','pilot.json'])
            g=np.load(td/'cuda.npz',allow_pickle=False)['graph'];c=np.load(td/'cpu.npz',allow_pickle=False)['graph'];finite=np.isfinite(g)&~np.eye(len(g),dtype=bool)
            row=load_json(td/'pilot.json')
            if next(v for v in summary_rows if v['job']==arm+'-brain' and v['readout']==r)!=row:raise ContractError('pilot summary contradicts saved child')
            cells=min(plan['models']['batch'],sparse.load_npz(out/(arm+'-brain')/'raw.npz').shape[0])
            if row.get('readout')!=r or row.get('job')!=arm+'-brain':raise ContractError('pilot job/readout identity differs')
            for device in ['cpu','cuda']:validate_telemetry(row[device],r,cells,device,resources=resources,batch=plan['models']['batch'])
            if row['cpu'].get('resources_sha256')!=row['cuda'].get('resources_sha256'):raise ContractError('CPU/CUDA checkpoint differs')
            if row.get('status')!='pilot_success' or row.get('panel_sha256')!=panels[arm]['panel_sha256'] or g.shape!=(len(panels[arm]['genes']),)*2 or c.shape!=g.shape or np.isinf(g).any() or np.isinf(c).any() or not np.array_equal(np.isfinite(g),np.isfinite(c)) or finite.sum()<plan['pilot']['min_comparable_pairs'] or not np.allclose(g[finite],c[finite],atol=plan['pilot']['cpu_cuda_atol'],rtol=plan['pilot']['cpu_cuda_rtol']):raise ContractError('corrective pilot graph verification failed')
            error=float(np.abs(g[finite]-c[finite]).max())
            if row.get('pairs')!=int(finite.sum()) or row.get('max_difference')!=error:raise ContractError('corrective pilot measured comparison drift')


def infer(plan,resources,out,ctx,check=None):
    prepared,panels=require_prepared(out,ctx,plan);require_pilot(plan,out,ctx,resources);outputs=[];rows=[];started=time.monotonic()
    cumulative_check=check
    def check():
        if cumulative_check:cumulative_check()
        if time.monotonic()-started>plan['models']['max_seconds']:raise ContractError('corrective inference budget exhausted')
        require_gpu(plan['gpu'])
    for job in prepared['jobs']:
        panel=panels[job.rsplit('-',1)[0]];cp=sparse.load_npz(out/job/'cp10k.npz');raw=sparse.load_npz(out/job/'raw.npz')
        for readout in plan['models']['readouts']:
            d=out/job/readout;d.mkdir(exist_ok=True);tc={**ctx,'job':job,'readout':readout,'panel':panel['panel_sha256']}
            if (d/'result.receipt.json').exists():receipt(d,'result.receipt.json',tc,['graph.npz','coverage.json']);r=load_json(d/'coverage.json')
            else:
                check();g,t=run_readout(readout,panel['genes'],cp,raw,resources,device='cuda',batch=plan['models']['batch'],seed=plan['seed'],gpu_policy=plan['gpu'],check=check);check()
                status=coverage_status(g,t,panel,plan['models']);r={**t,**status,'job':job,'panel_sha256':panel['panel_sha256']}
                np.savez_compressed(d/'graph.npz',graph=g);write_json(d/'coverage.json',r);finish(d,'result.receipt.json',tc,['graph.npz','coverage.json'])
            saved=np.load(d/'graph.npz',allow_pickle=False)['graph'];actual=coverage_status(saved,r,panel,plan['models'])
            if any(r.get(k)!=v for k,v in actual.items()):raise ContractError('cached corrective coverage drift')
            rows.append(r);outputs += [job+'/'+readout+'/'+f for f in ['graph.npz','coverage.json','result.receipt.json']]
            print('inferred',job,readout,'eligible',r['primary_statistics_eligible'],flush=True)
    if len(rows)!=plan['models']['registered_tasks']:raise ContractError('inference task accounting differs from registration')
    report={'status':'inference_complete','tasks':len(rows),'eligible':sum(r['primary_statistics_eligible'] for r in rows),'rows':rows,'source_context':ctx}
    write_json(out/'inference.json',report);outputs.append('inference.json');finish(out,'inference.receipt.json',ctx,outputs);require_inference(plan,out,ctx,resources);return report


def require_inference(plan,out,ctx,resources=None):
    prepared,panels=require_prepared(out,ctx,plan);require_pilot(plan,out,ctx)
    expected=['inference.json']+[job+'/'+readout+'/'+f for job in prepared['jobs'] for readout in plan['models']['readouts'] for f in ['graph.npz','coverage.json','result.receipt.json']]
    receipt(out,'inference.receipt.json',ctx,expected)
    summary=load_json(out/'inference.json')
    if summary.get('status')!='inference_complete' or summary.get('tasks')!=plan['models']['registered_tasks'] or summary.get('source_context')!=ctx or len(summary.get('rows',[]))!=plan['models']['registered_tasks']:
        raise ContractError('inference summary task identity incomplete')
    pairs={(r.get('job'),r.get('readout')) for r in summary['rows']}
    if pairs!={(j,r) for j in prepared['jobs'] for r in plan['models']['readouts']}:raise ContractError('inference summary missing or duplicate identities')
    for job in prepared['jobs']:
        p=panels[job.rsplit('-',1)[0]]
        for r in plan['models']['readouts']:
            d=out/job/r;tc={**ctx,'job':job,'readout':r,'panel':p['panel_sha256']};receipt(d,'result.receipt.json',tc,['graph.npz','coverage.json'])
            c=load_json(d/'coverage.json')
            if c.get('job')!=job or c.get('panel_sha256')!=p['panel_sha256']:raise ContractError('inference child job/panel differs')
            validate_telemetry(c,r,sparse.load_npz(out/job/'raw.npz').shape[0],'cuda',resources=resources,batch=plan['models']['batch'])
            pilotrow=load_json(out/'pilot'/job.rsplit('-',1)[0]/r/'pilot.json')
            if c.get('resources_sha256')!=pilotrow['cuda'].get('resources_sha256'):raise ContractError('inference checkpoint differs from verified panel pilot')
            g=np.load(d/'graph.npz',allow_pickle=False)['graph'];s=coverage_status(g,c,p,plan['models'])
            if any(c.get(k)!=v for k,v in s.items()):raise ContractError('source corrective coverage inconsistent')
            if next(v for v in summary['rows'] if v['job']==job and v['readout']==r)!=c:raise ContractError('inference summary contradicts child coverage')
    if summary.get('eligible')!=sum(r['primary_statistics_eligible'] for r in summary['rows']):raise ContractError('inference eligible tally differs')
    return prepared,panels


def analyze(plan,resources,out,ctx,check=None):
    prepared,panels=require_inference(plan,out,ctx,resources);slots=[];started=time.monotonic()
    cumulative_check=check
    def check():
        if cumulative_check:cumulative_check()
        if time.monotonic()-started>plan['statistics']['max_seconds']:raise ContractError('corrective statistics budget exhausted')
    for job in prepared['jobs']:
        panel=panels[job.rsplit('-',1)[0]];genes=panel['genes'];tf=np.array([genes.index(g) for g in panel['tf_genes']]);n=len(genes)
        ii,jj=np.repeat(tf,n),np.tile(np.arange(n),len(tf));ok=ii!=jj;ii,jj=ii[ok],jj[ok]
        a=dict(np.load(out/job/'graphs.npz',allow_pickle=False));g=a['proxy'];co=a['coexpression']
        structural=np.column_stack([a[k][jj] for k in ['peak_count','gene_length','detection','gc']]);controls=np.column_stack([structural,rankdata(co[ii,jj])])
        baseline_finite=bool(np.isfinite(co[ii,jj]).all())
        names=[];scores=[];coverage={}
        for r in plan['models']['readouts']:
            c=load_json(out/job/r/'coverage.json');coverage[r]=c
            if c['primary_statistics_eligible'] and baseline_finite:names.append(r);scores.append(np.load(out/job/r/'graph.npz',allow_pickle=False)['graph'])
        for degree in plan['statistics']['degree_modes']:
            basec=structural if degree=='none' else np.column_stack([structural,(g>0).sum(1)[ii],(g>0).sum(0)[jj]])
            br=partial_spearman(co[ii,jj],g[ii,jj],basec) if baseline_finite else None
            for null in plan['statistics']['null_modes']:
                d=out/job/('statistics-'+degree+'-'+null);d.mkdir(exist_ok=True);tc={**ctx,'job':job,'degree':degree,'null':null,'readouts':names}
                if names:
                    if (d/'result.receipt.json').exists():receipt(d,'result.receipt.json',tc,['rows.json']);rows=load_json(d/'rows.json')
                    else:
                        strata=np.zeros(n,int)
                        if null=='restricted':
                            ob=np.digitize(a['peak_count'],np.unique(np.quantile(a['peak_count'],[.2,.4,.6,.8])));gb=np.digitize(a['gc'],np.unique(np.quantile(a['gc'],[1/3,2/3])));strata=15*np.isin(np.arange(n),tf)+3*ob+gb
                        rows=batched_components(np.stack(scores),g,tf,controls,strata=strata,permutations=plan['statistics']['permutations'],seed=plan['seed'],degree_mode=degree,backend=plan['statistics']['backend'],gpu_policy=plan['gpu'],check=check,checkpoint=d/'draw-progress.json',checkpoint_context=tc)
                        for r,s,row in zip(names,scores,rows):
                            row.pop('nulls');fm=partial_spearman(s[ii,jj],g[ii,jj],basec)
                            row.update(job=job,readout=r,null_mode=null,comparable_fm=fm,comparable_baseline=br,comparable_difference=fm-br,concordance=float(spearmanr(s[ii,jj],co[ii,jj]).statistic),measured=True,correction_input_p=row['pIUT'])
                        write_json(d/'rows.json',rows);finish(d,'result.receipt.json',tc,['rows.json'])
                    if len(rows)!=len(names) or [r.get('readout') for r in rows]!=names:raise ContractError('statistics readout identities incomplete or duplicated')
                    for row in rows:
                        if row.get('job')!=job or row.get('degree_mode')!=degree or row.get('null_mode')!=null or row.get('permutations')!=plan['statistics']['permutations'] or row.get('backend')!=plan['statistics']['backend'] or row.get('measured') is not True:raise ContractError('statistics cached slot identity or draw registration differs')
                        if any(type(row.get(k)) not in (float,int) or not np.isfinite(row[k]) for k in ['observed','comparable_fm','comparable_baseline','comparable_difference','concordance']):raise ContractError('statistics cached effect or contrast nonfinite')
                        if not np.isclose(row['comparable_difference'],row['comparable_fm']-row['comparable_baseline'],rtol=1e-10,atol=1e-12):raise ContractError('statistics contrast arithmetic drift')
                        probs=[row.get(k) for k in ['pM','pD','pIUT','correction_input_p']]
                        if any(type(v) not in (float,int) or not np.isfinite(v) or not 0<v<=1 for v in probs) or row['pIUT']!=max(row['pM'],row['pD']) or row['correction_input_p']!=row['pIUT']:raise ContractError('statistics IUT or correction input inconsistent')
                        if any(abs(v*(plan['statistics']['permutations']+1)-round(v*(plan['statistics']['permutations']+1)))>1e-8 for v in probs[:2]):raise ContractError('statistics p is not plus-one Monte Carlo')
                    slots.extend(rows)
                for r in plan['models']['readouts']:
                    if r not in names:slots.append({'job':job,'readout':r,'degree_mode':degree,'null_mode':null,'measured':False,'observed':None,'pM':None,'pD':None,'pIUT':None,'permutations':0,'correction_input_p':1.,'reason':'incomplete_model_or_baseline_frozen_pair_coverage; no zero imputation','coverage':{k:coverage[r][k] for k in ['gene_coverage_fraction','evaluated_pair_coverage_fraction']}})
                print('analyzed',job,degree,null,'measured',len(names),flush=True)
    expected=len(prepared['jobs'])*len(plan['models']['readouts'])*len(plan['statistics']['degree_modes'])*len(plan['statistics']['null_modes'])
    expected_slots={(j,r,d,nm) for j in prepared['jobs'] for r in plan['models']['readouts'] for d in plan['statistics']['degree_modes'] for nm in plan['statistics']['null_modes']}
    actual_slots=[tuple(r[k] for k in ['job','readout','degree_mode','null_mode']) for r in slots]
    if len(slots)!=expected or len(actual_slots)!=len(set(actual_slots)) or set(actual_slots)!=expected_slots:raise ContractError('corrective registered family accounting differs')
    qb,ph=adjust([r['correction_input_p'] for r in slots]),adjust([r['correction_input_p'] for r in slots],'Holm')
    for row,q,p in zip(slots,qb,ph):row.update(q_IUT_BY=float(q),p_IUT_Holm=float(p))
    report={'status':'complete','source_context':ctx,'rows':slots,'registered_hypotheses':expected,'measured_hypotheses':sum(r['measured'] for r in slots),'limits':plan['limits'],'historical_results_replaced':False,'donor_intervals':None,'scope':'Post-review corrective fixed-panel conditional analysis; no independent donor CI, causal ground truth, gate validity or pretraining exclusion'}
    write_json(out/'results.json',report);finish(out,'analysis.receipt.json',ctx,['results.json']);return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--task',choices=['prepare','pilot','infer','analyze'],required=True);p.add_argument('--plan',type=Path,default=HERE/'corrective_plan.json');p.add_argument('--resources',type=Path,default=HERE/'corrective_resources.local.json');p.add_argument('--authorization',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
    plan=load_json(a.plan);resources=load_json(a.resources);authorize(plan,a.authorization,'corrective',a.task)
    if not a.out.resolve().is_relative_to(HERE):raise ContractError('source corrective outputs must remain in validation namespace')
    if not a.execute:print({'status':'planned_only','task':a.task});return
    ctx=identity(plan,resources)
    maximum=plan['pilot']['max_seconds'] if a.task=='pilot' else plan['models']['max_seconds'] if a.task=='infer' else plan['statistics']['max_seconds']
    with budget_stage(a.out,a.task,maximum) as check:
        result=globals()[a.task](plan,resources,a.out,ctx,check=check)
    print({'task':a.task,'status':result['status'],'out':str(a.out)})


if __name__=='__main__':main()
