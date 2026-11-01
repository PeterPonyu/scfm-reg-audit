#!/usr/bin/env python3
"""Generate the new main-text tables and in-text number macros from evidence/R3_RESULTS.json.

Every new number in the manuscript comes from this script, so text, tables and evidence agree.
"""
import json
import math
from pathlib import Path

W = Path(__file__).resolve().parents[1]
R = json.loads((W / 'evidence/R3_RESULTS.json').read_text())
OUT = W / 'source/parts'
JOBS = ['brain_selected_panel-brain', 'brain_selected_panel-pbmc', 'historical_panel_corrected_source-brain']
CTX = {'brain_selected_panel-brain': 'New / brain', 'brain_selected_panel-pbmc': 'New / PBMC',
       'historical_panel_corrected_source-brain': 'Old / brain'}
ORDER = ['geneformer_embed', 'geneformer_attention', 'random_init_geneformer', 'scgpt_encoder', 'scfoundation_encoder', 'uce_encoder']
SHORT = {'geneformer_embed': 'GF emb.', 'geneformer_attention': 'GF attn.', 'random_init_geneformer': 'Random GF',
         'scgpt_encoder': 'scGPT', 'scfoundation_encoder': 'scFoundation', 'uce_encoder': 'UCE'}
macros = {}


def num(v, nd=4):
    s = f'{v:.{nd}f}'
    return ('$-$' + s[1:]) if s.startswith('-') else s


def num_up(v, nd=4):
    """Upper bound at nd decimals, so that a stated 'at most' or 'below' holds for the raw value."""
    return num(math.ceil(v * 10 ** nd) / 10 ** nd, nd)


def mnum(v, nd=4):
    return f'{v:.{nd}f}'


def pct(v, nd=1):
    return f'{100 * v:.{nd}f}'


def macro(name, value):
    macros[name] = value


def write(name, text):
    (OUT / name).write_text(text)


# ---------------- Table 1: coverage ----------------
cs = R['common_support']
import sys
sys.path.insert(0, str(W / 'analysis'))
import r3lib as L
sup, gl = {}, {}
for _j in JOBS:
    _ctx = L.load_context(_j)
    sup[_j] = L.common_support(_ctx)[4]
    gl[_j] = L.common_support(_ctx, pair_readouts=L.GENE_LEVEL)[4]
if len(sup) == 3:
    reg = {'brain_selected_panel-brain': 16, 'brain_selected_panel-pbmc': 0, 'historical_panel_corrected_source-brain': 0}
    rows = []
    for j in JOBS:
        s, sgl = sup[j], gl[j]
        rows.append(f"{CTX[j]} & {reg[j]}/24 & {s['genes_S']:,}/{s['genes_panel']:,} & {s['tfs_S']}/{s['tfs_panel']} & "
                    f"{s['pairs_E']:,} ({pct(s['fraction_registered_pairs_retained'])}\\%) & "
                    f"{sgl['pairs_E']:,} ({pct(sgl['fraction_registered_pairs_retained'])}\\%) & 24/24 \\\\")
    write('t_coverage.tex', r"""\begin{table}[!htbp]
\centering\small
\caption{Registered admission and common support in the three source-matched contexts. The registered rule required every frozen non-self pair to be finite for every readout and for the comparator. The common support keeps panel genes with positive RNA variance and finite gene-level scores for every readout; its six-readout edges are also finite for the pair-level attention readout, so all readouts are compared on identical edges. Percentages are of the registered non-self TF--target pairs.}\label{tab:coverage}
\setlength{\tabcolsep}{3.2pt}
\begin{tabular}{@{}lrrrrrr@{}}\toprule
 & Registered & \multicolumn{2}{c}{Common support} & Six-readout & Gene-level & Common-support \\\cmidrule(lr){3-4}
Context & measured & Genes & TFs & edges & edges & measured \\\midrule
""" + '\n'.join(rows) + r"""
\bottomrule\end{tabular}
\end{table}
""")

# ---------------- Table 2: calibration ----------------
cal = R['calibration']
CAL_CTX = [('brain', 'Mixed-source (historical)'), ('pbmc', 'Old-panel PBMC (historical)'), ('genuine_brain', 'New / brain'),
           ('new_pbmc', 'New / PBMC'), ('old_brain', 'Old / brain'), ('bmmc_13272', 'Bone marrow 13272'), ('bmmc_28483', 'Bone marrow 28483')]
if all(c in cal for c, _ in CAL_CTX):
    def cell(c, null, sk, comp):
        r = next(x for x in cal[c][null] if x['scenario'] == sk)[comp]
        flag = '' if r['nominal_bradley_and_wilson'] else ('$^{\\ast}$' if r['rate'] < 0.025 else '$^{\\dagger}$')
        return f"{r['count']}{flag}"
    lines = []
    for c, name in CAL_CTX:
        cells = [name]
        for sk in ('iid', 'shock_stratifier_scaled', 'shock_detection_scaled'):
            for null in ('unrestricted', 'restricted'):
                cells.append(cell(c, null, sk, 'joint'))
        for sk in ('shock_stratifier_scaled', 'shock_detection_scaled'):
            cells.append(cell(c, 'restricted', sk, 'M'))
        for sk in ('shock_stratifier_scaled', 'shock_detection_scaled'):
            cells.append(cell(c, 'restricted', sk, 'D'))
        lines.append(' & '.join(cells) + r' \\')
    write('t_calibration.tex', r"""\begin{table}[!htbp]
\centering\small
\caption{Synthetic calibration of the component tests: rejections out of 200 synthetic readouts at raw $p\leq0.05$ with 999 transformations per component, in two historical contexts and five analysis contexts. The two historical contexts and genuine brain with the new panel were specified before the analyses; the new panel in PBMC, the old panel in brain and bone-marrow donors 13272 and 28483 were added after the family results. Joint: intersection--union test with unrestricted (U) or restricted (R) transformations under independent noise (Ind.), shocks scaled by the stratifying variables (Strata) and shocks scaled by RNA detection (Det.). The last four columns give the restricted gene-label (M) and target-shuffle (D) components under the two shock scenarios. A dagger marks rates above the nominal rule (Wilson 95\% interval containing 0.05 and rate between 0.025 and 0.075) and an asterisk conservative rates below 0.025. The historical contexts use the published synthetic seeds, and their unrestricted counts reproduce the published values exactly.}\label{tab:calibration}
\setlength{\tabcolsep}{2.6pt}
\begin{tabular}{@{}lrrrrrrrrrr@{}}\toprule
 & \multicolumn{6}{c}{Joint test} & \multicolumn{2}{c}{M, restricted} & \multicolumn{2}{c}{D, restricted} \\\cmidrule(lr){2-7}\cmidrule(lr){8-9}\cmidrule(lr){10-11}
 & \multicolumn{2}{c}{Ind.} & \multicolumn{2}{c}{Strata} & \multicolumn{2}{c}{Det.} & & & & \\\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}
Context & U & R & U & R & U & R & Strata & Det. & Strata & Det. \\\midrule
""" + '\n'.join(lines) + r"""
\bottomrule\end{tabular}
\end{table}
""")
    rs = [next(x for x in cal[c]['restricted'] if x['scenario'] == sk)['joint'] for c, _ in CAL_CTX for sk in ('shock_stratifier_scaled', 'shock_detection_scaled')]
    us = [next(x for x in cal[c]['unrestricted'] if x['scenario'] == sk)['joint'] for c, _ in CAL_CTX for sk in ('shock_stratifier_scaled', 'shock_detection_scaled')]
    macro('CalRok', str(sum(1 for r in rs if r['rate'] <= 0.075 and r['wilson95'][0] <= 0.05)))
    macro('CalUbad', str(sum(1 for r in us if r['rate'] > 0.075)))

# ---------------- Table 3: common-support results ----------------
six = cs['six']
if six.get('family'):
    rows = six['rows']

    def g(j, m, d, n):
        return next(r for r in rows if (r['job'], r['readout'], r['degree_mode'], r['null_mode']) == (j, m, d, n))
    lines = []
    for j in JOBS:
        lines.append(r'\multicolumn{7}{@{}l}{\textit{' + CTX[j] + r'}} \\')
        for m in ORDER:
            a, b = g(j, m, 'none', 'restricted'), g(j, m, 'own', 'restricted')
            au, bu = g(j, m, 'none', 'unrestricted'), g(j, m, 'own', 'unrestricted')
            minq = min(x['q_BY'] for x in (a, b, au, bu))
            minh = min(x['p_Holm'] for x in (a, b, au, bu))
            lines.append(f"\\quad {SHORT[m]} & {num(a['observed'], 5)} & {num(a['pIUT'])} & {num(b['observed'], 5)} & "
                         f"{num(b['pIUT'])} & {num(minq)} & {num(minh)} \\\\")
    f = six['family']
    fr = six['family_restricted_only']
    write('t_common.tex', r"""\begin{table}[!htbp]
\centering\small
\caption{Common-support results for all 72 hypotheses (9,999 transformations per component). $\rho$ is the incremental partial Spearman correlation on identical edges; joint $p$ is the restricted intersection--union probability. The last two columns give the smallest BY and Holm value over the four specifications of each row (no degree or own degree, restricted or unrestricted) after adjustment across all 72 hypotheses. Minimum raw joint $p$ """ + num(f['min_raw_pIUT']) + r"""; minimum BY """ + num(f['min_BY']) + r"""; minimum Holm """ + num(f['min_Holm']) + r"""; selections: BY """ + str(f['BY_selected']) + r""", Holm """ + str(f['Holm_selected']) + r""". Probabilities and adjusted values of 1.0000 are computed, not administrative inputs. GF, Geneformer; emb., embedding; attn., attention; BY, Benjamini--Yekutieli.}\label{tab:common}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}lrrrrrr@{}}\toprule
 & \multicolumn{2}{c}{No degree, restricted} & \multicolumn{2}{c}{Own degree, restricted} & \multicolumn{2}{c}{Family of 72} \\\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}
Readout & $\rho$ & Joint $p$ & $\rho$ & Joint $p$ & Min.\ BY & Min.\ Holm \\\midrule
""" + '\n'.join(lines) + r"""
\bottomrule\end{tabular}
\end{table}
""")
    macro('CSminp', num(f['min_raw_pIUT']))
    macro('CSminBY', num(f['min_BY']))
    macro('CSminHolm', num(f['min_Holm']))
    macro('CSnBY', str(f['BY_selected']))
    macro('CSnHolm', str(f['Holm_selected']))
    macro('CSRminp', num(fr['min_raw_pIUT']))
    obs = [r['observed'] for r in rows]
    macro('CSrhoMin', num(min(obs), 5))
    macro('CSrhoMax', num(max(obs), 5))
    glf = cs['gene_level'].get('family')
    if glf:
        macro('GLminp', num(glf['min_raw_pIUT']))
        macro('GLminHolm', num(glf['min_Holm']))
        macro('GLminBY', num(glf['min_BY']))
        macro('GLnBY', str(glf['BY_selected']))
        macro('GLnHolm', str(glf['Holm_selected']))

    # registered vs common support agreement (brain/new panel, four registered readouts)
    reg_path = W / 'analysis/resolution/results.json'
    if not reg_path.exists():  # layout of the public deposit
        reg_path = W.parent / 'resolution_9999/results/results.json'
    reg = json.loads(reg_path.read_text())
    regm = [r for r in reg['rows'] if r['measured']]
    diffs = [abs(r['observed'] - g('brain_selected_panel-brain', r['readout'], r['degree_mode'], r['null_mode'])['observed']) for r in regm]
    macro('CSregMaxDiff', num(max(diffs), 5))

    def seltext(fam):
        if fam['BY_selected'] == 0 and fam['Holm_selected'] == 0:
            return 'no hypothesis is selected by BY or Holm'
        items = []
        for r in fam['selected']:
            items.append(f"{SHORT.get(r['readout'], r['readout'])} ({CTX.get(r.get('job'), r.get('context', ''))}, "
                         f"{r.get('degree_mode', r.get('design', ''))} degree, {r.get('null_mode', 'unrestricted')}; "
                         f"$\\rho={mnum(r['observed'], 4)}$, Holm {num(r['p_Holm'])}, BY {num(r['q_BY'])})")
        return f"BY selects {fam['BY_selected']} and Holm {fam['Holm_selected']} hypotheses: " + '; '.join(items)
    macro('CSselText', seltext(f))
    macro('CSdOnly', str(sum(1 for r in rows if r['pD'] < 0.01 and r['pM'] > 0.1)))

    if glf:
        macro('GLselText', seltext(glf))

# ---------------- symmetric bootstrap macros ----------------
sb = R.get('symmetric_bootstrap')
if sb:
    co_none = [b['coexpression_rho'] for b in sb if b['degree_mode'] == 'none']
    macro('SymCoNoneMin', num(min(co_none), 3))
    macro('SymCoNoneMax', num(max(co_none), 3))
    neg = [(b['job'], r['readout']) for b in sb if b['degree_mode'] == 'none' for r in b['readouts'] if r['difference_ci95'][1] < 0]
    pos = [(b['job'], r['readout']) for b in sb if b['degree_mode'] == 'none' for r in b['readouts'] if r['difference_ci95'][0] > 0]
    macro('SymNegNone', str(len(neg)))
    macro('SymPosNone', str(len(pos)))
    nego = [(b['job'], r['readout']) for b in sb if b['degree_mode'] == 'own' for r in b['readouts'] if r['difference_ci95'][1] < 0]
    poso = [(b['job'], r['readout']) for b in sb if b['degree_mode'] == 'own' for r in b['readouts'] if r['difference_ci95'][0] > 0]
    macro('SymNegOwn', str(len(nego)))
    macro('SymPosOwn', str(len(poso)))

symc = cs['symmetric'].get('coexpression_rows', [])
if symc:
    pr = [r['pIUT'] for r in symc if r['degree_mode'] == 'none' and r['null_mode'] == 'restricted']
    if len(pr) == 3:
        macro('SymCoPmax', num(max(pr)))

# ---------------- grid ----------------
gr = R['grid']
if gr.get('pooled_180'):
    DES = ['none', 'own', 'outdeg_only', 'indeg_only', 'log_degree', 'rank_degree', 'spline_degree', 'strength',
           'external_fibro', 'external_other']
    DL = {'none': 'No degree', 'own': 'Own linear', 'outdeg_only': 'Out-degree only', 'indeg_only': 'In-degree only',
          'log_degree': 'Log degree', 'rank_degree': 'Rank degree', 'spline_degree': 'Spline degree', 'strength': 'Weighted degree',
          'external_fibro': 'Fibroblast-mix degree', 'external_other': 'Other-tissue degree'}
    HIST = {'none': 0, 'own': 7, 'outdeg_only': 1, 'indeg_only': 0, 'log_degree': 1, 'rank_degree': 1, 'spline_degree': 6,
            'strength': 0, 'external_fibro': 0, 'external_other': 0}
    lines = []
    for d in DES:
        per = [gr['per_context_design'][f'{j}|{d}'] for j in JOBS]
        lines.append(f"{DL[d]} & {HIST[d]}/13 & " + ' & '.join(f"{p['BY_selected']}/{p['Holm_selected']}" for p in per) +
                     f" & {num(min(p['min_Holm'] for p in per))} \\\\")
    p180 = gr['pooled_180']
    write('t_grid.tex', r"""\begin{table}[!htbp]
\centering\small
\caption{Selections under ten degree designs. The historical column counts IUT-BY selections of thirteen mixed-source and PBMC rows within each context and design. The source-matched columns give IUT-BY/IUT-Holm selections of six readouts within each context and design on common-support edges with unrestricted transformations and 9,999 draws; the last column is the smallest Holm value in that design across the three contexts. Pooling all 180 source-matched tests gives minimum BY """ + num(p180['min_BY']) + r""" and minimum Holm """ + num(p180['min_Holm']) + r""" (selections: BY """ + str(p180['BY_selected']) + r""", Holm """ + str(p180['Holm_selected']) + r"""). BY, Benjamini--Yekutieli; IUT, intersection--union test.}\label{tab:grid}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}lrrrrr@{}}\toprule
Degree design & Historical & New / brain & New / PBMC & Old / brain & Min.\ Holm \\\midrule
""" + '\n'.join(lines) + r"""
\bottomrule\end{tabular}
\end{table}
""")
    def gget(j, m, d):
        return next(r for r in gr['rows'] if (r['job'], r['readout'], r['design']) == (j, m, d))
    ratios = [gget(j, m, 'none')['nullM_sd'] / gget(j, m, 'own')['nullM_sd'] for j in JOBS for m in ORDER if m != 'random_init_geneformer']
    rratios = [gget(j, 'random_init_geneformer', 'none')['nullM_sd'] / gget(j, 'random_init_geneformer', 'own')['nullM_sd'] for j in JOBS]
    macro('GRwidthMin', f'{min(ratios):.1f}')
    macro('GRwidthMax', f'{max(ratios):.1f}')
    macro('GRwidthRandMax', f'{max(rratios):.1f}')
    sbd = gr['selection_by_design']
    macro('GRownBY', str(sbd['own']['BY']))
    macro('GRsplineBY', str(sbd['spline_degree']['BY']))
    macro('GRanySel', str(sum(v['BY'] + v['Holm'] for v in sbd.values())))
    persel = [(k.split('|')[0], k.split('|')[1], v) for k, v in gr['per_context_design'].items() if v['BY_selected'] or v['Holm_selected']]
    if not persel:
        macro('GRselText', 'no design selects any row within its context-by-design family')
    else:
        macro('GRselText', 'selections occur within ' + str(len(persel)) + ' context-by-design families: ' +
              '; '.join(f"{CTX[j]}, {d} ({v['BY_selected']} BY, {v['Holm_selected']} Holm)" for j, d, v in persel))

    macro('GRminBY', num(p180['min_BY']))
    macro('GRminHolm', num(p180['min_Holm']))
    macro('GRnBY', str(p180['BY_selected']))
    macro('GRnHolm', str(p180['Holm_selected']))

# ---------------- split ----------------
sp = R['split']
if sp.get('family_48'):
    rows = sp['rows']

    def gs(f, m, d, n='restricted'):
        return next(r for r in rows if (r['fold'], r['readout'], r['degree_mode'], r['null_mode']) == (f, m, d, n))
    lines = []
    for m in ORDER:
        cells = [SHORT[m]]
        for f in (1, 2):
            for d in ('none', 'own'):
                r = gs(f, m, d)
                cells.append(num(r['observed'], 5))
                cells.append(num(r['pIUT']))
        lines.append(' & '.join(cells) + r' \\')
    f48 = sp['family_48']
    prep = sp['prepared']
    write('t_split.tex', r"""\begin{table}[!htbp]
\centering\small
\caption{Disjoint discovery and evaluation samples. In each fold a panel selected from nine sample identifiers is evaluated on the other nine; values are incremental partial $\rho$ and restricted joint $p$ on the evaluation half's common support (9,999 transformations). Across all 48 hypotheses (two folds, six readouts, two degree modes, two transformation classes) the minimum raw joint $p$ is """ + num(f48['min_raw_pIUT']) + r""", the minimum BY value """ + num(f48['min_BY']) + r""" and the minimum Holm value """ + num(f48['min_Holm']) + r""" (selections: BY """ + str(f48['BY_selected']) + r""", Holm """ + str(f48['Holm_selected']) + r"""). GF, Geneformer; emb., embedding; attn., attention; BY, Benjamini--Yekutieli.}\label{tab:split}
\setlength{\tabcolsep}{3.4pt}
\begin{tabular}{@{}lrrrrrrrr@{}}\toprule
 & \multicolumn{4}{c}{Fold 1 (discovery A, evaluation B)} & \multicolumn{4}{c}{Fold 2 (discovery B, evaluation A)} \\\cmidrule(lr){2-5}\cmidrule(lr){6-9}
 & \multicolumn{2}{c}{No degree} & \multicolumn{2}{c}{Own degree} & \multicolumn{2}{c}{No degree} & \multicolumn{2}{c}{Own degree} \\\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}
Readout & $\rho$ & $p$ & $\rho$ & $p$ & $\rho$ & $p$ & $\rho$ & $p$ \\\midrule
""" + '\n'.join(lines) + r"""
\bottomrule\end{tabular}
\end{table}
""")
    macro('SPminp', num(f48['min_raw_pIUT']))
    macro('SPminBY', num(f48['min_BY']))
    macro('SPminHolm', num(f48['min_Holm']))
    macro('SPnBY', str(f48['BY_selected']))
    macro('SPnHolm', str(f48['Holm_selected']))
    macro('SPoverlapA', f"{prep['1']['overlap_with_registered_brain_panel']:,}")
    macro('SPoverlapB', f"{prep['2']['overlap_with_registered_brain_panel']:,}")
    macro('SPtfA', str(prep['1']['panel_tfs']))
    macro('SPtfB', str(prep['2']['panel_tfs']))
    agree = sum(1 for m in ORDER for d in ('none', 'own') if gs(1, m, d)['observed'] * gs(2, m, d)['observed'] > 0)
    macro('SPsignAgree', str(agree))
    absr = [abs(gs(f, m, d)['observed']) for f in (1, 2) for m in ORDER for d in ('none', 'own')]
    macro('SPabsMax', num_up(max(absr), 4))
    rnd = [abs(gs(f, 'random_init_geneformer', d)['observed']) for f in (1, 2) for d in ('none', 'own')]
    macro('SPrandMax', num_up(max(rnd), 4))
    def selt(fam):
        if fam['BY_selected'] == 0 and fam['Holm_selected'] == 0:
            return 'no hypothesis is selected by BY or Holm'
        return f"BY selects {fam['BY_selected']} and Holm {fam['Holm_selected']}: " + '; '.join(
            f"{SHORT[r['readout']]} (fold {r['fold']}, {r['degree_mode']} degree, {r['null_mode']}; $\\rho={mnum(r['observed'], 4)}$, Holm {num(r['p_Holm'])})" for r in fam['selected'])
    macro('SPselText', selt(f48))


# ---------------- bmmc ----------------
bm = R['bmmc']
if bm.get('family_150'):
    FIVE = ['geneformer_embed', 'random_init_geneformer', 'scgpt_encoder', 'scfoundation_encoder', 'uce_encoder']
    summ = {s['readout']: s for s in bm['donor_summary']}
    lines = []
    for m in FIVE:
        s = summ[m]
        di, pb = s['donor_interval'], s['paired_baseline_interval']
        lines.append(f"{SHORT[m]} & {num(di['point'], 5)} & [{num(di['ci95'][0], 5)}, {num(di['ci95'][1], 5)}] & "
                     f"{num(pb['point'], 5)} & [{num(pb['ci95'][0], 5)}, {num(pb['ci95'][1], 5)}] & "
                     f"{num(s['paired_baseline_signflip_p'])} \\\\")
    f150 = bm['family_150']
    write('t_bmmc.tex', r"""\begin{table}[!htbp]
\centering\small
\caption{Five bone-marrow evaluation donors on the common finite gene set (""" + f"{bm['common']['genes_common']:,}" + r""" genes, """ + str(bm['common']['tfs_common']) + r""" TFs). Donor means use the reference-degree restricted effect; the paired difference compares each model with co-expression under identical structural and reference-degree controls. Intervals are equal-donor percentile bootstrap intervals (2,000 resamples); the exact paired sign-flip test has a smallest attainable two-sided value of 0.0625 with five donors. Across all 150 donor $\times$ readout $\times$ degree $\times$ transformation hypotheses (9,999 draws) the minimum raw joint $p$ is """ + num(f150['min_raw_pIUT']) + r""", minimum BY """ + num(f150['min_BY']) + r""" and minimum Holm """ + num(f150['min_Holm']) + r""" (selections: BY """ + str(f150['BY_selected']) + r""", Holm """ + str(f150['Holm_selected']) + r"""). GF, Geneformer; emb., embedding; BY, Benjamini--Yekutieli.}\label{tab:bmmc5}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}lrrrrr@{}}\toprule
 & \multicolumn{2}{c}{Donor-mean $\rho$} & \multicolumn{2}{c}{Model minus co-expression} & Sign-flip \\\cmidrule(lr){2-3}\cmidrule(lr){4-5}
Readout & Mean & 95\% interval & Mean & 95\% interval & $p$ \\\midrule
""" + '\n'.join(lines) + r"""
\bottomrule\end{tabular}
\end{table}
""")
    macro('BMminp', num(f150['min_raw_pIUT']))
    macro('BMminBY', num(f150['min_BY']))
    macro('BMminHolm', num(f150['min_Holm']))
    macro('BMnBY', str(f150['BY_selected']))
    macro('BMnHolm', str(f150['Holm_selected']))
    macro('BMgenes', f"{bm['common']['genes_common']:,}")
    macro('BMtfs', str(bm['common']['tfs_common']))
    pts = [summ[m]['donor_interval']['point'] for m in FIVE]
    macro('BMmeanMin', num(min(pts), 4))
    macro('BMmeanMax', num(max(pts), 4))
    zero_in = sum(1 for m in FIVE if summ[m]['donor_interval']['ci95'][0] <= 0 <= summ[m]['donor_interval']['ci95'][1])
    macro('BMzeroIn', str(zero_in))
    dpts = [summ[m]['paired_baseline_interval']['point'] for m in FIVE]
    macro('BMdiffMin', num(min(dpts), 4))
    macro('BMdiffMax', num(max(dpts), 4))
    dneg = sum(1 for m in FIVE if summ[m]['paired_baseline_interval']['ci95'][1] < 0)
    dpos = sum(1 for m in FIVE if summ[m]['paired_baseline_interval']['ci95'][0] > 0)
    macro('BMdiffNeg', str(dneg))
    macro('BMdiffPos', str(dpos))
    macro('BMsignMin', num(min(summ[m]['paired_baseline_signflip_p'] for m in FIVE), 4))
    exc = [len(v) for v in bm['common']['excluded_by_donor'].values()]
    macro('BMexcMin', str(min(exc)))
    macro('BMexcMax', str(max(exc)))
    def selb(fam):
        if fam['BY_selected'] == 0 and fam['Holm_selected'] == 0:
            return 'no hypothesis is selected by BY or Holm'
        return f"BY selects {fam['BY_selected']} and Holm {fam['Holm_selected']}: " + '; '.join(
            f"{SHORT[r['readout']]} (donor {r['donor']}, {r['degree_mode']} degree, {r['null_mode']}; $\\rho={mnum(r['observed'], 4)}$, Holm {num(r['p_Holm'])})" for r in fam['selected'])
    macro('BMselText', selb(f150))


# ---------------- training exposure (Table 6) ----------------
scf = json.loads((W / 'evidence/scfoundation_verification.json').read_text())
scf_official = scf.get('official_sha256') is not None
if scf_official and not scf.get('official_file_matches_published_hash'):
    raise SystemExit("the downloaded scFoundation checkpoint does not match its published checksum; re-fetch it first")
scf_same = scf_official and bool(scf.get('local_identical_to_official_submodel'))
scf_part = ('the cell-embedding sub-model' if scf_same and scf['local_identical_to_official_submodel'][0] == 'cell'
            else 'a sub-model')
if scf_same:
    n = scf['comparisons'][scf['local_identical_to_official_submodel'][0]]['exactly_equal']
    scf_weight = (f"Evaluated file equals the published checksum of its source mirror; all {n} tensors equal "
                  f"{scf_part} of the developers' published checkpoint, whose checksum was confirmed.")
    scf_limit = 'Exact-cell membership is not verified.'
elif scf_official:
    print('WARNING: the evaluated scFoundation weights are not tensor-identical to the developers\' checkpoint')
    scf_weight = ("Evaluated file equals the published checksum of its source mirror; no sub-model of the developers' "
                  "published checkpoint has identical tensors.")
    scf_limit = "Identity with the developers' release is not established; exact-cell membership is not verified."
elif scf.get('local_equals_mirror_file'):
    scf_weight = ("The evaluated file is byte-identical (SHA-256) to the published model file of a third-party mirror "
                  "distributed as the scFoundation cell-embedding model.")
    scf_limit = "Tensor identity with the developers' own checkpoint and exact-cell membership are not verified."
else:
    scf_weight = None
if scf_weight:
    write('t_exposure.tex', r"""\begin{table}[!htbp]
\centering\small
\caption{Training exposure of the evaluated checkpoints to the genuine-brain source study. Study and sample inventory inclusion is distinct from weight identity and exact-cell membership. No entry establishes a never-pretrained holdout; the random-initialized Geneformer readout is unexposed by construction.}\label{tab:exposure}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}p{.17\linewidth}p{.47\linewidth}p{.30\linewidth}@{}}\toprule
Model & Positive inventory and weight evidence & Remaining limit \\\midrule
Geneformer V2 & Genecorpus-104M lists the GSE174367 publication and 61,770 cells. Evaluated weights and dictionaries match the recorded repository revision. & Exact membership of the evaluated nuclei is not verified. \\
scFoundation & All 19 source RNA GSMs, including all 18 sample identifiers used here and Sample-101; 45,508 retained cells. """ + scf_weight + r""" & """ + scf_limit + r""" \\
scGPT & No complete checkpoint-specific exclusion inventory established. & Study, sample and cell overlap cannot be excluded. \\
UCE (four-layer) & No complete checkpoint-specific exclusion inventory established. & Study, sample and cell overlap cannot be excluded. \\
\bottomrule\end{tabular}
\end{table}
""")

# ---------------- abstract sentence ----------------
if all(k in macros for k in ('CSminHolm', 'GRminHolm', 'SPminHolm', 'BMminHolm')):
    fam_sel = {'common': (cs['six']['family']['BY_selected'], cs['six']['family']['Holm_selected']),
               'grid': (gr['pooled_180']['BY_selected'], gr['pooled_180']['Holm_selected']),
               'split': (sp['family_48']['BY_selected'], sp['family_48']['Holm_selected']),
               'bmmc': (bm['family_150']['BY_selected'], bm['family_150']['Holm_selected'])}
    positives = []
    for key, fam in (('common', cs['six']['family']), ('grid', gr['pooled_180']), ('split', sp['family_48']), ('bmmc', bm['family_150'])):
        positives += [(key, r) for r in fam['selected'] if r['observed'] > 0]
    if sum(a + b for a, b in fam_sel.values()) == 0:
        macro('AbsMain', 'Across 72 common-support hypotheses with 9,999 permutations, no readout showed alignment '
              'beyond co-expression and structural covariates (minimum Holm-adjusted $p=' + macros['CSminHolm'] + '$), '
              'and no selection arose under ten degree designs, with gene panels selected on disjoint samples, '
              'or across five bone-marrow donors.')
    else:
        macro('AbsMain', 'Across 72 common-support hypotheses with 9,999 permutations the minimum Holm-adjusted $p$ was '
              + macros['CSminHolm'] + '; across all new families ' + str(sum(b for a, b in fam_sel.values())) +
              ' Holm and ' + str(sum(a for a, b in fam_sel.values())) + ' BY selections occurred, ' +
              (str(len(positives)) + ' of them positive.' if positives else 'none of them positive.'))

if scf_same:
    macro('ScfDA', f"The evaluated scFoundation weights are tensor-identical to {scf_part} of the developers' published "
          "checkpoint, and byte-identical to the model file of the third-party mirror from which they were obtained.")
elif scf_official:
    macro('ScfDA', "The evaluated scFoundation weights are byte-identical to the model file of a third-party mirror "
          "distributed as the scFoundation cell-embedding model but are not tensor-identical to any sub-model of the "
          "developers' published checkpoint.")
else:
    macro('ScfDA', "The evaluated scFoundation weights are byte-identical to the model file of a third-party mirror distributed as the scFoundation cell-embedding model; tensor identity with the developers' own checkpoint has not been verified.")

write('r3_numbers.tex', '% result macros\n' +
      '\n'.join(f'\\newcommand{{\\{k}}}{{{v}}}' for k, v in sorted(macros.items())) + '\n')
print('tables:', sorted(p.name for p in OUT.glob('t_*.tex')), '| macros:', len(macros))
