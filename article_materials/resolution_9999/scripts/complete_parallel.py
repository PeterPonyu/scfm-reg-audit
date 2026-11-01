#!/usr/bin/env python3
"""Scheduling-only continuation of disjoint fixed streams on the same remote host.

The original runner and statistical implementation remain unchanged. An interrupted
sequential runner must have exited before launch. Each worker owns one group directory;
existing checkpoint arrays and RNG states are retained, including any extension draws.
"""
import concurrent.futures
import hashlib
import os
from pathlib import Path
import time

import numpy as np
from run_resolution import inputs
from screg_validation.batched_null import batched_components
from screg_validation.common import digest, load_json, sha256, write_json


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'inputs'
OUTPUT = ROOT / 'output'


def execute(group):
    started = time.monotonic()
    plan = load_json(INPUT / 'posthoc-plan.json')
    for name, expected in plan['input_sha256'].items():
        assert sha256(INPUT / name) == expected, name
    pilot = load_json(OUTPUT / 'pilot.json')
    assert pilot['passed'] and pilot['plan_sha256'] == digest(plan)
    a, tf, ii, jj, structural, controls, names, scores = inputs(INPUT)
    degree, null = group.split('-')
    strata = np.zeros(len(a['proxy']), int)
    if null == 'restricted':
        ob = np.digitize(a['peak_count'], np.unique(np.quantile(a['peak_count'], [.2, .4, .6, .8])))
        gb = np.digitize(a['gc'], np.unique(np.quantile(a['gc'], [1/3, 2/3])))
        strata = 15 * np.isin(np.arange(len(a['proxy'])), tf) + 3 * ob + gb
    prefix = load_json(INPUT / ('prefix-' + group + '.json'))
    context = {'plan_sha256': digest(plan), 'group': group, 'readouts': names,
               'prefix_sha256': sha256(INPUT / ('prefix-' + group + '.json'))}
    target = OUTPUT / group
    target.mkdir(exist_ok=True)
    checkpoint = target / 'draw-progress.json'
    if not checkpoint.exists():
        state = {k: v for k, v in prefix.items() if k != 'state_sha256'}
        state['context_sha256'] = digest(context)
        write_json(checkpoint, {**state, 'state_sha256': digest(state)})
    before = load_json(checkpoint)
    assert before['context_sha256'] == digest(context)
    assert before['state_sha256'] == digest({k: v for k, v in before.items() if k != 'state_sha256'})
    assert before['next_draw'] < 9999 and not (target / 'rows.json').exists()
    write_json(target / 'parallel-start-checkpoint.json', before)
    last = [0.0]
    def progress():
        elapsed = time.monotonic() - started
        if elapsed > 21600:
            raise RuntimeError('Fixed continuation wall-clock budget exceeded')
        if elapsed - last[0] > 60:
            print('PROGRESS', group, load_json(checkpoint)['next_draw'], 'of 9999;',
                  round(elapsed, 1), 'seconds', flush=True)
            last[0] = elapsed
    backend = 'cuda' if degree == 'own' else 'cpu'
    rows = batched_components(scores, a['proxy'], tf, controls, strata=strata,
                              seed=20261001, degree_mode=degree, permutations=9999,
                              backend=backend, check=progress, checkpoint=checkpoint,
                              checkpoint_context=context)
    after = load_json(checkpoint)
    assert after['next_draw'] == 9999
    for component in ['M', 'D']:
        assert after['nulls'][component][:before['next_draw']] == before['nulls'][component]
        assert after['nulls'][component][:999] == prefix['nulls'][component]
    for name, row in zip(names, rows):
        row.pop('nulls')
        row.update(readout=name, null_mode=null)
    write_json(target / 'rows.json', rows)
    receipt = {'group': group, 'start_draw': before['next_draw'], 'end_draw': 9999,
               'all_preexisting_draws_preserved_exactly': True, 'seconds': time.monotonic() - started,
               'pid': os.getpid(), 'backend': backend,
               'starting_checkpoint_sha256': sha256(target / 'parallel-start-checkpoint.json'),
               'final_checkpoint_sha256': sha256(checkpoint)}
    print('FINISHED', group, receipt['seconds'], 'seconds', flush=True)
    return receipt


if __name__ == '__main__':
    started = time.time()
    groups = [degree + '-' + null for degree in ['own', 'none']
              for null in ['restricted', 'unrestricted']
              if not (OUTPUT / (degree + '-' + null) / 'rows.json').exists()]
    assert groups, 'Nothing remains; do not overwrite a completed run'
    print('CONTINUE', groups, 'with at most three disjoint stream workers', flush=True)
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        receipts = list(pool.map(execute, groups))
    write_json(OUTPUT / 'parallel-execution.json', {
        'scheduling_only_change': True, 'statistical_implementation_unchanged': True,
        'reason': 'CPU-dominated execution. The existing one-GPU-process guard is retained: the own-degree continuation uses CUDA, and degree-free continuations use the CPU backend validated by the original portability pilot.',
        'plan_sha256': digest(load_json(INPUT / 'posthoc-plan.json')),
        'scheduler_sha256': sha256(Path(__file__)), 'original_runner_sha256': sha256(ROOT / 'run_resolution.py'),
        'start_unix': started, 'finish_unix': time.time(), 'groups': receipts,
        'no_adaptive_draw_count_or_data_dependent_stopping': True})
    print('ALL STREAMS COMPLETE; original runner can assemble the fixed families.', flush=True)
