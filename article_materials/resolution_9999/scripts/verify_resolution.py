#!/usr/bin/env python3
"""Recount saved component tails and independently reconstruct BY/Holm adjustments.

Uses Python's standard library only; does not import the experiment implementation.
The supplied root contains the completed sidecar and original/ with its bound plan,
original 72-row results and four original 999-draw prefixes.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def correction(probabilities, method):
    """Explicit ordered-test definitions, independent of vectorized runner code."""
    m = len(probabilities)
    order = sorted(range(m), key=lambda i: probabilities[i])
    harmonic = math.fsum(1 / k for k in range(1, m + 1))
    result = [None] * m
    for rank, index in enumerate(order):
        if method == 'BY':
            value = min(probabilities[order[j]] * m * harmonic / (j + 1)
                        for j in range(rank, m))
        elif method == 'Holm':
            value = max(probabilities[order[j]] * (m - j)
                        for j in range(rank + 1))
        else:
            raise ValueError(method)
        result[index] = min(1.0, value)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--out', type=Path, help='folder for the report files; without it the report is only printed')
    args = p.parse_args()
    root = args.root
    original = root / 'original'
    # Files edited for publication keep their run-time digests in publication_edits.json.
    edits_file = root.resolve().parents[1] / 'publication_edits.json'
    prefix = 'resolution_9999/results/'
    edits = ({k[len(prefix):]: v for k, v in read(edits_file).items() if k.startswith(prefix)}
             if edits_file.exists() else {})
    plan = read(original / 'posthoc-plan.json')
    plan_digest = edits.get('original/posthoc-plan.json', {}).get('original_digest') or digest(plan)
    old = read(original / 'results.json')
    current = read(root / 'results.json')
    done = read(root / 'completion.json')
    pilot = read(root / 'pilot.json')
    assert current['status'] == 'complete' and current['permutations'] == 9999
    assert current['post_hoc'] is True
    assert done['all_four_groups_complete'] and done['prefix_preserved']
    assert pilot['passed']
    assert len(pilot['groups']) == 4
    pilot_errors = [group['direct_projection_max_abs'] for group in pilot['groups']]
    for group in pilot['groups']:
        for backend in ['cpu', 'cuda']:
            pilot_errors.extend(group['checks'][backend].values())
    assert max(pilot_errors) < 1e-10
    assert all(x['plan_sha256'] == plan_digest for x in [current, done, pilot])
    assert done['results_sha256'] in (sha(root / 'results.json'), edits.get('results.json', {}).get('original_sha256'))
    for filename, expected in done['outputs'].items():
        assert expected in (sha(root / filename), edits.get(filename, {}).get('original_sha256')), filename
    for path in original.iterdir():
        if path.name in plan['input_sha256']:
            assert plan['input_sha256'][path.name] in (sha(path), edits.get('original/' + path.name, {}).get('original_sha256')), path
    assert len(current['rows']) == len(old['rows']) == 72
    measured = [r for r in current['rows'] if r['measured']]
    assert len(measured) == 16
    keys = ['job', 'readout', 'degree_mode', 'null_mode']
    def identity(row):
        return tuple(row[k] for k in keys)
    old_by_key = {identity(row): row for row in old['rows']}
    assert len(old_by_key) == 72
    assert set(old_by_key) == {identity(row) for row in current['rows']}
    readouts = ['geneformer_embed', 'random_init_geneformer', 'scgpt_encoder', 'uce_encoder']
    counted = []
    groups = []
    tail_boundary_distances = []
    for degree in ['own', 'none']:
        for null in ['restricted', 'unrestricted']:
            group = degree + '-' + null
            prefix = read(original / ('prefix-' + group + '.json'))
            saved = read(root / group / 'draw-progress.json')
            assert prefix['next_draw'] == 999 and saved['next_draw'] == 9999
            for state in [prefix, saved]:
                assert state['state_sha256'] == digest({k: v for k, v in state.items() if k != 'state_sha256'})
            context = {'plan_sha256': plan_digest, 'group': group, 'readouts': readouts,
                       'prefix_sha256': sha(original / ('prefix-' + group + '.json'))}
            assert saved['context_sha256'] == digest(context)
            for component in ['M', 'D']:
                assert saved['nulls'][component][:999] == prefix['nulls'][component]
                assert len(saved['nulls'][component]) == 9999
                assert all(len(row) == 4 and all(math.isfinite(x) for x in row)
                           for row in saved['nulls'][component])
            resumed = root / group / 'parallel-start-checkpoint.json'
            if resumed.exists():
                at_switch = read(resumed)
                assert at_switch['state_sha256'] == digest({k: v for k, v in at_switch.items() if k != 'state_sha256'})
                for component in ['M', 'D']:
                    assert saved['nulls'][component][:at_switch['next_draw']] == at_switch['nulls'][component]
            for i, name in enumerate(readouts):
                row = next(r for r in measured if r['degree_mode'] == degree
                           and r['null_mode'] == null and r['readout'] == name)
                before = old_by_key[identity(row)]
                assert row['permutations'] == 9999
                assert abs(row['observed'] - before['observed']) < 1e-10
                for field in ['pM', 'pD', 'pIUT', 'q_IUT_BY', 'p_IUT_Holm', 'observed']:
                    assert row['registered_999'][field] == before[field]
                counts = {}
                for component in ['M', 'D']:
                    tail_boundary_distances.append(min(abs(abs(x[i]) - abs(row['observed'])) for x in saved['nulls'][component]))
                    counts[component] = sum(abs(x[i]) >= abs(row['observed'])
                                            for x in saved['nulls'][component])
                    probability = (counts[component] + 1) / 10000
                    assert probability == row['p' + component], (group, name, component)
                    old_count = sum(abs(x[i]) >= abs(before['observed'])
                                    for x in prefix['nulls'][component])
                    assert (old_count + 1) / 1000 == before['p' + component]
                assert row['pIUT'] == max(row['pM'], row['pD'])
                assert row['correction_input_p'] == row['pIUT']
                counted.append({**{k: row[k] for k in keys}, 'observed': row['observed'],
                                'tail_M_9999': counts['M'], 'tail_D_9999': counts['D'],
                                'registered_pIUT_999': before['pIUT'],
                                **{k: row[k] for k in ['pM', 'pD', 'pIUT', 'q_IUT_BY',
                                     'p_IUT_Holm', 'q_IUT_BY_measured16', 'p_IUT_Holm_measured16']}})
            groups.append({'group': group, 'draws_per_component': 9999,
                           'readouts': 4, 'original_prefix_exact': True,
                           'continuation_backends': sorted({r['backend'] for r in measured
                               if r['degree_mode'] == degree and r['null_mode'] == null})})
    for row in current['rows']:
        if not row['measured']:
            assert row == old_by_key[identity(row)], identity(row)
            assert row['correction_input_p'] == 1
    summaries = {}
    for family, rows, suffix in [('registered72', current['rows'], ''), ('measured16', measured, '_measured16')]:
        probabilities = [r['correction_input_p'] for r in rows]
        summaries[family] = {}
        for method, field in [('BY', 'q_IUT_BY'), ('Holm', 'p_IUT_Holm')]:
            independent = correction(probabilities, method)
            error = max(abs(x - r[field + suffix]) for x, r in zip(independent, rows))
            assert error < 1e-14, (family, method, error)
            summaries[family][method] = {'minimum': min(independent),
                                         'selected_below_0_05': sum(x < .05 for x in independent),
                                         'independent_max_abs_difference': error}
    old_measured = [r for r in old['rows'] if r['measured']]
    old_exploratory = {method: min(correction([r['pIUT'] for r in old_measured], method))
                       for method in ['BY', 'Holm']}
    assert old_exploratory == {'BY': 1.0, 'Holm': .432}
    report = {'passed': True, 'verification': 'Independent standard-library tail recount and ordered correction definitions',
              'groups': groups, 'component_probabilities_recounted': 32, 'registered_slots': 72,
              'measured_slots': 16, 'unavailable_rows_unchanged': 56,
              'total_saved_component_statistics': 4 * 2 * 9999 * 4,
              'minimum_IUT_9999': min(r['pIUT'] for r in measured), 'adjustments': summaries,
              'registered_999_measured_only_minima': old_exploratory,
              'observed_max_abs_difference': max(abs(r['observed'] - old_by_key[identity(r)]['observed']) for r in measured),
              'observed_absolute_max': max(abs(r['observed']) for r in measured),
              'minimum_distance_from_any_tail_boundary': min(tail_boundary_distances),
              'portability_pilot_max_abs_error': max(pilot_errors),
              'results_sha256': sha(root / 'results.json')}
    if args.out:  # the checked folder itself is never written
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / 'independent-verification.json').write_text(json.dumps(report, indent=2) + '\n')
        with (args.out / 'resolution-comparison.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(counted[0]))
            writer.writeheader()
            writer.writerows(counted)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
