"""Bind completed fresh replays into a proof-chain manifest (no solver).

This checks replay records and their dependencies; the full arithmetic replay is
coverage14.py --replay --fresh-prices, not this manifest check.
"""
import hashlib
import json
from pathlib import Path
from lift11 import atomic
from maximal14 import guard


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    records = []
    guard5 = read('out/coverage14_D5_guard.json')
    assert sha(guard5['path']) == guard5['sha256']
    d5 = read('out/replay14_D5_corrected.json')
    assert d5['status'] == 'EXACT_D5_CORRECTED_TABLE_REPLAYED'
    assert d5.get('fresh_prices') is True
    assert d5['ledger_sha256'] == guard5['sha256'] and d5['ledger_untouched']
    assert d5['projection_table_sha256'] == sha('out/union_tables_q5.json')
    ledger5 = read(guard5['path'])
    assert ledger5['counts'] == dict(pruned_certified=5198, surviving=0, open=0)
    assert len(d5['results']) == 2
    for item in d5['results']:
        assert item['detail']['rhs'] < 0 and item['detail']['min_column'] >= 0
        assert sha(item['certificate']) == item['certificate_sha256']
    records.append(dict(D=5, cases=5198, ledger=guard5['path'],
                        ledger_sha256=guard5['sha256'],
                        replay='out/replay14_D5_corrected.json'))
    for D, total in ((6, 51244), (7, 274418)):
        path = f'out/coverage14_D{D}.json'
        ledger = read(path)
        replay_path = f'out/replay14_D{D}_final.json'
        replay = read(replay_path)
        expected = dict(pruned_certified=total, surviving=0, open=0)
        assert ledger['counts'] == replay['counts'] == expected
        assert replay['status'] == 'EXACT_COVERAGE14_REPLAYED'
        assert replay['fresh_prices'] is True and replay['tree_cases'] == total
        assert replay['ledger_sha256'] == sha(path)
        assert len({x['node'] for x in replay['certificates']}) == len(replay['certificates'])
        assert sum(x['cases'] for x in replay['certificates']) == total
        for item in replay['certificates']:
            assert sha(item['certificate']) == item['sha256']
            assert item['rhs'] < 0 and item['min_column'] >= 0
            assert ledger['nodes'][item['node']]['certificate'] == item['certificate']
            cert = read(item['certificate'])
            if cert.get('domain'):
                from domain14 import check_record,membership
                check_record(cert)
                assert membership(D,ledger['nodes'][item['node']]['prefix'],cert['domain'])==item['cases']
            assert cert['projection_table_sha256'] == sha('out/union_tables_q5.json')
            if cert.get('maximal_only'):
                assert cert['reduction_sha256'] == guard(D, cert['config'])
            if cert.get('star_bounds'):
                from star14 import guard as star_guard
                assert cert['star_sha256'] == star_guard(D, cert['config'])
            if cert.get('rounded_budget'):
                from roundbudget14 import guard as rounded_guard
                assert cert['rounded_budget_sha256'] == rounded_guard()
        dp = read(f'out/dp14_D{D}.json')
        assert dp['status'] == 'EXACT_DP14_VERIFIED' and dp['mismatches'] == 0
        assert dp['types'] == {6: 119104, 7: 224438}[D]
        records.append(dict(D=D, cases=total, ledger=path, ledger_sha256=sha(path),
                            replay=replay_path, replay_sha256=sha(replay_path)))
    reduction = read('out/reduction14.json')
    assert reduction['status'] == 'EXACT_REDUCTION14_VERIFIED'
    assert sum(x['verified'] for x in reduction['results']) == 82
    assert reduction['projection_table_sha256'] == sha('out/union_tables_q5.json')
    source_files = sorted(set(Path('.').glob('*.py')) | set(Path('.').glob('*.cpp'))
                          | set(Path('.').glob('*.sh')))
    evidence_files = set()
    for pattern in ('*14_D[67]_prices.json', '*14*validation*.json', 'validation14*.json',
                    'reduction14*.json', 'dp14_D[67].json', '*14_D[67]_independent.json',
                    'point14_projective_q*.json', 'nonfano14_projective_q*.json', 'rich14_projective_q*.json'):
        evidence_files.update(Path('out').glob(pattern))
    for pattern in ('union14_*.json', 'geometry14_D[67].json', 'types14_D[67].json', 'raw14_D[67].json'):
        evidence_files.update(Path('out').glob(pattern))
    for filename in ('out/lift11_prices.json','out/lift13_pair_prices.json',
                     'out/lift14_D5_fresh_check.json','out/lift13_pair_fresh_check.json',
                     'out/replay14_D5_corrected.json','out/coverage14_D5_guard.json',
                     'out/union_tables_q5.json','out/degree_independent14_tables.json',
                     'out/cert14_lower46.log',
                     'solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json',
                     'solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        assert Path(filename).is_file();evidence_files.add(Path(filename))
    assert 'certificate OK' in Path('out/cert14_lower46.log').read_text()
    assert read('out/lift14_D5_fresh_check.json')['mismatches']==0
    assert read('out/lift13_pair_fresh_check.json')['mismatches']==0
    geometric=read('out/geometric_inputs14_check.json')
    assert geometric['status']=='EXACT_GEOMETRIC_INPUTS14_REPLAYED' and geometric['line_systems']==46656
    evidence_files.add(Path('out/geometric_inputs14_check.json'))
    for filename,expected_sha in geometric['input_sha256'].items():
        assert sha(filename)==expected_sha;evidence_files.add(Path(filename))
    from geom import AG
    from planecount_check import witness
    upper_witnesses=[]
    ag=AG(3,5)
    for filename in ('solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=read(filename);assert (data['n'],data['q'])==(3,5)
        points={ag.pidx[tuple(p)] for p in data['points']};assert len(points)==53
        lines=witness(ag,points);assert len(lines)==31
        upper_witnesses.append(dict(file=filename,sha256=sha(filename),points=53,witness_lines=lines))
    record = dict(status='COVERAGE14_PROOF_CHAIN_MANIFEST_CHECKED',
                  method='Binds completed full fresh arithmetic replays; not a substitute for them',
                  q=5, excluded_n=46, bracket=[47, 53],
                  total_cases=sum(x['cases'] for x in records), degrees=records,
                  upper_witnesses=upper_witnesses,
                  reduction_certificates=82,
                  reduction_sha256=sha('out/reduction14.json'),
                  projection_table_sha256=sha('out/union_tables_q5.json'),
                  sources={str(p): sha(p) for p in source_files},
                  evidence_sha256={str(p): sha(p) for p in sorted(evidence_files)})
    atomic('out/final14_manifest.json', record)
    print('COVERAGE14_PROOF_CHAIN_MANIFEST_CHECKED', record['total_cases'], flush=True)


if __name__ == '__main__':
    main()
