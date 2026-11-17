"""NOT part of the independent checker (same role as dump_research_rownames.py, for the D = 6/7 session-XIV roots): dumps the
research model's row list (name, sense, rhs, nnz) and column-family counts for the model a given certificate was built on,
so that the positional multipliers can be mapped to row names. Imports lp/src on purpose.
Run from lp/src:  .venv/bin/python ../../audit/dump_research_rownames14.py <cert.json> <tag>"""
import os, sys, json, collections
os.environ.update(PC6_SECTIONS='1', PC6_EXACT_K='4', PC6_F6_RMIN='5', PC6_DMH='1', PC6_UNMARKED='1')
sys.path.insert(0, '.')
from coverage14 import build
cert = json.load(open(sys.argv[1])); tag = sys.argv[2]
flags = dict(nonlinear=cert.get('nonlinear', False), cut_scope=cert.get('cut_scope', 'all'), nonlinear_scope=cert.get('nonlinear_scope', 'all'),
             maximal_only=cert.get('maximal_only', False), pair_bounds=cert.get('pair_bounds', False), degree_hull=cert.get('degree_hull', False),
             triple_scope=cert.get('triple_scope', 'none'), compact_objectives=cert.get('compact_objectives', False), point_pairs=cert.get('point_pairs', False),
             plane_budget=cert.get('plane_budget', False), point_triples=cert.get('point_triples', False), star_bounds=cert.get('star_bounds', False),
             rounded_budget=cert.get('rounded_budget', False), nonfano=cert.get('nonfano', False), extra=cert.get('extra', False), rich=cert.get('rich', False),
             point_mix=cert.get('point_mix', False), round_degree=cert.get('round_degree', False), point_subset=cert.get('point_subset', False))
if 'domain' in cert:   # domain certificate: restricted raw-vector universe + explicit domain row, built by the research domain14 module
    from domain14 import build as dbuild, check_record
    flags = check_record(cert); lp, info, types = dbuild(cert['D'], cert['config'], cert['domain'], flags)
else:
    lp, info, types = build(cert['D'], cert['config'], tuple(cert['prefix']), **flags)
rows = [(name, s, int(b), len(d)) for (d, b, s, name) in lp.rows]
pref = collections.Counter(''.join(ch for ch in nm.split(':')[0] if not ch.isdigit() and ch not in ',') for nm in lp.names)
json.dump(dict(cert=sys.argv[1], D=cert['D'], config=cert['config'], prefix=cert['prefix'], flags=flags, ncols=len(lp.names), col_prefix_counts=dict(pref),
               ntypes=len(types), nmult=len(cert['multipliers']), rows=rows), open(f'../../audit/out/research_rownames_{tag}.json', 'w'))
print('DUMPED', tag, len(lp.names), len(lp.rows), 'multipliers', len(cert['multipliers']), dict(pref), flush=True)
