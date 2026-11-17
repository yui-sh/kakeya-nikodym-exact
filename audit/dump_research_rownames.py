"""NOT part of the independent checker. Bookkeeping helper: the D=5 root certificates are positional multiplier lists,
so this dumps the research model's row list (name, sense, rhs) once, for the position -> row-name map used by
check_root.py. It imports lp/src (research code) on purpose and is the only file in audit/ that does so.
Run from lp/src with the research venv:  cd lp/src && .venv/bin/python ../../audit/dump_research_rownames.py 515|617"""
import os, sys, json, collections
os.environ.update(PC6_SECTIONS='1', PC6_EXACT_K='4', PC6_F6_RMIN='5', PC6_DMH='1', PC6_UNMARKED='1')
sys.path.insert(0, '.')
from coverage11 import build, node_lp
cfg = sys.argv[1]; scope = 'core' if cfg == '515' else 'pair'
base, info, types = build(scope, cfg, []); lp = node_lp(base, types, cfg, [])
rows = [(name, s, int(b), len(d)) for (d, b, s, name) in lp.rows]
pref = collections.Counter(''.join(ch for ch in nm if not ch.isdigit() and ch not in ',') for nm in lp.names)
json.dump(dict(config=cfg, scope=scope, ncols=len(lp.names), col_prefix_counts=dict(pref), ntypes=len(types),
               rows=rows), open(f'../../audit/out/research_rownames_{cfg}.json', 'w'))
print('DUMPED', cfg, len(lp.names), len(lp.rows), dict(pref), flush=True)
