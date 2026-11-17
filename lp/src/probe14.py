"""Isolated root probe: never writes the campaign ledger; import only after exact replay."""
import argparse,json,fcntl,re,signal
from pathlib import Path
import coverage14 as c
def request_stop(*_):
    c.stop=True
signal.signal(signal.SIGTERM,request_stop)
signal.signal(signal.SIGINT,request_stop)
ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True,choices=[6,7]);ap.add_argument('--config',required=True,choices=['515','617']);ap.add_argument('--seconds',type=int,default=3600);ap.add_argument('--nonlinear',action='store_true');ap.add_argument('--cut-scope',choices=['all','core','pair','coordinate'],default='all');ap.add_argument('--nonlinear-scope',choices=['all','pair','singletons'],default='all');ap.add_argument('--tag',default='');ap.add_argument('--prefix',default='');ap.add_argument('--maximal-only',action='store_true');ap.add_argument('--solver',choices=['ipm','simplex'],default='ipm');ap.add_argument('--pair-bounds',action='store_true');ap.add_argument('--degree-hull',action='store_true');ap.add_argument('--triple-scope',choices=['none','all','pair','singletons'],default='none');ap.add_argument('--parallel-probe',action='store_true');ap.add_argument('--compact-objectives',action='store_true');ap.add_argument('--point-pairs',action='store_true');ap.add_argument('--plane-budget',action='store_true');ap.add_argument('--point-triples',action='store_true');ap.add_argument('--star-bounds',action='store_true');ap.add_argument('--rounded-budget',action='store_true');ap.add_argument('--nonfano',action='store_true');ap.add_argument('--extra',action='store_true');ap.add_argument('--rich',action='store_true');ap.add_argument('--point-mix',action='store_true');ap.add_argument('--round-degree',action='store_true');ap.add_argument('--point-subset',action='store_true');a=ap.parse_args();prefix=tuple(map(int,a.prefix.split(','))) if a.prefix else ();suffix=('_'+a.tag) if a.tag else ('_p'+'_'.join(map(str,prefix)) if prefix else '')
assert not a.tag or re.fullmatch('[A-Za-z0-9_]+',a.tag)
assert not a.parallel_probe or a.tag
c.CERT_SUFFIX=('_'+a.tag) if a.tag else ''
lock=open(f'out/probe14_D{a.D}_{a.config}'+(suffix if a.parallel_probe else '')+'.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
ledger=json.load(c.path(a.D).open());node=c.nid(a.config,prefix)
assert node in ledger['nodes']
assert ledger['nodes'][node]['status'] not in ('LP-certified','joint-LP-certified','exact-primal')
output=Path(f'out/probe14_D{a.D}_{a.config}{suffix}.json');c.path=lambda D:output
c.solve_node(ledger,node,a.seconds,a.solver,a.nonlinear,a.cut_scope,a.nonlinear_scope,a.maximal_only,a.pair_bounds,a.degree_hull,a.triple_scope,a.compact_objectives,a.point_pairs,a.plane_budget,a.point_triples,a.star_bounds,a.rounded_budget,a.nonfano,a.extra,a.rich,a.point_mix,a.round_degree,a.point_subset);c.check_tree(ledger);c.refresh(ledger)
print('PROBE14_COMPLETE',a.D,node,ledger['nodes'][node]['status'],flush=True)
