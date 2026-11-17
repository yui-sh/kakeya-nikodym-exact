"""Exact parallel-class rounding linked to the complete raw degree-vector hull.

For each degree subset S and t>=2, five parallel sections partition its M
points. Hence max(0,ceil((M-5(t-1))/t)) <= sum floor(h_i/t) <= floor(M/t).
Sum over all 31 classes. Local plane costs depend only on h, so their prices
are exact constants on every assignment of every admitted full type.
"""
import json,hashlib,argparse
from collections import Counter,defaultdict
from pathlib import Path
from functools import lru_cache
from geometry14 import freeze
from lift11 import atomic,key
FAMILY='parallel-degree-subset-rounding-v1'
def subsets(D):return sorted({(j,) for j in range(1,D+1)}|{tuple(range(1,j+1)) for j in range(2,D)})
def lower(M,t,q):return max(0,-(-(M-q*(t-1))//t))
def upper(M,t,q):return M//t
@lru_cache(None)
def costs(h,D):return tuple(sum(h[j-1] for j in s)//t for s in subsets(D) for t in range(2,26))
def table(D,fresh=False):
    # Formula prices are recomputed, never loaded from a filtered catalogue.
    if fresh:costs.cache_clear()
    path=Path(f'out/rounddegree14_D{D}_prices.json');stored=json.load(open(f'out/types14_D{D}.json'));digest=hashlib.sha256();hist=set()
    for raw in stored['types']:
        t=freeze(raw);h=t[2];vals=costs(h,D);hist.add(h)
        digest.update(key((t,vals)).encode());digest.update(b'\n')
    result=dict(D=D,family=FAMILY,types=len(stored['types']),histograms=len(hist),prices=len(stored['types'])*len(subsets(D))*24,sha256=digest.hexdigest(),formula='mu_a(S,t)=floor(sum_{j in S} h_j/t), constant on every admitted assignment',subsets=subsets(D))
    if fresh and path.exists():assert json.load(path.open())==json.loads(key(result))
    atomic(path,result);print('ROUNDDEGREE14_ALL_TYPE_PRICES',result,flush=True);return result

def validate(D):
    from geom import AG
    from structure import planes_of
    from planecount_check import witness
    records=[]
    for fn in ('solutions/best_n3_q3_13.json','solutions/best_n3_q4_28_dvir.json','solutions/best_n3_q5_53.json','solutions/bc_n3_q5_53.json'):
        data=json.load(open(fn));ag=AG(3,data['q']);q=ag.q;S={ag.pidx[tuple(p)] for p in data['points']};W=witness(ag,S);degree=Counter(p for l in W for p in ag.lines[l][1]);m=Counter(degree[p] for p in S);classes=defaultdict(list)
        for normal,offset,ps in planes_of(ag):classes[key(normal)].append(tuple(sum(degree[p]==j for p in S&ps) for j in range(1,D+1)))
        assert len(classes)==q*q+q+1 and all(len(v)==q for v in classes.values())
        checks=0
        for subset in subsets(D):
            M=sum(m[j] for j in subset)
            for rows in classes.values():
                hs=[sum(h[j-1] for j in subset) for h in rows];assert sum(hs)==M
                for t in range(2,26):
                    v=sum(h//t for h in hs);assert lower(M,t,q)<=v<=upper(M,t,q);checks+=2
        from planecount5 import LP
        lp=LP();actual_h=[h for rows in classes.values() for h in rows]
        types=[(0,sum(h),h,0,(),None) for h in actual_h]
        for i in range(len(types)):lp.var(f'N{i}')
        mv=tuple(m[j] for j in range(1,D+1));qn='Q14:'+','.join(map(str,mv));lp.var(qn)
        add_cuts(lp,{},types,D,[mv],False,q)
        values=[1]*len(lp.names)
        for row,b,sense,name in lp.rows:
            lhs=sum(c*values[j] for j,c in row.items());assert lhs<=b if sense=='<' else lhs==b,(fn,name,lhs,b)
        rec=dict(file=fn,class_bound_checks=checks,generated_model_checks=len(lp.rows),violations=0);records.append(rec);print('ROUNDDEGREE14_VALIDATE',rec,flush=True)
    # Exhaust all five-part allocations through total 18, independently of geometry.
    checks=0
    def rec(xs,left):
        nonlocal checks
        if len(xs)==4:
            parts=xs+[left];M=sum(parts)
            for t in range(2,10):
                v=sum(x//t for x in parts);assert lower(M,t,5)<=v<=upper(M,t,5);checks+=2
            return
        for x in range(left+1):rec(xs+[x],left-x)
    for M in range(19):rec([],M)
    atomic(f'out/rounddegree14_D{D}_validation.json',dict(D=D,family=FAMILY,status='ZERO_VIOLATIONS',results=records,exhaustive_five_part_bound_checks=checks))

def add_cuts(lp,info,types,D,candidates,check=True,q=5):
    if check:
        v=json.load(open(f'out/rounddegree14_D{D}_validation.json'));assert (v['D'],v['family'],v['status'])==(D,FAMILY,'ZERO_VIOLATIONS')
        assert all(r.get('generated_model_checks',0)>0 for r in v['results'])
        prices=json.load(open(f'out/rounddegree14_D{D}_prices.json'));assert (prices['D'],prices['family'])==(D,FAMILY)
        assert prices['subsets']==json.loads(key(subsets(D)))
    assert all('Q14:'+','.join(map(str,m)) in lp.col for m in candidates)
    count=0
    for subset in subsets(D):
        bytotal=defaultdict(list)
        for m in candidates:bytotal[sum(m[j-1] for j in subset)].append(lp.col['Q14:'+','.join(map(str,m))])
        marginal={M:lp.var('RD14:'+','.join(map(str,subset))+':'+str(M)) for M in sorted(bytotal)}
        lp.add([(col,1) for col in marginal.values()],'=',1,f'XIV rounded degree normalization {key(subset)}')
        for M,col in marginal.items():lp.add([(col,1)]+[(j,-1) for j in bytotal[M]],'=',0,f'XIV rounded degree marginal {key(subset)} {M}')
        totals=[sum(a[2][j-1] for j in subset) for a in types]
        for t in range(2,26):
            row=[(lp.col[f'N{i}'],v//t) for i,v in enumerate(totals)]
            lp.add(row+[(col,-upper(M,t,q)*(q*q+q+1)) for M,col in marginal.items()],'<',0,f'XIV rounded degree upper {key(subset)} {t}')
            lp.add([(col,-v) for col,v in row]+[(col,lower(M,t,q)*(q*q+q+1)) for M,col in marginal.items()],'<',0,f'XIV rounded degree lower {key(subset)} {t}');count+=2
    print('ROUNDDEGREE14_CUTS',D,count,flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--D',type=int,required=True);ap.add_argument('--fresh',action='store_true');a=ap.parse_args();validate(a.D);table(a.D,a.fresh)
