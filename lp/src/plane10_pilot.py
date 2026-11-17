"""Bounded coordinate CP-SAT gauge. No timeout is an exclusion."""
import json,time
from collections import Counter
from ortools.sat.python import cp_model
from geom import AG
from structure import planes_of
from plane10 import geometries
from rmax_model import build

def pilot(short,limit=300):
    cfg,mvec=((6,17),(0,3,27,12,4)) if short=='617' else ((5,15),(0,4,24,15,3))
    ag=AG(3,5);a2=AG(2,5);g=next(g for g in geometries() if tuple(g['signature'][:2])==cfg and g['signature'][3]==cfg[1])
    forced=[ag.line_of[frozenset(ag.pidx[(*a2.points[p],0)] for p in a2.lines[l][1])] for l in g['lines']]
    t0=time.monotonic();model,x,y=build(ag,forced,cfg[0],46,D=5,ident=False,proj=False,mvec=mvec)
    # Redundant plane degree identities and conditional Lemma C on every incidence.
    dp=[model.NewIntVar(0,5,f'pd{p}') for p in range(125)]
    for p in range(125):model.Add(dp[p]==sum(y[l] for l,(_,ps) in enumerate(ag.lines) if p in ps))
    incidence=0
    for _,_,ps in planes_of(ag):
        ls=[l for l,(_,pts) in enumerate(ag.lines) if pts<=ps]
        rr=sum(y[l] for l in ls);kk=sum(x[p] for p in ps)
        model.Add(sum(dp[p] for p in ps)==25+5*rr)
        for l in ls:
            model.Add(sum(dp[p] for p in ag.lines[l][1])>=25+5*rr-5*(kk-5)).OnlyEnforceIf(y[l]);incidence+=1
    proto=model.Proto();size=dict(variables=len(proto.variables),constraints=len(proto.constraints),lemma_C_rows=incidence,build_s=time.monotonic()-t0)
    print(short,'MODEL',size,flush=True)
    solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=limit;solver.parameters.num_workers=1;solver.parameters.random_seed=10;solver.parameters.log_search_progress=True
    st=solver.Solve(model)
    result=dict(leaf=short,config=cfg,m=mvec,status=solver.StatusName(st),wall_s=solver.WallTime(),branches=solver.NumBranches(),conflicts=solver.NumConflicts(),**size)
    if st in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        result['points']=[ag.points[p] for p in range(125) if solver.Value(x[p])];result['witness_lines']=[l for l in range(775) if solver.Value(y[l])]
    json.dump(result,open(f'out/plane10_pilot_{short}.json','w'),indent=2);print('PILOT_RESULT',result,flush=True)
if __name__=='__main__':
    pilot('617');pilot('515');print('PILOT_DONE',flush=True)
