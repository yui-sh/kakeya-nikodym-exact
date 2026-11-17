"""Numerical row projection; certificates are lifted to the unchanged full LP.

Drop all rows involving X/Y/Z/G variables and then discard unused columns.
This is a relaxation. A feasible projection is explicitly NOT full-model
feasibility. Any Farkas vector is padded by zero on dropped rows and verified
against every column of the original model before normal certificate saving.
Arguments and output checkpoints are those of probe14.py.

Optional XIV_Q_SAMPLE restricts numerical degree-hull columns. That experiment
is NOT a relaxation. Every omitted column is checked exactly, with a full-model
repair attempted; a failure only saves additional canonical Q seeds and remains
open. Numerical feasibility is never promoted to a full-prefix primal.
"""
import runpy,os,ast,json,random
from types import MethodType
import coverage14 as c
from planecount5 import LP

BaseMaster=c.Master14

class ProjectedMaster(BaseMaster):
    def __init__(self,lp,method='phase1'):
        self.full_lp=lp
        self.core_projection=os.environ.get('XIV_CORE_PROJECTION')=='1'
        self.drop_prefixes=tuple(os.environ.get('XIV_DROP_PREFIXES','XYZG'))
        assert self.drop_prefixes and set(self.drop_prefixes)<=set('XYZG')
        excluded={j for j,n in enumerate(lp.names) if n.startswith(self.drop_prefixes)}
        if self.core_projection:
            coreg=set()
            for row,_,_,name in lp.rows:
                if name.startswith('F8 Y1 ') and (ast.literal_eval(name[len('F8 Y1 '):])[:2] in ((5,15),(6,17)) or os.environ.get('XIV_CORE_EXTRA')=='1' and ast.literal_eval(name[len('F8 Y1 '):])[0]==6):
                    coreg.update(j for j in row if lp.names[j].startswith('G'))
            corez=set()
            for row,_,_,name in lp.rows:
                if name.startswith('F8 Y2 ') and coreg.intersection(row):
                    assert {j for j in row if lp.names[j].startswith('G')}<=coreg
                    corez.update(j for j in row if lp.names[j].startswith('Z'))
            excluded={j for j,n in enumerate(lp.names) if n.startswith(('X','Y')) or n.startswith('G') and j not in coreg or n.startswith('Z') and j not in corez}
            print('THIN14_CORE_VARIABLES',len(coreg),len(corez),flush=True)
        projected=[];self.keep_rows=[];self.row_signs=[]
        for i,(row,b,s,name) in enumerate(lp.rows):
            removed=[v for j,v in row.items() if j in excluded]
            sign=1
            if removed:
                if not self.core_projection:continue
                if min(removed)>=0:s='<'
                elif s=='=' and max(removed)<=0:s='<';sign=-1
                else:continue
            kept={j:sign*v for j,v in row.items() if j not in excluded and v}
            target=sign*b
            is_bound=name.startswith(('XIV upper','XIV lower'))
            if not is_bound and s=='<' and target>=0 and all(v<=0 for v in kept.values()):continue
            if s=='=' and not kept and target==0:continue
            projected.append((kept,target,s,name))
            self.keep_rows.append(i);self.row_signs.append(sign)
        used=sorted({j for row,_,_,_ in projected for j in row})
        self.omitted_q=set()
        sample=int(os.environ.get('XIV_Q_SAMPLE','0'))
        if sample:
            assert sample>0
            qs=[j for j in used if lp.names[j].startswith('Q14:')]
            assert qs
            vectors={j:tuple(map(int,lp.names[j][4:].split(','))) for j in qs}
            selected=set(random.Random(140007).sample(qs,min(sample,len(qs))))
            for d in range(len(vectors[qs[0]])):
                selected.add(min(qs,key=lambda j:vectors[j][d]))
                selected.add(max(qs,key=lambda j:vectors[j][d]))
            seedfile=os.environ.get('XIV_Q_SEED_FILE')
            if seedfile:
                seeds=set(json.load(open(seedfile))['columns'])
                selected.update(j for j in qs if lp.names[j] in seeds)
            self.omitted_q=set(qs)-selected
            used=[j for j in used if j not in self.omitted_q]
            self.q_seeds=[lp.names[j] for j in sorted(selected)]
            print('THIN14_RESTRICTED_Q',len(selected),'of',len(qs),'full-column verification mandatory',flush=True)
        small=LP();col={j:small.var(lp.names[j]) for j in used}
        for row,b,s,name in projected:
            small.add([(col[j],v) for j,v in row.items() if j in col],s,b,name)
        # The unchanged final degree-bound rows are needed by Master.solve.
        # Retain all final bounds even when implied by nonnegativity, because
        # the shared numerical Master updates these final twelve rows.
        # The filtering above exempts these bound rows.
        assert self.keep_rows[-12:]==list(range(len(lp.rows)-12,len(lp.rows)))
        small.matrices=MethodType(c.fast_matrices,small)
        if hasattr(lp,'canonical14'):
            small.canonical14={name:lp.canonical14[name] for name in small.names if name in lp.canonical14}
            small.model14=lp.model14
        self.small=small
        print('THIN14_PROJECTED',self.drop_prefixes,len(lp.names),len(lp.rows),'to',len(small.names),len(small.rows),flush=True)
        super().__init__(small,method)

    def solve(self,lp,seconds):
        status,data,detail=super().solve(self.small,seconds)
        detail['numerical_projection']='core G/Z with sign-valid row projection' if self.core_projection else 'drop rows involving '+','.join(self.drop_prefixes)+'; remove unused columns'
        if status=='LP-certified':
            y=[0]*len(lp.rows)
            for i,sign,v in zip(self.keep_rows,self.row_signs,data):y[i]=sign*v
            if self.omitted_q:
                from plane10_cert import combo,bump
                from lift11 import atomic
                cc,rhs=combo(lp,y)
                negative=[j for j,v in enumerate(cc) if v<0]
                assert set(negative)<=self.omitted_q
                if negative:
                    by,bc=bump(lp,{}, {})[:2]
                    repair=max((-cc[j]+bc[j]-1)//bc[j] for j in negative)
                    yy=[a+repair*b for a,b in zip(y,by)]
                    repaired,rrhs=combo(lp,yy)
                    if rrhs<0 and min(repaired)>=0:
                        y=yy;detail['full_Q_repair']=repair
                    else:
                        additions=sorted(negative,key=lambda j:cc[j])[:128]
                        path='out/qseed14'+c.CERT_SUFFIX+'.json'
                        atomic(path,dict(columns=sorted(set(self.q_seeds+[lp.names[j] for j in additions])),
                               status='unfinished restricted numerical model',negative_full_columns=len(negative),
                               restricted_rhs=rhs,full_repair_rhs=rrhs,
                               reason='Candidate fails omitted degree-hull columns; no coverage credit'))
                        print('THIN14_Q_PRICING_FAILED',len(negative),rhs,rrhs,path,flush=True)
                        return 'open',None,dict(detail,reason='restricted Q candidate fails full-column check',q_seeds=path,negative_Q_columns=len(negative))
                detail['restricted_Q_columns']=len(self.q_seeds)
            detail.update(c.verify(lp,y))
            print('THIN14_FULL_COLUMN_VERIFIED',detail,flush=True)
            return status,y,detail
        if status=='numerical-feasible':
            from lift11 import atomic
            selected=[dict(m=list(map(int,n[4:].split(','))),weight=float(data[j])) for j,n in enumerate(self.small.names) if n.startswith('Q14:') and data[j]>1e-9]
            selected.sort(key=lambda x:-x['weight'])
            record=dict(status='diagnostic numerical projection only; not an exact primal',degree_vectors=selected,detail=detail)
            atomic('out/projected14'+c.CERT_SUFFIX+'.json',record)
            detail['positive_degree_vectors']=selected[:32]
            return 'open',None,dict(detail,reason='projected/restricted numerical feasibility only; full prefix not solved')
        return status,data,detail

if __name__=='__main__':
    c.Master14=ProjectedMaster
    runpy.run_path('probe14.py',run_name='__main__')
