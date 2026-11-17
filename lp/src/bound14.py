"""Redundant finite column bounds, derived using original integer rows only.

Each bound is rounded UP to a power of two. Numerical certificates are still
checked on the unmodified original model; these bounds never create proof rows.
"""
from plane10_cert import bump

def upper_bounds(lp):
    by,bc,B=bump(lp,{}, {})
    assert B>0 and min(bc)>0
    assert all(y>=0 for y,row in zip(by,lp.rows) if row[2]=='<')
    def upper(a,b):
        assert a>=0 and b>0
        if a==0:return 0
        q=(a+b-1)//b
        result=1<<((q-1).bit_length())
        assert result*b>=a
        return result
    bounds=[upper(B,a) for a in bc];sources=[None]*len(bounds)
    positive=0
    for index,(row,rhs,sense,name) in enumerate(lp.rows):
        if rhs>=0 and row and all(v>=0 for v in row.values()):
            positive+=1
            for j,v in row.items():
                if v>0:
                    previous=sources[j]
                    old_rhs,old_coef=(B,bc[j]) if previous is None else (lp.rows[previous][1],lp.rows[previous][0][j])
                    if rhs*old_coef<old_rhs*v:
                        bounds[j]=upper(rhs,v);sources[j]=index
    assert all(type(x)==int and x>=0 for x in bounds)
    return bounds,dict(method='positive original rows and exact positive Farkas bump; upper powers of two',bump_rhs=B,positive_rows=positive,zero_columns=bounds.count(0),maximum=max(bounds)),sources

def repair_bounds(lp,y,sources):
    from plane10_cert import combo
    from collections import defaultdict
    columns,_=combo(lp,y);groups=defaultdict(list)
    for j,i in enumerate(sources):
        if i is not None and columns[j]<0:groups[i].append(j)
    for i,js in sorted(groups.items(),key=lambda item:(lp.rows[item[0]][1],item[0])):
        row=lp.rows[i][0]
        multiplier=max(0,max((-columns[j]+row[j]-1)//row[j] for j in js))
        if multiplier:
            y[i]+=multiplier
            for j,v in row.items():columns[j]+=multiplier*v
    return y
