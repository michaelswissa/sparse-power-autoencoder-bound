"""Exact finite-support falsification of the proof. Python 3 standard library only.

All expectations, projections and comparisons use Fraction/integer arithmetic.
Irrational-root inequalities are compared by raising nonnegative sides to powers.
No sampling of input states, floating tolerance, network, model training or GPU.
"""
from fractions import Fraction as F
from itertools import product
from collections import Counter
import argparse, json, random, time, platform
from pathlib import Path

COUNTS = Counter()
BRANCHES = Counter()
STATES = 0

def check(name, lhs, rhs):
    COUNTS[name] += 1
    if lhs > rhs:
        raise AssertionError((name, str(lhs), str(rhs)))

def equal(name, lhs, rhs):
    COUNTS[name] += 1
    if lhs != rhs:
        raise AssertionError((name, str(lhs), str(rhs)))

def basis(rows):
    out, pivots = [], []
    for row in rows:
        v = list(map(F, row))
        for old, pivot in zip(out, pivots):
            z = v[pivot]
            v = [a-z*b for a,b in zip(v,old)]
        pivot = next((i for i,a in enumerate(v) if a), None)
        if pivot is not None:
            z = v[pivot]
            out.append([a/z for a in v]); pivots.append(pivot)
    return out

def inverse(a):
    n = len(a)
    b = [[F(v) for v in row]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for i in range(n):
        k = next(k for k in range(i,n) if b[k][i])
        b[i],b[k] = b[k],b[i]
        z = b[i][i]; b[i] = [v/z for v in b[i]]
        for k in range(n):
            if k!=i:
                z=b[k][i]; b[k]=[v-z*w for v,w in zip(b[k],b[i])]
    return [row[n:] for row in b]

def projection(w, n):
    b = basis(w); r = len(b)
    if not r: return [[F(0)]*n for _ in range(n)],r
    g = [[sum(x*y for x,y in zip(u,v)) for v in b] for u in b]
    gi = inverse(g)
    p = [[sum(b[a][i]*gi[a][c]*b[c][j] for a in range(r) for c in range(r))
          for j in range(n)] for i in range(n)]
    return p,r

def inputs(n, support, p):
    dist = {0:1-p}
    for x,w in support: dist[x] = dist.get(x,F(0))+p*w
    dist = [(x,w) for x,w in dist.items() if w]
    for state in product(dist, repeat=n):
        x = tuple(a for a,b in state)
        weight = F(1)
        for a,b in state: weight *= b
        yield x,weight

def moments(a, support, p, m):
    global STATES
    n=len(a); c=F(0); di=[F(0)]*n; zhi=[F(0)]*n; zlo=[F(0)]*n
    mass=F(0); centered=[F(0)]*n
    k=m*2**(m-2)
    for x,weight in inputs(n,support,p):
        STATES += 1; mass += weight
        for i,row in enumerate(a):
            y=sum(v*t for v,t in zip(row,x)); z=y-row[i]*x[i]; aa=row[i]*x[i]
            check('mean_value_pointwise',abs(y**m-z**m),k*(abs(aa)**m+abs(aa)*abs(z)**(m-1)))
            c += weight*x[i]*y**m
            di[i] += weight*y**(2*m)
            zhi[i] += weight*z**(2*m)
            zlo[i] += weight*z**(m-1)
            centered[i] += weight*x[i]*z**m
    equal('probability_mass',mass,1)
    for v in centered: equal('centering',v,0)
    return c,di,zhi,zlo

def run_case(w,n,support,p,m,scale):
    # Integral Gram entries permit efficient exact state enumeration.
    a=[[sum(row[i]*row[j] for row in w) for j in range(n)] for i in range(n)]
    proj,r=projection(w,n); ell=[proj[i][i] for i in range(n)]
    d=[a[i][i] for i in range(n)]; q=[sum(v*v for v in row) for row in a]
    equal('projection_trace',sum(ell),r)
    for i in range(n):
        check('leverage_nonnegative',0,ell[i]); check('leverage_at_most_one',ell[i],1)
        check('diagonal_leverage',d[i]**2,ell[i]*q[i])
        for j in range(n):
            equal('projection_range',sum(proj[i][k]*a[k][j] for k in range(n)),a[i][j])
    sm=sum(v**m for v in d); t=sum(v**(2*m) for row in a for v in row)
    qq=sum(v**m for v in q); h=sum(d[i]**(2*m-2)*q[i] for i in range(n))
    check('rank_cauchy',sm**2,r*h)
    check('holder_geometry',h**m,t**(m-1)*qq)
    check('combined_geometry',sm**(2*m),r**m*t**(m-1)*qq)
    mu2=sum(weight*x**2 for x,weight in support)
    muhi=sum(weight*x**(2*m) for x,weight in support)
    mumid=sum(weight*x**(m+1) for x,weight in support)
    c,di,zhi,zlo=moments(a,support,p,m); dd=sum(di)
    first=p*mumid*sm; remainder=p*mu2*sum(d[i]*zlo[i] for i in range(n))
    k=m*2**(m-2)
    check('correlation_decomposition',abs(c),k*(first+remainder))
    for i in range(n):
        check('moment_expansion_row',p*muhi*sum(v**(2*m) for v in a[i]),di[i])
        check('variance_jensen_row',(p*mu2*q[i])**m,di[i])
        check('conditional_jensen',zhi[i],di[i])
        check('moment_monotonicity',zlo[i]**(2*m),di[i]**(m-1))
        check('diagonal_from_variance',d[i]**(2*m)*(p*mu2)**m,ell[i]**m*di[i])
    check('moment_expansion_global',p*muhi*t,dd)
    check('variance_jensen_global',p**m*mu2**m*qq,dd)
    check('interpolation',p**(2*m-1)*muhi**(m-1)*mu2**m*t**(m-1)*qq,dd**m)
    check('first_term_squared',first**(2*m)*muhi**(m-1)*mu2**m,r**m*dd**m*mumid**(2*m)*p)
    check('remainder_squared',remainder**2,p*mu2*r*dd)
    improvement=2*c-dd
    if dd:
        check('complete_square',improvement*dd,c*c)
        # C^2 <= 2 K^2 r D (a^2 p^(1/m) + p mu2).
        delta=c*c-2*k*k*r*dd*p*mu2
        if delta>0:
            BRANCHES['combined_correlation_positive_delta']+=1
            check('combined_correlation_bound',delta**m*muhi**(m-1)*mu2**m,
                  (2*k*k*r*dd*mumid*mumid)**m*p)
        else:
            BRANCHES['combined_correlation_nonpositive_delta']+=1
            COUNTS['combined_correlation_bound']+=1
        delta=improvement-2*k*k*r*p*mu2
        if delta>0:
            BRANCHES['final_positive_delta']+=1
            check('final_bound',delta**m*muhi**(m-1)*mu2**m,(2*k*k*r*mumid*mumid)**m*p)
        else:
            BRANCHES['final_nonpositive_delta']+=1
            COUNTS['final_bound']+=1
        check('tensor_upper_optimized',c*c,p*mu2*r**m*dd)
        check('baseline_upper_optimized',c*c,p*mu2*n*dd)
        if c>0:
            # Bracket the true optimal scalar by exact bisection; fresh evaluation
            # ensures final-loss checks include positive near-optimal improvements.
            lo,hi=F(0),F(1)
            while hi**m<c/dd: hi*=2
            for _ in range(12):
                mid=(lo+hi)/2
                if mid**m<c/dd: lo=mid
                else: hi=mid
            cs,ds,_,_=moments([[lo*v for v in row] for row in a],support,p,m)
            equal('optimized_scale_C',cs,lo**m*c)
            equal('optimized_scale_D',sum(ds),lo**(2*m)*dd)
            ii=2*cs-sum(ds)
            COUNTS['optimized_strictly_positive']+=1
            assert ii>0
            delta=ii-2*k*k*r*p*mu2
            if delta>0:
                BRANCHES['optimized_final_positive_delta']+=1
                check('optimized_final_bound',delta**m*muhi**(m-1)*mu2**m,
                      (2*k*k*r*mumid*mumid)**m*p)
            else:
                BRANCHES['optimized_final_nonpositive_delta']+=1
                COUNTS['optimized_final_bound']+=1
    else:
        equal('zero_D',improvement,0)
    # Fresh enumeration at positive rational scale, not a formula-only check.
    if scale:
        scaled=[[F(v)*scale for v in row] for row in a]
        cs,ds,_,_=moments(scaled,support,p,m)
        equal('homogeneity_C',cs,c*scale**m)
        equal('homogeneity_D',sum(ds),dd*scale**(2*m))
    # Linear formula independently enumerated.
    if m==3:
        linear=sum(weight*(2*sum(x[i]*sum(a[i][j]*x[j] for j in range(n)) for i in range(n))
                   -sum(sum(a[i][j]*x[j] for j in range(n))**2 for i in range(n)))
                   for x,weight in inputs(n,support,p))
        equal('linear_trace',linear,p*mu2*(2*sum(d)-sum(q)))
        check('linear_rank_bound',linear,p*mu2*r)

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',default='results.json')
    args=parser.parse_args(); start=time.perf_counter(); rng=random.Random(20261002)
    supports=[[(1,F(1,2)),(-1,F(1,2))],
              [(1,F(1,4)),(-1,F(1,4)),(3,F(1,4)),(-3,F(1,4))],
              [(0,F(1,2)),(2,F(1,4)),(-2,F(1,4))]]
    matrices=[(1,[[0]]),(1,[[1]]),(2,[[1,0]]),(2,[[1,-1]]),
              (3,[[1,1,-1],[2,2,-2]]),(3,[[1,0,0],[0,1,0],[0,0,1]]),
              (4,[[1,2,-3,0],[0,-1,2,0]])]
    for n in (2,3,4):
        for r in range(1,n+1):
            matrices.append((n,[[rng.randint(-2,2) for _ in range(n)] for _ in range(r)]))
    cases=0
    for ix,(n,w) in enumerate(matrices):
        for s,support in enumerate(supports):
            for p in (F(0),F(1,100),F(1,3),F(1)):
                for m in (3,5,7):
                    scale=F(1,3) if ix in (1,3,6) and s==0 and p==F(1,3) else None
                    run_case(w,n,support,p,m,scale); cases+=1
    # Mutation controls: plausible but false strengthenings must be caught.
    c,di,_,_=moments([[1,1],[1,1]],supports[0],F(1),3)
    assert c>2  # C is not bounded by p mu_(m+1) S_m alone.
    assert sum(di)>4  # D is not equal to the pure-term moment lower bound.
    # A non-centered law invalidates conditional Jensen: Y=1+(-1)=0.
    assert abs(-1)**6>abs(1-1)**6
    result={'status':'PASS','seed':20261002,'cases':cases,'moment_enumeration_state_visits':STATES,
            'assertions':sum(COUNTS.values()),'by_check':dict(sorted(COUNTS.items())),
            'mutation_controls_detected':3,'seconds':round(time.perf_counter()-start,3),
            'final_comparison_branches':dict(BRANCHES),
            'count_notes':'State visits include repeated moment enumerations and exclude separate linear enumerations. Checks count easy implication branches; three controls are direct counterexamples, not harness mutations.',
            'python':platform.python_version(),'arithmetic':'exact rational, no tolerance',
            'matrices':len(matrices),'supports':3,'powers':[3,5,7],'p':['0','1/100','1/3','1'],
            'scope':'finite small-case falsification; not proof, not optimization'}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
