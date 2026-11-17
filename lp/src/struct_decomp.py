# structured point n=1, v=(1/4,0): preimage of grid pt (X,Y) is (X-Y-1/4, Y); F = x^2 - y
import sys
for q in [101,103,109,107]:
    sq={a*a%q for a in range(q)}; chi=lambda t:0 if t%q==0 else (1 if t%q in sq else -1)
    S=sorted(sq); w=pow(4,-1,q)
    tot=zero=axis=0; ax_in=0
    for X in S:
        for Y in S:
            x=(X-Y-w)%q; F=(x*x-Y)%q; inn = F in sq
            tot+=inn; zero+=(F==0)
            if X==0 or Y==0: axis+=1; ax_in+=inn
    G=len(S)**2
    # also: Hasse-style - is F(X,Y) expressible symmetric?
    print(q,q%4,'G',G,'in',tot,'in-G/2',tot-G/2,'zeros',zero,'axis pts',axis,'axis in',ax_in)
