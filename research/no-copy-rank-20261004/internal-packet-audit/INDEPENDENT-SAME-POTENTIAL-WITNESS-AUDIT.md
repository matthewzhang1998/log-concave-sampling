# Independent audit of the two same-potential matrix witnesses

2026-10-04. Independent algebra check of the selected K0 expectation counterexamples in INTERNAL-K3-ANCESTOR-TERMINAL-COMPARISON-CYCLE.md. This checks the covariance and heat calculations, not an all-rank producer or full paired K0−K1 noncancellation theorem.

## Findings

For the shared potential

    V(y)=y*T y/2−(delta/k^2)cos(k y1),
    T=[[.5,.1],[.1,.5]], delta=.1, k=2,

its gradient is the displayed g and its Hessian is T+delta cos(k y1)N. Since the eigenvalues of T are .4,.6, the claimed global [.3I,.7I] sandwich is valid. The orthogonal X,Y transformation of independent VS,VU is exact, and direct substitution gives the stated native terminal phase. Both factors come from this same potential.

At the zero coarse root, with q=.63125,beta=.00625 and L=exp[−(.8^2+.00625^2)/2],

    Cov(Ct,C0)
      >=L{exp[−(q^2+1)/2](cosh q−1)
                          −beta(1+exp(−1/2))}
      =.0670260701223.

The original minus independent-terminal-bank K0 expectation is EXACTLY delta^2 Cov(Ct,C0)N. Its (1,1) gap is therefore at least .0006702607012. The main checker obtains approximately .0007483063179 by one-dimensional Gaussian quadrature. The proof does not depend on that quadrature.

For independent physical-argument heat of width v=.5 at each frozen query label, alpha=exp(−1/2), and the (1,2) difference is EXACTLY

    .01(1−alpha) E Ct.

The bound E Ct>=L(exp(−q^2/2)−beta)=.5904236146 yields at least .0023231359014. This calculation requires additive heat at the frozen x(W),t0(W) physical arguments; it does not describe perturbing W and recomputing the nonlinear t0.

## Independent check of the owned-root integral

Let R,V,V' be independent full records, W=hR+vV, W'=hR+vV', h^2=.75,v^2=.25. For the linear terminal phase Q and ancestor phase P=kx1, the numbers

    A=Var Q=4.1540625,
    B=Var P=4,
    D=Cov(Q,P)=2.525

are correct. The two Gaussian identities give

    E cos Q cos P =exp[−(A+B)/2]cosh D,
    E_R[E(cos Q|R)E(cos P|R)]
                 =exp[−(A+B)/2]cosh(h^2D).

Hence their difference is .0489757113032. Replacing cos Q by cos(Q+beta sin P) changes each of the two defining terms by at most beta. This proves

    E_R Cov(Ct,C0|R)>=.0364757113032,

and the SAME-root split-bank K0 expectation gap is at least .0003647571130 N. Averaging the owned root later does not repair this particular independent-factor regeneration.

A slightly sharper error estimate uses |E(C0|R)|<=exp(−v^2B/2)=exp(−1/2) and replaces 2beta by beta(1+exp(−1/2)). It gives at least .0003893489468 N. The main note deliberately uses the simpler conservative bound.

## Scope qualification retained

These are K0 branch tests. The full six-query K3 embedding is supplied explicitly in the main note; Z,r,epsilon and q1,q2,q4,q6 remain actual nodes even though they do not enter this one branch's covariance formula. No cancellation statement for K0−K1 follows from the branch calculation alone.

The results refute two specific target-preservation claims: independently regenerating the terminal half of a same-record word, and independently heating two physical primitive arguments after freezing their original labels. They do not rule out regenerating the WHOLE coupled word, a new joint packet comparison, or another complete E-mean producer.


## Supplement: minimal-label inverse and literal primitive seed

For M=I, positive a,epsilon,s, and globally strongly convex C1 gradient maps, the main note's terminal-label inverse is exact:

    x0=g0^(-1)(p0), x1=g0^(-1)(p1), U=(x0−cS)/s,
    Delta=(x1−x0)/a,
    Z=[gt^(-1)(gt(S)−Delta)−S]/epsilon.

Global invertibility follows from strong convexity/coercivity, not a new oracle assumption. The forward and inverse substitutions cancel exactly. The tuple (S,p0) already recovers S,U; (S,p0,p1) recovers the complete native record, and hence the complete fine record at fixed coarse root and positive clock width. These identities concern the actual terminal construction labels, not arbitrary fine observers or an approximate reconstruction bound. The diagnostic checks both directions for the explicit smooth same-potential fixture; the maximum forward-recovery error is below 4e-13.

LOW30's generic order-two pair seed leaves a choice of finite tail realization. The main note now declares its specific choice J=D−s0 I_K, with I_K using degree-K Lagrange derivative weights at zero. Since E I_K approximates (E Df)p, E J approximates exactly m_tail=m−s0(E Df)p. Expanding the two N2 calls gives the displayed finite signal. Counting two D differences, K+1 interpolation nodes and one shared filter baseline gives K+6 original affine primitive VALUES before exact aliases. This is one legal concrete instantiation, not a universal count for every LOW30 tail implementation. The four fills and the order-rho^2 comparison from unit covariance to the selected stationary covariance remain stated explicitly.
