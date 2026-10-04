# The natural paired-potential lift has an unbounded third-jet row

This NEW bounded screen concerns the exact K3 source c8b20bf4. It excludes one specific natural full-gradient lift under the admitted bounded-Hessian contract. It does not exclude other lifts or an averaged orientation estimate.

## 1. Exact paired potential and its companions

In independent query coordinates(S,x,Z), write Delta=gt(S)-gt(S+epsilon Z), x1=x+aM*Delta, t0=S+aM g0(x), t1=S+aM g0(x1). Consider the analytical paired potential

    Psi=r[Vt(t0)-Vt(t1)].

For D=Ht(S)-Ht(S+epsilon Z), H01=H0(x1), Hplus=Ht(S+epsilon Z), direct chain rule gives

    grad_S Psi=E-r a² D M H01 M*gt(t1),
    grad_x Psi=r a[H0(x)M*gt(t0)-H01 M*gt(t1)],
    grad_Z Psi=r a² epsilon Hplus M H01 M*gt(t1).

These signs and orders are correct. They already require original Hessian actions as output VALUES, which the raw oracle model does not provide. Even if they are treated analytically, the full gradient need not have a uniform first-radius bound: differentiating the x companion exposes original third jets multiplying an unmarked terminal VALUE.

## 2. Exact SAME-potential scalar native family

Fix m=1/2,d=1/10 and use the SAME primitive at both nodes,

    g(x)=m x+d sin(kx)/k,
    V(x)=m x²/2-d cos(kx)/k².

Then .4<=g'(x)<=.6 globally, g(0)=0, and V is strongly convex. Let the scalar coupling M be1. Choose an integer L>=8 and set

    a=r=1/(2L), epsilon=a^.9,
    k=2pi N, N an arbitrarily large integer,
    S=1, x=0, Z=-L/(N epsilon).

This is an actual finite source record. In the original orthogonal coordinates it corresponds to the fixed finite U=-c/s. The displayed Z tends to zero as N grows, rather than escaping to a large-radius query.

Since kS=2pi N and k(S+epsilon Z)=2pi(N-L), the two sine values vanish exactly. Therefore

    Delta=L/(2N), x1=a Delta=1/(4N)=pi/(2k),
    t0=1, t1=1+a C/k,   C=pi/4+.1.

In particular

    g''(x)=0,  g''(x1)=-d k,
    g'(x)=.6, g'(x1)=.5,
    g'(t0)=.6, g'(t1)=.5+.1cos(aC),
    g(t1)=.5+[.5aC+.1sin(aC)]/k >=.5.

The physical copied-feedback derivative D is zero at this record, but that does not remove the x-companion defect.

## 3. The first radius of the full gradient diverges

At fixed S,Z, x1 differs from x by a constant. Hence

    Psi_xx=ra[g''(x)g(t0)+a g'(x)^2g'(t0)
                 -g''(x1)g(t1)-a g'(x1)^2g'(t1)]
           =ra[.1k g(t1)
                +a{(.6)^3-(.5)^2(.5+.1cos(aC))}].

The bracketed second term is positive. Consequently

    Psi_xx >= .05 r a k.                              (1)

For fixed admissible r,a,epsilon, k can grow arbitrarily while EVERY original gradient first bound and all exact native source contracts remain fixed. Thus no uniform Lipschitz bound for grad Psi follows from those contracts. Paired use of the SAME potential does not cancel the dangerous third-jet row in this candidate.

The original E itself has complete first O(r), recorded curl O(ra) and pointwise |E|<=r a² epsilon |Z|, uniformly in k. More strongly, scalar monotonicity gives its actual centered Gaussian energy

    .4³ r a² epsilon <= e <= .6³ r a² epsilon.

The lower bound follows from D_ZE>=.4³ r a² epsilon and conditional Gaussian integration by parts; the upper bound follows from the pointwise Lipschitz chain. Thus the gradient-lift failure is not an artifact of replacing the actual E mark by a foreign large-energy source.

## 4. Scope

Equation(1) excludes the natural full paired-potential gradient as an admissible bounded-first VALUE source under the present oracle/regularity contract. It does not rule out a different finite VALUE lift, a source-qualified weak coefficient/current estimate, or a row bound on the Gaussian output-swap defect of the actual E. In particular the independently proved scalar direct-Gram program never executes grad Psi and is unaffected.
