# Independent check of the paired coherent-drift regeneration gap

2026-10-04. The proposed paired-word upgrade is valid, with the finite-N and clock qualifications stated in PAIRED-K3-COHERENT-DRIFT-REGENERATION-GAP.md.

The source's explicit Hessian formula and its established actual VALUE/displacement Lp limits imply directly

    A0,00−A1,00 -> 2beta cos(x0),
    T0,00,T1,00 -> lambda+beta cos(S0).

The radial Hessian terms vanish at the ancestor coordinates because their first components stay tight and their radial norms diverge. They vanish at the terminal coordinates because their O(N) first components divided by their O(sqrt(d_N)) radii are O(a_N). Periodicity removes the integer 2pi N terminal phase, and the t0-to-t1 displacement is o_Lp(1). No differentiated approximation is used. Every Hessian entry stays bounded by .6.

For either true or split bank assignment, decompose the tested entry exactly as

    T0(B)[A0(A)−A1(A)]+[T0(B)−T1(B)]A1(A).

The second term vanishes in L1 and the first converges by boundedness and the marginal L1 limits. Independence between approximation errors and the other factors is unnecessary.

The true Gaussian coordinate correlation is Cov(S_A,0,x_A,0)=c. For two complete banks independent conditional on the same root with clock t, it is Cov(S_T,0,x_A,0)=tc. The identity E[cos J cos K]=exp(−1)cosh(Cov(J,K)) for centered unit Gaussian J,K therefore proves the exact limiting difference

    2beta^2 exp(−1)[cosh(c)−cosh(tc)]>0,  t<1.

Conditional-bank independence and the tower property identify this with the difference AFTER the owned-root average. Positivity at every individual fixed root is neither required nor established.

All approximation bounds depend only on standard-marginal L1 errors, so convergence is uniform in t in [0,1]. This justifies fixed positive finite-clock averages and the ordinary unit-clock integral; integrating cosh(tc) gives sinh(c)/c. A nonvanishing amount of weight away from t=1 is needed for a uniform positive lower bound for varying finite clock rules.

For each fixed positive fine width, a finite N0 exists beyond which the gap exceeds half its positive limit. No certified numerical threshold is supplied. If the fine width vanishes with N, additive uniform convergence alone does not give a relative positive gap. Numerical native-family checks are diagnostics only.

The conclusion excludes this particular paired-factor regeneration identity. It proves no covariance/orientation lower bound solely from the selected-word discrepancy, and no impossibility for an additional compensating current or another complete producer.
