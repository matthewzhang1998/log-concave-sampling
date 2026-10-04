# Paired K0−K1 regeneration fails on the actual coherent-drift native family

2026-10-04. This addendum upgrades the separate K0 screen to a rigorous PAIRED selected-word counterexample, using the held exact native coherent-drift family. It proves a nonzero averaged word gap after independent terminal-factor regeneration, shared-root integration, and ordinary positive clock averaging. It does not infer a covariance/orientation lower bound from that word gap alone.

## 1. Exact source family and coupling under test

Use the actual source from c-native-twin/coherent-drift/NATIVE-COHERENT-DRIFT-ORIENTATION-LOWER-BOUND.md. There is ONE potential in each dimension,

    R(y)=sqrt(1+|y|^2), Phi_d(y)=y0(R(y)−1),
    V_d(y)=lambda |y|^2/2+delta Phi_d(y)−beta cos(y0),
    lambda=1/2, delta=beta=1/40,
    g_d=grad V_d, H_d=Dg_d,
    .4 I <= H_d <= .6 I.

Set c=s=1/sqrt(2), M=e0e0*,

    a_N=r_N=N^(−5/9), epsilon_N=N^(−1/2)=a_N^(.9),
    d_N=nearest integer to [2pi N/(delta a_N)]^2.

For sufficiently large integers N, a_N<=1/16. On the actual complete W=(S,U,Z), use the six native original gradient queries and define

    x=cS+sU,
    h_N=a[g_d(S)0−g_d(S+epsilon Z)0],
    x1=x+h_N e0,
    t0=S+a g_d(x)0 e0,
    t1=S+a g_d(x1)0 e0,
    A0=H_d(x), A1=H_d(x1), T0=H_d(t0), T1=H_d(t1),
    K_i=T_i M A_i.

Because M has rank one,

    (K_i)00=(T_i)00 (A_i)00.                            (1)

These are analytical Hessian labels of the SAME finite original source, not executed Hessian VALUES.

For a covariance clock t in [0,1], draw independent complete standard triples R,V_A,V_T and put

    W_A=sqrt(t)R+sqrt(1−t)V_A,
    W_T=sqrt(t)R+sqrt(1−t)V_T.

Conditional on the SAME owned R, the two complete fine banks are independent. Each bank separately rebuilds the complete original K3 graph, preserving all its within-bank terminal/ancestor/feedback aliases. The tested wrong substitution breaks only the cross-factor alias:

    Q_N = E[(K0−K1)(W_A)]00,
    Qsplit_N(t)
       =E[(T0(W_T) M A0(W_A)−T1(W_T) M A1(W_A))]00.    (2)

Both expectations already include the owned-root average. The same W_T is used for its two terminal factors; the same W_A is used for its two ancestor factors. Thus the calculation does not manufacture a failure by independently unpairing i=0 and i=1.

## 2. Direct Hessian limits, with no differentiated approximation

The explicit first Hessian entry is

    H_d(y)00=lambda+delta[3y0/R(y)−y0^3/R(y)^3]
                         +beta cos(y0).                (3)

The held source proves, on standard complete marginal W, in every fixed Lp,

    h_N -> −pi,
    a g_d(x)0−2pi N ->0,
    d0_N:=g_d(x)0−g_d(x1)0
               ->lambda pi+2beta sin(x0).              (4)

In particular d0_N=O_Lp(1), and t1=t0−a d0_N e0. These are VALUE and displacement limits of the actual source. We do not differentiate any convergence assertion in (4).

Instead substitute directly into (3). At x and x1, the first coordinates are O_Lp(1), while their radii diverge like sqrt(d_N), so their radial Hessian term tends to zero in probability. It is uniformly bounded by an absolute constant. The cosine at x1 tends to cos(x0−pi)=−cos(x0). Thus, for every fixed finite p,

    (A0)00 -> lambda+beta cos(x0),
    (A1)00 -> lambda−beta cos(x0).                      (5)

At each terminal, the first coordinate is S0+2pi N+o_Lp(1), up to the additional o_Lp(1) displacement a d0_N for t1. The unshifted bulk S coordinates retain Gaussian radius of order sqrt(d_N), and N/sqrt(d_N)=O(a_N). Hence t_i,0/R(t_i)->0 in probability. Again the radial term in (3) vanishes by bounded domination, and cosine periodicity plus (4) gives

    (T0)00,(T1)00 -> lambda+beta cos(S0).                (6)

The original global Hessian bound .6 controls every factor in (5)–(6). Therefore these limits are Lp limits, and products converge without an unproved derivative-error estimate. The limits remain valid marginally on BOTH W_A and W_T for every clock coupling t, because both complete records have the original standard Gaussian law.

## 3. Exact limiting paired gap after root integration

Write x_A,0=c S_A,0+s U_A,0. Equations (1),(5),(6) imply

    (K0−K1)(W_A)00
       ->2beta[lambda+beta cos(S_A,0)]cos(x_A,0),

    split paired word00
       ->2beta[lambda+beta cos(S_T,0)]cos(x_A,0).        (7)

The Gaussian pairs (S_A,0,x_A,0) and (S_T,0,x_A,0) have unit marginal variances and covariances c and tc, respectively. Thus

    lim Q_N
      =2beta lambda exp(−1/2)+2beta^2 exp(−1) cosh(c),
    lim Qsplit_N(t)
      =2beta lambda exp(−1/2)+2beta^2 exp(−1) cosh(tc).

The paired mean-word discrepancy is therefore

    Delta_N(t):=Q_N−Qsplit_N(t) -> G(t),
    G(t)=2beta^2 exp(−1)[cosh(c)−cosh(tc)].              (8)

It is strictly positive for every t<1. In particular, for fine width v=1/2, t=3/4,

    G(3/4)=0.0000536367688883... >0.                    (9)

For every fixed t<1 there is a finite N0(t) such that each actual finite native source N>=N0(t) has Delta_N(t)>=G(t)/2>0. This is a rigorous finite-family existence statement. The present addendum does not provide a certified numerical N0; the diagnostic finite N choices do not substitute for one.

The common lambda contribution cancels between original and split targets, but the beta^2 shared Gaussian correlation term does not. This rules out the possibility that K0 and K1 automatically cancel the factor-regeneration error on the general actual native family.

At the physical B=ra E[K0−K1] normalization, the discrepancy in B00 is ra Delta_N(t), of the full ra word scale. With the source's actual C_anc=[cI,sI,0], the (S0,U0) selected-skew entry is s B00, so this linear selected-skew coefficient also changes. This last observation is only a linear coefficient identity; it is not a claim about a quadratic covariance or an orientation contraction.

## 4. Uniformity in the clock and positive clock averages

The convergence in (8) is uniform over t in [0,1]. To see this, bound each product error by a bounded factor times the marginal L1 error of one Hessian entry. All four marginal errors in (5)–(6) have the same standard-input law regardless of t. Their common bound tends to zero. No convergence of a conditional Hessian derivative or changing-clock field is needed.

Consequently:

- For every fixed t0<1, Delta_N(t)>=G(t0)/2>0 for all t<=t0 once N is sufficiently large.
- Any FIXED finite positive clock sum with at least one positive-weight node t<1 has limiting gap sum_j w_j G(t_j)>0.
- More generally, a family of finite positive clock sums with bounded total mass and a fixed positive amount of weight at t<=t0<1 has a uniform positive limiting lower bound. Mere presence of a node whose weight tends to zero is insufficient.
- For the ordinary unit covariance-clock integral,

    lim integral_0^1 Delta_N(t) dt
      =2beta^2 exp(−1)[cosh(c)−sinh(c)/c]
      =0.0000805426949365... >0.                       (10)

Here the integral convergence follows from the same uniform error bound (or bounded dominated convergence). Thus neither K0/K1 pairing, owned-root averaging, nor a fixed positive covariance-clock average removes this explicit word gap.

If all clock weights approach the zero-fine-width endpoint t=1, the limiting gap can shrink. Equation (8) alone is not a uniform relative lower bound for arbitrary t_N->1. A source-dependent adaptive clock rule must be checked under its actual mass distribution; it is not automatically covered by the fixed-clock statement.

## 5. What is and is not excluded

This exact same-potential native family excludes the claimed target identity obtained by regenerating the terminal complete bank independently of the old ancestor bank, even when i=0,1 pairing and the same coarse root are respected. Regenerating the WHOLE coupled K0−K1 on one common fresh record remains legal and preserves its distribution; it also preserves the internal joint query dependence that the original primitive-pair composition has not compiled.

Each tested program is finite at every N. A complete source occurrence uses the original six gradient VALUES and a first/adjoint sweep uses its original HVP sites. A conservative pair of complete banks costs twelve original VALUES before exact semantic aliases, with dimension d_N, Gaussian generation, known rank-one actions and all precision costs charged. The proof introduces no inverse-gradient oracle, Hessian-as-VALUE instruction, numerical observer, or new source potential. The dimension grows with N, so no fixed-dimensional quantitative obstruction or dimension-free query bill is claimed.

The primitive comparison-order audit and its sufficient small new interface remain unchanged: a new joint selected-skew service may retain only captured theta, the owned R and the incoming p while integrating all fine records internally; another complete E-mean producer may bypass that service. This addendum does not require arbitrary fine observers at a successful new service boundary.

Finally, a nonzero selected-word difference does NOT by itself establish a nonzero error in B B*, the full physical mixed-adjoint orientation, a particular whole law, or every future correction program. The held coherent-drift source proves its actual orientation by a separate energy/isometry argument. That separate theorem is not reproduced or inferred from (8) here.

## 6. Source pins and diagnostics

The source-qualified ingredients are:

- c-native-twin/coherent-drift/NATIVE-COHERENT-DRIFT-ORIENTATION-LOWER-BOUND.md, §§1–4: actual potential, dimensions, Hessian formula, and the three Lp limits used in (4).
- p-native-twin/THIRD-ITERATE-TWO-NODE-WEAK-TWIN-PORT.md: original finite six-query K3 source, protected/native rows and scales.
- p-native-twin/MATRIX-K3-SAME-QUERY-TWO-SIDED-WORD-CONTRACT.md: definition of K0−K1 and its same-query selected-skew role.
- INTERNAL-K3-ANCESTOR-TERMINAL-COMPARISON-CYCLE.md: exact regeneration operation being audited, original VALUE/first/caller cost rules, and the sufficient smaller new service boundary.

The adjacent checker uses the exact radial sufficient statistics from a six-by-six correlated Gaussian Wishart Gram rather than allocating the enormous ambient vectors. It checks the explicit Hessian identities and paired mean diagnostics, together with deterministic Gaussian formulas for the limiting gap and clock integral. These finite samples diagnose the actual family; equations (3)–(10), not Monte Carlo estimates, prove the result.
