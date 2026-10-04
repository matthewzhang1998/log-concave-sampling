# Scalar native hidden-column true Gram from a mixed VALUE square

New bounded source result for the d=1 native twin c8b20bf4. This replaces the
preliminary scalar skew-pair completion by an actual one-energy VALUE source.
It relies on the proved one-incoming Hilbert cubic-padding estimate in LOW30
`t30:lem:actions`; it does not assert a higher-order padding estimate.

## 1. Exact scalar orientation and the conditional roots

Use the orthogonal master coordinates(S,U,Z) from the native-twin port. The
physical square lift is f=(E,0,0), with recorded row selecting S. At a whole
OU clock c^2+v^2=1 write

    J_f(X)=[j_S(X),j_U(X),j_Z(X);0,0,0;0,0,0].

All j_l are scalar. Hence the physical true orientation is EXACTLY

    O_E=integral w(t) E_X[j_U(X)^2+j_Z(X)^2] dt.         (1)

The S forward square already equals its true square. No feedback term is
omitted from(1). The U-column bound is C ra<=kappa; the Z-column bound is
C ra^2 epsilon<=kappa. The native source's complete first remains A.

Fix a clock and coarse X. For l=U or Z, expose all OTHER fine coordinates in
a coefficient record R. Define the literal scalar chart

    h_(l,R)(z)=[E(other fine coordinates,c X_l+v z)
                           -E(other fine coordinates,c X_l)]/v. (2)

The order of coordinates in this display is the original(S,U,Z) order.
Every E is the COMPLETE same finite native callback, with its feedback and
twin aliases unchanged. Both branches of E are kept together.
The chart is scalar, has active first at most kappa, complete old-root first
C A/v, and conditional centered energy whose L2 over R,X is at most C e/v.
Its derivative mean B_(l,R)=E_z h'_(l,R) satisfies

    E_R B_(l,R)=j_l(X).

For the two factors below draw TWO independent complete conditional records
R1,R2 at the SAME coarse X. Using one identical R twice would target
E_R B_R^2, which is generally the wrong coefficient.

## 2. Literal mixed square with no small readout division

Let I1(p,z1) be the standard K-node common-Gaussian first-chaos response of
h_(l,R1), using its SAME z1 at every sign/node. It is not divided by kappa.
Let beta=1/2 and define

    R2(w,z2)=[h_(l,R2)(beta w+sqrt(1-beta^2)z2)
             -h_(l,R2)(-beta w+sqrt(1-beta^2)z2)]/(2 beta).

Use independent z0,z1,z2 and emit

    S_(l,j)(p)= [R2(z0+I1(p,z1),z2)
                         -R2(z0-I1(p,z1),z2)]/2.        (3)

All four outer calls share z0,z2,I1,R2 and the same finite source version.
Their anchors cancel exactly. The source is pointwise odd in p. Both E
branches inside every h call remain paired. No derivative matrix is queried.

For analysis only normalize each h by kappa and I1 by kappa. Then(3) is the
LOW30 square formula with source radius kappa and incoming padding sigma=kappa,
i.e. effective mu=kappa^2. Its executable form(3) contains no inverse-kappa
readout. The one-incoming Hilbert proof uses one inner source energy and the
outer source's uniform first bound; it therefore applies to these two
independent conditional source copies as well as to identical copies.
It gives, after the private roots are averaged in the stated order,

    E_private S_(l,j)(p)=E_R2 B_R2 E_R1 B_R1 p + error,
    ||error||_(L2_X,p)<=C[kappa^3 e/v
                                      +kappa e 4^(-K)/v],
    ||S_(l,j)||_Lp<=Lambda_p kappa e_p/v.               (4)

The conditional product in the first line is j_l(X)^2. There is no transpose
ambiguity in dimension one. The first-filter error is paid at its actual
outer kappa factor. The cubic estimate is an averaged one-Hilbert statement;
one must not replace it by a pointwise high-moment Taylor bound on I1.

## 3. Actual complete first and all-owned-root aggregation

The source's public p first is at most Lambda kappa^2. Its own inner-source
paths gain the outer kappa. Its outer private/root paths are bounded directly
by the two VALUE terms, without differentiating their small difference.
Consequently its COMPLETE private first, including coarse X and both R banks,
is safely Lambda A/v. Its physical caller is the actual original caller
bound divided by v (and fixed/log factors), not a derivative of the mean.

Define one source on a fresh scalar public G by

    S(G)=sum_j w_j [S_(U,j)(G)+S_(Z,j)(G)].             (5)

Every coarse X, both conditional R banks, all filters and innovations in(5)
belong INSIDE one complete private set, independent of the old consumer and G.
The actual positive sums sum w/v<=Lambda then give

    energy(S)<=Lambda kappa e,
    full private first(S)<=Lambda A,
    public first(S)<=Lambda kappa^2,
    physical caller(S)<=Lambda sqrt(A)

when the original caller is canonical. The reference mean is the CONSTANT
O_clock G, with calibration error Lambda kappa^3 e plus assigned finite
filter/source/clock errors. There is no remaining coefficient-root mixture:
all roots are averaged in this conditional comparison at G fixed.

The omitted orientation clock can be assigned its own delta kappa e budget.
For an A kappa e target choose delta<=c A with its actual fixed/log allocations.
All clock floors/finite first-response query widths and source tolerances are
listed before execution. The weighted first sum, rather than a worst-node
1/v_min bound, is the actual complete return.

## 4. Positive negative-covariance consumer at a new independent public

Choose fixed b0,zeta0>0 and emit

    b0 G - q S(G)/(2 b0) + zeta0 Z_keep,
    b0^2+zeta0^2=v_new,

with a fixed untouched keep and all records independent of the old reserve.
At G fixed, conditional-private Gaussian Riesz and the keep price the actual
source fluctuation by its complete private first times its energy:

    Lambda A kappa e.

Calibrate the conditional mean to O_clock G afterward. The remaining positive
coherent square has norm at most C kappa^3 e: op O_E<=C kappa^2 and its HS
bound is C kappa e, so O_E^2 has HS bound C kappa^3 e. Numerical coefficient
errors and their actual readouts are separately included. Thus this source
can target N(0,v_new-q O_E) with error

    Lambda[A kappa e+kappa^3 e+delta kappa e]+e_num.

The mean/calibration and the whole private comparison preserve the original
same-source callback. No law retaining the actual source-zero carrier is
assumed. This is a one-energy correction with actual signal energy kappa e.

The original forward reserve still has its own error (A^(7/3)e at the simple
LOW30 choice mu=A^(4/3)). That row is not improved merely by the new orientation
source. A stronger forward provider must be separately supplied if the total
negative-Cov kernel is claimed at the finer A kappa e scale.

## 5. Literal work, primitive graph, and scope

With K first-filter nodes, each(3) uses at most(2K+4) COMPLETE E evaluations;
its h anchors cancel in the shared differences. The safe two-column bill is
sum_j(4K_j+8) Q_E(K_source), plus known Gaussian/row arithmetic. For the
unrolled two-node twin, Q_E(K_source)<=4K_source+2 is a safe original-gradient
VALUE bound before the declared identical-query cache savings. Every source
version and all actual anchor/caller paths remain charged. No inverse-A
replica bank or inverse-kappa sampling multiplier is introduced.

FIRST/adjoint sweeps follow the actual finite VALUE circuit, using original
HVPs at its recorded queries and paid replay when necessary. Source precision
is amplified by the actual1/v and first-filter factors. Matrix/derivative
oracles are absent.

The outer readout in(5) is w_j, not sqrt(w_j). Original output-nearest source
coefficients therefore sum through w/v, giving O(A) mass. Inner and outer
charts at a given column use the SAME v; the receiving query cancels its own
divisor. For the U column the cross-source primitive edge is O(r); for the Z
column it additionally carries the original epsilon or a epsilon path. The
native feedback edges a and all moving anchors remain. A full retained host
still needs its explicit carrier, original Gaussian rows and physical terminal
proof; alternatively the new source signal has the stated kappa e deletion
budget, which is genuinely smaller than an A e covariance signal when
kappa<A.

This result is scalar. In general dimension the conditional hidden column is
a nonsymmetric matrix, and(3) targets j_l^2 rather than j_l j_l*. The same-query
adjoint/true-Gram completion is the unresolved general-dimensional port.
