# The rotation constraint defect is an actual projected-gradient offspring

New finite source construction, 2026-10-04. This resolves the source type of the NEW lambda defect in the variance-preserving auxiliary path3953238199cdc91934e689015a17ba01c24e9dd27871aae6fbf5d155b111a2c6. It does not by itself prove the whole nonlinear ambient endpoint comparison.

## 1. The exact record and defect

Let M be any known contraction, including singular M. Put

    A0=a sigma M, ||A0||<=1/2,
    D=(I-A0 A0*)^(1/2), D2=(I-A0* A0)^(1/2).

Use the REGRESSED independent standard blocks(S',U,Z,V'). The old physical Gaussian is

    S=D S'-A0 V'.

The original actual K3 output E'=E(S',U,Z) keeps its full finite query/feedback graph and ignores V'. Define

    H_Z(s)=gt(s)-gt(s+epsilon Z),
    L_lambda=ra M*[H_Z(S')-H_Z(S)].                    (1)

This is exactly the ambient lambda defect rM*(d'-d), not an averaged coefficient or a substituted reference. It has the shared S',Z records required by the old E' endpoint.

## 2. Four original affine queries on the COMPLETE master tape

Add one unused independent d-dimensional Gaussian Zpad. Set

    W=(S',U,Z,V',Zpad),
    C1=[I,0,0,0,0], C2=[I,0,epsilon I,0,0],
    C3=[D,0,0,-A0,0], C4=[D,0,epsilon I,-A0,0],
    signs=(+1,-1,-1,+1).

The U and padding columns are exactly unread by the new source. U is nevertheless retained on the complete common tape for E'. Both C1,C3 have unit row covariance and C2,C4 have row covariance(1+epsilon^2)I. The latter two are the actual twin queries, not rescaled standard queries.

Choose beta=1/2 and define the preliminary physical readout

    B0=beta[M*,0,0,-M* A0(I+D2)^(-1)]                 (2)

on the first four blocks. Since

    A0(I+D2)^(-1)A0*=I-D,

one has the EXACT identity

    B0 C_i*=beta M*, for i=1,2,3,4.                  (3)

In (3) the unused padding column of C_i is simply omitted when B0 is used. There is no inverse of M, A0, a, sigma or epsilon.

For ||A0||<=1/2,

    B0 B0* <= beta^2(1+||A0||^2)I <=5I/16.

Let F=(I-B0B0*)^(1/2) be a KNOWN fill and put B=[B0,F]. Then BB*=I and(3) remains true on the complete master tape. The fill has a fixed numerical gap at least11/16. Known dense matrix arithmetic is charged when M is not a cheap scalar-block map.

## 3. Literal genuine-gradient source and exact readout

Define the 5d-output source

    G_lambda(W)=(ra/beta) sum_i signs_i C_i* gt(C_i W). (4)

It uses four original gradient VALUES, their exact recorded query arrays and known row actions. It is the genuine full gradient of

    (ra/beta) sum_i signs_i Vt(C_i W).

Potential VALUES and original Hessian actions are not needed to execute(4). The U/padding outputs are zero. The physical identity is exact at every actual record:

    B G_lambda(W)=ra M*[gt(S')-gt(S'+epsilon Z)
                            -gt(S)+gt(S+epsilon Z)]
                  =L_lambda(W).                       (5)

Thus the new offspring is a projected-gradient source on the FULL standard common tape. It is not merely a gradient in V' conditional on already averaging away the old endpoint's roots.

Its complete first is bounded by

    ||DG_lambda||op <=(ra/beta) sum_i ||C_i||^2
                       <=12ra,                         (6)

for epsilon<=1. Both physical row/column bounds follow by the actual known B. No derivative convergence argument or inverse-width normalizer enters. With all geometry fixed in theta, the direct caller is at most C ra L_(gt,theta), using the four original caller paths. At the canonical scales this is O(A^1.5) if the original direct caller is O(A^-1/2).

The source is zero at the full fresh origin. Indeed sum_i signs_i C_i=0, so a common original-gradient anchor cancels exactly even before imposing gt(0)=0. This cancellation must use the same recorded origin and finite row version in all four calls; it creates no moving-anchor debt.

## 4. Whole-source energy and the companion blocks

Writing H=H_Z, the three possibly nonzero blocks of the signed sum in(4) are

    S' block: H(S')-D H(S),
    Z block: epsilon[gt(S+epsilon Z)-gt(S'+epsilon Z)],
    V' block: A0*H(S).                                 (7)

These are real source outputs. In particular, when Z=0 the physical defect L_lambda vanishes, but the Z output in(7) need not vanish. The whole gradient must not be assigned a pointwise E mark by looking only at(5).

Under the original unit first bound,

    |H(S')-H(S)|<=2|S'-S|, |H(S)|<=epsilon |Z|.

The exact Gaussian pair satisfies Cov(S'-S)=2(I-D), and ||I-D||<=C a^2 sigma^2. Therefore(7) gives the conservative ACTUAL full-source profile

    ||G_lambda||_p <= C_p sqrt(d) r a^2 sigma           (8)

for fixed p, epsilon<=1 and sigma in the stated guard. The physical L_lambda has the same upper envelope. The full first remains O(ra), not the stronger energy scale in(8). This is a new explicit offspring profile, not automatically a factor times the actual centered energy e of E'.

For the projection-M subclass, a sharper comparison to the actual original mark is possible using the separately proved rank-aware source bounds; it is not assumed for arbitrary singular M here. A consumer may instead retain this four-leaf gradient body at its actual outward mass O(ra), which is smaller than the original O(r) force mass. That retention statement still requires the literal parent/mean-composition paths when the body is inserted into another program; it is not inferred from the covariance bound alone.

## 5. Work, numerical maps and observer boundary

A standalone call to(4) uses four original VALUES. A joint evaluation of the complete E' and L_lambda uses eight: the six original K3 calls at S', plus gt(S) and gt(S+epsilon Z). Existing same-query feedback calls may be cached; a changed argument rebuilds them. FIRST/adjoint work uses the corresponding actual original HVP sites and known matrix sweeps. No differentiated HVP output is executed.

D,D2,(I+D2)^(-1),F are known matrix functions on fixed spectral gaps. Their finite polynomial/linear-algebra approximations have explicit operator-error budgets after all actual source readouts, source profiles and Gaussian dimension factors. A finite approximation is not silently substituted into the exact identities(3),(5); its known-matrix and query-law errors are paid coherently. The inverse(I+D2)^(-1) has a numerical norm bound, so this implementation creates no inverse-sigma sampling factor or small-singular-value cost for M. Dense known-matrix actions may still be nontrivial and remain charged.

The complete source relationship is now available to a same-record projected-gradient/mixed consumer: E' and L_lambda are functions of this SAME W. The conditional covariance frames from39532381 remain valid. If a proof conditions on(S',U,Z), however, the full-tape gradient theorem for(4) cannot average those observed roots without a matching joint endpoint comparison. Their shared dependence has not disappeared.

This construction closes the lambda offspring's source-type gate. It does not establish that the original ambient tuple is stationary, nor that conditional auxiliary normalization removes the remaining original-root orientation. Those are still distinct consumer questions.

## 6. Diagnostics

The adjacent checker passes3,500 noncommuting/singular-matrix cases covering
the known coisometry/fill, exact raw lambda readout, full-gradient Hessian,
all companion blocks and unread U/padding columns. These diagnostics do not
assert the still-separate nonlinear endpoint-current join.
