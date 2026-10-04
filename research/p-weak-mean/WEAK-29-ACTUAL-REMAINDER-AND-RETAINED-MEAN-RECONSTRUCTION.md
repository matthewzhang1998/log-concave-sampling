# Weak-mark 2.9 proxy, actual remainder, and the complete two-arm mean

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

**Status: reconstructed mathematical certificate, new version, 2026-10-04.**
This document is a new reconstruction, not a byte-identical copy of the former
weak-family proof. It records an executable local constructor relative to
explicit old-source and finite-compiler interfaces. The local mean, variance,
empirical comparison, and cost algebra below are proved here. The retained
host return has additional graph premises stated in Section 8. The verification
list in Section 11 distinguishes checked source text from premises still to be
reinstantiated for a particular call.

The historical weak-family pin was
`f4b9a9fdf27adf5f57eda2d376151428c8dc09c634ab0ce99b140d0d699eb724`;
its retained-source supplement was
`003c3bb972ca2bb4b7818ff1bdf04ecff50fac6756d156dcadf4762fef887ad9`.
Those identifiers are provenance only. They do not identify this document.

## 1. Source foundation and conventions

The source foundation used here is:

1. LOW30, SHA-256
   `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`.
   Relevant labels are `t30:lem:finite-columns`, `t30:eq:twin`,
   `t30:eq:proxy-readout`, `t30:prop:proxy`, `b27:compiler:mean`,
   `b27:compiler:paths`, `b27:compiler:serial`, and
   `b27:compiler:restoration`.
2. The explicitly completed OLD-R7-on-OLD-R13 hidden provider, SHA-256
   `38cf61f69f80371e0e0607446189213841a404fa9a89082b95fabdd1ca74a521`.
   It uses OLD R13 force/history occurrences, the completed OLD R7 local
   statistic, and the completed OLD R7 known-center decoder in that order.
3. The historical retained hidden-force certificate, SHA-256
   `38c3a9c1c5a0ce6211aaa5b7461465dc06ed8103e88bbe579fb1739064625c6e`,
   and its same-version batch certificate, SHA-256
   `cff71ce50596264c93b81b88dc8daf9426bbb438f1773e211e00b741f3912ce9`.
   The numerical 2.8/3.8 batch in that second certificate is not reused as a
   2.9 batch. Its dependency order, source identities, and restoration rules
   are the relevant foundation; the proxy is rebuilt at the mark below.

Fix the physical caller theta before sampling any new private record. Fix the
potential, geometry, row maps, source order, moment list, precision schedule,
and all integer/compiler choices before differentiating theta. Dimension of
the physical output is d; Dim is the common dimension factor in the inherited
moment bounds. All constants denoted Lambda may depend on a fixed compiler
order and fixed generation and polynomially on the declared public logarithms.
They are not asserted uniform in an increasing requested order.

Let 0<A<1 be the local heat/radius envelope. The canonical physical force has
amplitude r=A. We also state the amplitude-r formulas where helpful. All
source moments are conditional on the actual entering caller, and subsequently
integrated against its proved moment profile. A small Lp error is not a
pointwise bound at every possible caller.

## 2. The literal source batch: actual, retained, proxy, remainder

Start with a COMPLETE finite old force program

    F_actual(theta; W0) in R^d.

It includes every old history, observation, decoder, finite iteration, and
moving zero/anchor value it reads. Its retained STATE is constructed on the
same record before applying the final original-gradient readout. Denote the
resulting force by F_ret. For the completed CW7 input, the retained certificate
provides, at physical amplitude r,

    ||F_actual-F_ret||_Lp <= Lambda_p sqrt(Dim) r A^4.8,
    Lip F_actual + Lip F_ret <= Lambda r,
    caller(F_actual), caller(F_ret) <= Lambda r/sqrt(A).

The actual old source may have a weaker full recorded curl than F_ret. No curl
of F_actual is needed in the empirical branch below.

First center the ACTUAL finite retained graph node by node. If a displayed
acyclic portion is p_i(W)=g_i(b_i+C_i W+sum_(j<i) M_ij p_j(W)), compute its
same-version zero trajectory p_i^0 and q_i^0=b_i+sum M_ij p_j^0. Then

    z_i(W)=g_i(q_i^0+C_i W+sum M_ij z_j(W))-g_i(q_i^0)

is exactly p_i(W)-p_i^0. Each original-gradient difference is again a gradient
and is exactly zero at its zero argument. Store all q_i^0 and caller paths.
Old admitted whole blocks stay in their admitted coordinates. Finite iteration
versions stay finite; an exact fixed point is not substituted into F_ret.
Use zero-preserving retained substitutions so F_ret(0)=F_actual(0), or retain
the separately bounded actual-to-retained origin offset.

Apply the protected-row/signed-twin construction to this exact centered graph at

    J=29/10, epsilon_tw=A^(9/10),
    v_prot=19/15, A_prot=A^(34/15).

Here v_prot is the protected-row exponent in LOW30, not a variance allocation.
The buffer and twin errors both have relative force mark A^J. Indeed
1+3 v_prot/2=29/10 and epsilon_tw A^2=A^(29/10).
The protected projection has gamma >= A^(19/15)/Lambda. Set

    beta=(gamma^(-1)+epsilon_tw^(-2))^(-1/2),
    B=beta [P_rec Pi/gamma, -I/epsilon_tw].

Then BB*=I_d and beta is A^(9/10) up to the actual logarithmic factors.
The original recorded row P_rec and the proxy row B are different maps.

Unroll the symmetric signed system by the specified finite simultaneous
Picard iteration, starting from zero and using the centered primitive differences above. At zero
input each actual iterate is zero. Let C_tw
be the full affine row of that system and p_K its complete plus/minus record.
The raw finite gradient-lift readout is

    G_raw,K=(r/beta) C_tw* p_K,
    H_raw,K=B G_raw,K = r p_terminal,+,K.

The exact fixed-point reference G_raw,* is a genuine full gradient. G_raw,K
is the actual finite VALUE map and need not be a gradient. Every original
primitive query and finite iteration in G_raw,K is executed.

For the common-origin convention, store the actual finite zero trajectory
p_K(theta;0), and define

    G_fin(W) = (r/beta) C_tw* [p_K(W)-p_K(0)]
                  + B* F_actual(theta;0),
    H_fin(W)=B G_fin(W).

Use the identically translated exact reference G_* with p_* in place of p_K.
The additive vector is a gradient of a linear function of W, so G_* remains a
full gradient. It is a captured caller label, not a zero derivative in theta.
This gives the pointwise identity H_fin(theta;0)=F_actual(theta;0).

**Origin qualification.** Centering a proxy is not justified merely by a
Gaussian Lp approximation. Its error satisfies exactly

    E_fin(W)=[F_actual(W)-H_raw,K(W)]
                 -[F_actual(0)-H_raw,K(0)].

Thus the deterministic zero-trajectory error must also be bounded at the
actual caller. The direct centered-graph construction above removes this debt when the
actual-to-retained substitutions preserve the same origin: redo the native
row-perturbation/resolvent comparison on that graph, where only Gaussian rows
and centered twin perturbations change. The origin identity alone does not
prove its approximation order. The new reconstructed finite-origin certificate
with SHA-256 2bc878f9219188b73f4e05c57f84bbf81452f92d09c11d598f21f0746c36f93a
prints this exact procedure and its separate source premises. If only an uncentered proxy estimate is
available, add the explicit origin debt

    e_origin,p = ||F_actual(theta;0)-H_raw,K(theta;0)||_Lp(theta)

and do not drop it. In the intended centered native construction, removal of
a Gaussian row and the twin perturbation have zero VALUE price at the all-zero
fresh tape; finite zero restoration is still separately charged. This must be
verified on the chosen source transcript, rather than inferred from its law.

On the master record W=(W0,Z_tw), with the original source ignoring Z_tw, form

    E_fin = F_actual-H_fin
          = (F_actual-F_ret) + (F_ret-H_fin).

This is a pointwise same-record identity. The first term is present. No
independent proxy draw, retained-force substitution, or deletion of the twin
bank is performed inside E_fin. After parameters are chosen, freeze these
versions throughout every occurrence below.

## 3. Finite first, profile, and approximation bounds

The finite terminal recurrence and the identity
C_tw B*=beta e_terminal,+ give, independently of derivative convergence,

    ||DG_fin|| <= Lambda r/beta,
    ||B DG_fin|| + ||DG_fin B*|| <= Lambda r,
    ||D_theta G_fin|| <= Lambda r/(beta sqrt(A)),
    ||B D_theta G_fin|| <= Lambda r/sqrt(A).

The translated zero path is included in the two caller rows. The finite
terminal formula proves the recorded-P_rec curl of H_fin; no B*E_fin curl is
assumed or needed here.

At r=A the full lifted radius is Lambda A^.1, the physical first is Lambda A,
and the full/physical caller bounds are Lambda A^(-.4), Lambda sqrt(A).
The actual remainder has

    Lip E_fin <= Lambda A,
    caller(E_fin) <= Lambda sqrt(A),
    e_E,p := ||E_fin||_Lp
       <= Lambda_p sqrt(Dim)(A^3.9+A^5.8)
             + e_num,p + e_origin,p.

Absorb the A^5.8 term when A is small. The displayed e_E,p is UNcentered and
includes the complete actual old tail. The numerical/origin rows are assigned
before the nominal A^3.9 notation is used. An energy estimate for a difference
of marginal laws would not prove this statement.

The complete finite calls satisfy their original actual query-domain and
caller-moment contracts. Finite proxy substitution is controlled in VALUE at
those calls. It does not permit differentiating G_fin-G_* as a small error.

## 4. The complete positive mean program

Choose fixed positive numerical variances v_H=v_E=1/2. Let m be a fixed
arbitrary-order genuine-gradient compiler order (called B in LOW30). With
independent complete branch records, execute

    F0=F_actual(theta;0),  G_c=G_fin-B*F0,
    Y_H = F0 + sqrt(v_H) B N_m(G_c/sqrt(v_H); Omega_H),
    Y_E = (1/N) sum_{i=1}^N E_fin(theta; Wi) + sqrt(v_E) Z_E,
    Y = Y_H+Y_E.

Here N_m is the literal VALUE compiler with each source call replaced by the
fixed finite G_c implementation, using the stored F0 and its full caller path.
For an already zero-anchored source G_c=G_fin. The explicit F0 readout avoids
assuming a small profile for a nonzero deterministic anchor. Wi are independent COMPLETE E_fin records
conditional only on theta. A Wi evaluates both F_actual and H_fin on the same
padded record. Z_E and Omega_H are independent of all Wi. All old aliases
inside one complete source record remain unchanged.

The N_m output lies in the full lifted output space; B is applied afterward.
Because BB*=I_d, the Gaussian reference of Y_H has covariance v_H I_d. No
normalizer is applied directly to the non-gradient E_fin. No queried unknown
mean, covariance, Hessian tensor, or covariance square root appears.

Analyze Y_H first with G_*-B*F0; restore every finite VALUE occurrence at its actual
query, and then return the target mean to E H_fin. If eps_G bounds the relevant
restorations and V_m is the complete absolute source-value path certificate,
the two steps cost at most a known constant times

    V_m eps_G + ||E(H_fin-H_*)|| <= (V_m+1) eps_G,

with actual fixed variance/readout constants retained. If query-dependent
path envelopes are random, use the full sum of their L2 products, not a product
of marginal means. The joint finite-restoration lemma in LOW30 supplies that
form. The target is therefore N(E H_fin,v_H I_d), with error

    eps_H <= Lambda_m sqrt(n_H) (Lambda A^.1)^m
                       + eps_H,rest + eps_H,filters + eps_H,rows.

All terms refer to the same finite source batch. In particular
E F_actual=E H_fin+E E_fin exactly.

## 5. The empirical arm: proof of law, energy, and first

At fixed theta let S=E_fin(W), mu=E S, ell=Lip S, and
T_N=N^(-1) sum_i(S_i-mu). Gaussian independence gives

    ||T_N||_2 <= ||S-mu||_2/sqrt(N),
    ||D T_N||_op <= ell/sqrt(N).

For the interpolation X_t=mu+t T_N+sqrt(v_E)Z_E, Gaussian Riesz duality on the
complete Wi banks gives

    d/dt E phi(X_t)
      = t E[(R T_N)(D T_N)* : D^2 phi(X_t)].

The coefficient is independent of Z_E. One integration against that unchanged
keep gives an L2 velocity bound

    t ||S-mu||_2 ell/(sqrt(v_E) N).

The Gaussian Riesz L2 contraction and conditional Jensen justify the bound.
Integrating t from zero to one proves

    W2(Law Y_E, N(mu,v_E I_d))
       <= ell ||S-mu||_2/(2 sqrt(v_E) N)
       <= ell e_E,2/(2 sqrt(v_E) N).                 (E-law)

This result uses no curl, high derivative, clipping, or independence of the
coordinates inside S. It requires independent complete copies after theta is
fixed and the explicit numerical keep v_E.

The actual, uncentered emitted signal satisfies for each required fixed p

    ||N^(-1) sum_i S_i||_p
       <= ||mu|| + C_p ||S-mu||_p/sqrt(N) <= C_p e_E,p.

Its fresh first is ell/sqrt(N); its theta first is at most the original caller
bound. The latter is a coherent average and receives no sqrt(N) improvement.
A bias common to every copy, including numerical source bias, also receives
no such improvement. The exact finite E_fin identity, rather than an unbiased
approximation to some different E, determines mu.

By coupling the two independent branches and adding them,

    W2(Law Y,N(E F_actual,I_d))
       <= eps_H + ell e_E,2/(2 sqrt(v_E)N).          (mean-law)

The entire output is a positive probability kernel with fixed covariance
shares. Neither variance is obtained by subtracting an unknown matrix.

## 6. Precision schedule and the real copy power

For a requested normalized absolute allowance Lambda sqrt(Dim) A^Q, first fix
the finite structural moment list and actual active dimension n_H. If
n_H<=Lambda Dim A^(-kappa_H), a sufficient strict fixed-order prior condition is

    m>=2,  m/10-kappa_H/2 > Q,
    equivalently m > 10Q+5 kappa_H together with m>=2,

with additional fixed/logarithmic slack for all constants and the genuine
compiler's m-dependent small-radius guard. When using the ordinary-profile argument below, assign its absolute mean
error at least order A even if Q<1: the conservative choice is the same
inequality with Q replaced by max(Q,1). This is an absolute error grade;
the actual LOW30 compiler starts at order m=2. A fixed prior floor
cannot be driven to arbitrary precision at fixed A just by more Picard steps.

Choose N from the actual assigned bound, for example

    N=max(1,ceil(ell_bound e_E,2,bound/(2 sqrt(v_E) eps_E))).

At nominal e_E=Lambda sqrt(Dim) A^3.9 this is a real inverse-heat copy cost
A^(-k) with any strict k>max(0,Q-4.9), allowing logarithmic allocations. It is
not free and is not claimed sublinear in requested order.

The geometric finite-source error q^K with q<=3/8 is made smaller than the
allocated VALUE tolerance AFTER all original source-domain envelopes, path
sums, inverse consuming widths, and scalar weights are enumerated. A sufficient
form is K>=log(C_query/path/eps_G)/log(8/3), rounded upward. Rebuild the entire
matched finite batch if this changes a proxy or old provider version. Include
both the random query and deterministic-zero restoration; never refine only
one of two copies of the same subtracted source.

The real-valued finite program and its numerical encoding are distinct:
known Gaussian rows/coisometries, finite filters and arithmetic have separately
allocated encoding errors. Enforce their covariance identities or pay their
explicit errors. Uniform smooth first bounds apply to the specified real
VALUE circuit; a discontinuous rounded implementation does not inherit those
bounds by assertion.

## 7. Actual physical first and ordinary profile

The full source radius of the genuine-gradient prior is A^.1, but every
source-dependent physical output path ends in B DG_fin or B D_theta G_fin.
The executed scalar-block affine maps commute with the physical readout at
that last source occurrence. LOW30 `b27:raw:split` and `b27:compiler:paths`
therefore give

    fresh residual first(Y_H) <= Lambda_m A,
    physical caller(Y_H) <= Lambda_m sqrt(A).

The finite G_fin firsts used here are those of Section 3. A Gaussian-law error
is not used as a derivative estimate. Fixed captured zero values have their
same original caller paths. The empirical branch has the bounds in Section 5.
Consequently Y has fresh nonlinear first Lambda A and physical caller
Lambda sqrt(A), apart from its explicit known identity-covariance carrier.

After subtracting the explicit F0 and known Gaussian carrier, the ordinary
physical nonlinear profile is Lambda_p sqrt(Dim) A, under the inherited
ordinary source profile ||F_actual-F0||_p<=Lambda_p sqrt(Dim)A. To see the
dimension-safe point: a d-output residual with operator first Lambda A has
Hilbert Jacobian at most Lambda sqrt(d)A, irrespective of its input tape size.
Gaussian Hilbert Poincare bounds its centered Lp profile; its actual mean is
controlled by E(H_fin-F0), that actual source profile, and the already assigned
mean error. Add the empirical e_E,p profile. Do not replace this proof by a
radial bound using sqrt(the complete polynomial-sized tape).

## 8. Full retained graph and conditional re-entry

The auxiliary kept program retains the WHOLE finite Y_H graph, its exact
known Gaussian carrier, every original VALUE occurrence, all moving anchors
and actual zero computations. It keeps sqrt(v_E) Z_E and deletes only the
empirical signal N^(-1) sum E_fin(Wi). The deletion is measured at the actual
chronological caller and costs C_p e_E,p, not its smaller law error. In the
canonical common-origin implementation it also vanishes at the full fresh
origin. If that convention has not been verified, retain an explicit offset
and certify the remaining signal separately.

The graph return needs more than this deletion inequality. For its native
RAW/hidden host use the following literal admissions:

(a) Expand the scalar-affine N_m transcript, but retain previously admitted
symmetric implicit blocks in their admitted coordinates. At the output-nearest
proxy occurrence simplify B G_fin to the physical terminal VALUE. Its outgoing
original-gradient mass is O(A). Every deeper full-G edge has positive small
radius A^.1 (or the compiler's fixed positive powers). Apply the finite
controlled-block balancing to this fixed enumerated graph, retaining all
original rows and offsets. Numerical Picard iterates are not relabeled as
independent gauge levels.

(b) At a phase insertion restore its actual physical sqrt(A), time, harmonic
and source-readout factors. Compare state outputs directly at actual callers,
then transport through the explicit kept small-A map. The decoder's inherited
one-half coefficient has a bounded geometric resolvent. Replace centers by
the finite kept-center construction; do not differentiate the deleted full
polynomial-tape source during this transport.

(c) Under those host rows, the deletion costs a physical STATE error
Lambda_p sqrt(Dim) sqrt(a) a^3.9. Earlier retained-force replacements have
an additional physical predecessor and four fixed kept-center substitutions
are higher order. Thus the new retained mark is K=3.9. This is a state proof;
one does not invert a force estimate to obtain it.

(d) After attachment, the final force is still an original-gradient terminal
at that kept state, with actual readout r/sqrt(a). Its known leading Gaussian
row is sqrt(a)P_new; its retained state remainder has first sqrt(a)a. The
terminal square lift therefore has curl O(ra) on the COMPLETE owned record.
Old P_rec, B, and P_new are not identified.

(e) The independent empirical keep is unread by its source evaluations. In a
hidden decoder the final completed known-center child supplies a fixed part
of the total carrier; earlier statistic/observation/decoder records ignore
its final protected increment. For proxy mark J the necessary protected heat
is a^((2J+1)/3). Weak J=2.9 requires cutoff exponent strictly above 34/15;
all fixed J<3.9 are covered by a cutoff strictly above 44/15, with the original
first-crossing/lower-ratio and logarithmic margins retained.

Items (a)-(e) are the retained-source verification to perform for each fixed
completed host. The historical CW7 certificate supplies the OLD input source;
it is not by itself a completed new high-order known/hidden host. A local
arbitrary-Q mean does not alone establish a new hidden law or sampler rate.

## 9. Active tape, full tape, and tiny amplitudes

After the empirical signal is deleted from the auxiliary kept graph, every Wi
record is unread, provided all retained centers were replaced literally and
all retained offsets have certified read sets. Its independent Z_E remains.
The next gradient source may use that literal active selector. If E_keep
selects its active coordinates, its padding back to the actual complete record
is exactly

    G_pad(W)=E_keep* G(E_keep W),   B_pad=B E_keep.

This is an identity on the kept graph; it does not compress the actual sampled
F_actual or E_fin program. The actual E still evaluates F_actual and its
padded H_fin on the same full W. Every empirical call remains in the executed
work and full tape. Nonlinear H compiler copies remain active.

If the fixed-generation retained graph has polynomial-logarithmically many
blocks, kappa_H=0 is legal for its next gradient prior. If that count has an
inverse-A factor, retain the resulting kappa_H in Section 6 instead.

For a known small signed scalar c, |c|<=1, a literal amplitude-scaled occurrence
may use the canonical source directly and emit

    Y_H,c=c Y_H + sqrt(v_H(1-c^2)) Z_extra,
    E_c=c E_fin.

This has the desired mean and variance without division by c. All numerical
mean errors scale |c| automatically. The original compiler schedule stays at
the canonical A envelope. Structural c=0 is a direct Gaussian. This wrapper
avoids assuming a nonlinear compiler is homogeneous merely because its input
was multiplied by c. Arbitrary matrix readouts require the separately declared
componentwise row/variance adapter.

## 10. Full VALUE, FIRST, arithmetic, and clipping cost

Let Q_F be the complete original-gradient VALUE count for F_actual, Q_G the
count for one G_fin evaluation including its finite system and required zero
path, and Q_H the corresponding B G_fin physical count. A safe complete
remainder bill is Q_E<=Q_F+Q_H, with all internal old calls expanded at their
actual heat labels. The two-arm bill is

    Q_total <= Q_m(L) Q_G + N(Q_F+Q_H) + Q_known,

plus any separately stored zero evaluation not already included. It is
additive. There is no N^2 term, and no claim that the empirical copies evaluate
only retained F_ret. Conversely, an empirical E copy does not execute N_m.

LOW30 supplies a finite symbolic enumeration with Delta=1/16:

    Q_2(L)<=c_2(1+L)^k_2,
    Q_(b+Delta)(L)<=c_b(1+L)^k_b[1+Q_b(c_b(1+L))].

Thus for each fixed m, Q_m(L)<=C_m(1+L)^K_m. Constants and exponents come from
the full queue enumeration: all trees, signs, filters, positive clocks,
replicas, nested source calls and readouts. No uniform growing-m claim is
made. Without actually enumerating a chosen m, the symbolic recurrence is the
available quantitative bill; it is not a printed numerical call count.

One actual forward or adjoint sweep in one entering direction costs a constant
times the expanded original call count in original HVPs. Reverse accumulation
adds all parent contributions before reusing a stored child. Discarded primal
records require paid replay. Neither a completed normalizer derivative nor an
HVP derivative is a new oracle. Gaussian sampling and known linear-map actions
also cost work: record their actual dimensions/actions in Q_known; a dense
coisometry is not assumed to have a free application. Full storage includes
all N complete E records unless streaming is chosen and its reverse replay is
charged. In a serial RAW host use every old call at its actual heat; do not
replace its cost by a same-heat shortcut.

**Clipping accounting.** The empirical arm uses all original VALUES directly
and has no clipping operation. The inherited N_m program in the cited compiler
is an executed finite affine/VALUE transcript with finite filters and positive
clock rules. The operator-ball projection in LOW30's polynomial-comparison
proof is expressly an analytical hybrid, not an executed gate; no derivative
or evaluation of that projection is charged as a primitive. Its tail and
filter calibration enter eps_H. LOW30's separate terminal covariance chapter
has a clipped-parent TV argument; that terminal constructor is not invoked in
this two-arm mean. Any clipping or finite machine truncation introduced by a
particular implementation must instead be listed and charged, with its VALUE,
first and tail consequences. This certificate licenses no unlisted clipping.

## 11. Verification status of this reconstruction

Directly checked against historical/current source text:
- the completed old provider and its retained/actual source distinction;
- protected/twin normalization, beta exponent and finite physical columns;
- exact-gradient reference versus finite VALUE chronology;
- original mean compiler, complete path and serial recurrence statements;
- analytical-only operator gates versus executed empirical/mean transcript;
- the explicit local two-arm mean/variance, buffered empirical proof, first,
  energy and additive-count algebra in Sections 4-7 and 10.

Required to instantiate before promoting this to an unconditional new host:
- the exact chosen finite source/zero/anchor list and its deterministic-origin
  debt, caller profiles, precision floors and uniform finite firsts;
- actual enumeration of the selected compiler order and B-dependent guards;
- the complete kept-center/controlled-block/protected-terminal verification
  of Section 8 on the selected new known/hidden program;
- active and full dimension records and all inherited heat-labelled costs.

These are explicit source contracts, not conclusions drawn from a marginal
Gaussian approximation. This new reconstruction preserves the former local
constructor's quantitative content without claiming that an old audit hashes
or an old global theorem automatically applies to its rewritten bytes.
