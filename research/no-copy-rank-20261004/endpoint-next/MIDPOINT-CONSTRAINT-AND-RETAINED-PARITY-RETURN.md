# The opposite word is an exact nonlinear input port, not an ordinary mixed covariance

New exact source identities and a restricted-generator obstruction, 2026-10-04. These statements concern the actual K3 genealogy. They do not rule out a new joint constrained-input channel. Throughout this note the covariance field stays a field on its retained roots.

## Result first

There is a finite six-VALUE representation of the entire physical E through an opposite-mode baseline and a nonlinear input. Its chain rule gives exactly the terminal-difference/common word and terminal-sum/opposite word. Shrinking the nonlinear input by changing the auxiliary width does not improve the latter word: the actual readout and incoming first cancel that apparent gain exactly.

A second obstruction is structural. The fresh opposite baseline is odd in its auxiliary Gaussian V, whereas E and the ancestor contrast ignore V. Consequently every retained ordinary mixed covariance field between either of those old sources and the baseline is odd in the retained coarse V root, and is identically zero when that coarse root is zero. The missing terminal-sum/ancestor-difference word is independent of V and nonzero at the terminal-flat native fixture. Thus simply applying the already-closed mixed-covariance service to these sources does not realize the needed field, even if a full-gradient baseline lift were granted for free.

The needed operation must couple the ancestor contrast into the baseline's selected V port. That coupling is exactly the nonlinear input constraint displayed below. Neither a constant covariance replacement nor ordinary conditional mixed fields have removed it.

## 1. Exact same-record midpoint program

Use the actual finite K3 record

    x0=cS+sU,
    Delta=gt(S)-gt(S+epsilon Z),
    x1=x0+a M*Delta,
    p0=g0(x0), p1=g0(x1),
    t0=S+aM p0, t1=S+aM p1,
    E=r[gt(t0)-gt(t1)].

Put q=p0-p1, pbar=(p0+p1)/2, tbar=S+aM pbar. Fix any h>0 and set sigma=h/a and kappa0=ra. Define the actual baseline VALUE source, on one fresh standard V,

    N_h(W,V)=kappa0/(2h)
          [gt(tbar(W)+hMV)-gt(tbar(W)-hMV)],
    mu_h(W)=a q(W)/(2h)=q(W)/(2sigma).

Here W=(S,U,Z). The full source call to N_h uses the original two feedback terminal VALUES, two original ancestor VALUES, and two displayed terminal VALUES, for six original gradient VALUES before exact-key caching. mu_h alone needs the first four. At V=mu_h the final terminal query arrays are literally t0,t1. Therefore

    E(W)=2sigma N_h(W,mu_h(W))                         (1)

is an exact finite-program identity with the original aliases and finite versions. It is not an identity only in expectation. Changing an argument rebuilds the complete original feedback source; only exact semantic matches may be cached.

N_h is odd in V and zero at V=0, for every fixed original W. Its fresh V energy and first obey

    ||N_h(W,V)||_p <= C_p kappa0 ||M V||_p,
    Lip_V N_h <= kappa0 L ||M||,

under the original terminal first bound L. These are actual unmarked body bounds. The root/caller bound is instead

    Lip_W N_h <= C kappa0/h,

because D_W tbar is bounded but Ht(tbar+hMV)-Ht(tbar-hMV) need not be small uniformly in the C2 class. Every original external caller contributes its literal query path, including the kappa0/h readout. The source needs only original VALUES; a first/adjoint sweep needs only original HVPs at those queries.

For M=I, N_h is a genuine gradient in V conditional on W. This fact alone does not make it a full-tape genuine gradient or admit the nonlinear incoming value mu_h.

## 2. The exact word split is the chain rule of the input constraint

Hold S,Z fixed and differentiate with respect to the common ancestor coordinate x=x0, so x1=x+aM*Delta. At the literal constrained evaluation V=mu_h, write Ti=Ht(ti), Ai=H0(xi). Then

    D_V N_h = (kappa0/2)(T0+T1)M,
    D_x mu_h = a(A0-A1)/(2h),
    D_x tbar = aM(A0+A1)/2,
    D_tbar N_h = kappa0(T0-T1)/(2h).

The two paths through (1) give

    2sigma (D_V N_h)(D_x mu_h)
       = (kappa0/2)(T0+T1)M(A0-A1),                   (2)

    2sigma (D_tbar N_h)(D_x tbar)
       = (kappa0/2)(T0-T1)M(A0+A1).                   (3)

Their sum is exactly kappa0(T0 M A0-T1 M A1). No matrix commutes and no transpose is inserted for free.

In (2), the 1/h in the smaller input mu_h, the h/a in the outer readout, and the selected baseline coefficient cancel exactly. Choosing h=A^delta, or any other positive finite width, does not introduce an extra heat power in the returned missing word. This is a statement about the full physical path, not merely the size of one intermediate VALUE.

## 3. The ancestor-input source is genuinely smaller, but its caller restores the scale

The actual source q has an important bounded structure. Put C=[cI,sI,0] and J_Delta=D_W Delta. Direct differentiation gives

    D_W q=(A0-A1)C-a A1 M*J_Delta.

The square lift C*q has a symmetric first term C*(A0-A1)C. Hence, for ||M||<=1 and epsilon<=1,

    Lip q <=C,
    ||Curl(C*q)||op <=C a,
    Lip mu_h <=C a/h,
    ||Curl(C*mu_h)||op <=C a^2/h.                     (4)

These statements concern the actual firsts, not derivatives of a small approximation. Under the original first bounds the literal VALUE envelope is

    |mu_h| <= C a^2 epsilon |Z|/h.                    (5)

For the nondegenerate M=I subclass with both original Hessians in a fixed positive numerical sandwich, the already-proved native auxiliary lower bound makes (5) a bound by the actual original centered energy:

    ||mu_h||_p <= C_p e/(r h),                        (6)

with the specified Gaussian moment factor. No inverse-singular-value assertion is made for general M. Outside that subclass keep the nominal source profile or a separately proved projected actual-energy mark.

The smaller source (4) is therefore a legitimate near-gradient VALUE object. This is useful partial structure, and it is not an impossibility result for a compiler that can genuinely consume it. But the actual port in (1) has incoming Lipschitz at most 2sigma kappa0 = 2r h. Multiplying (6) by this incoming factor returns only O(e), with no new positive heat power. At the terminal-flat fixture the corresponding selected word equality (2) is exact and nonzero.

Changing normalization cannot turn this into an energy-gain edge by itself. Any approximation to mu_h must be propagated through 2r h, and any approximation to N_h through 2sigma. Its numerical budget must therefore be selected after both readouts are fixed. Unknown or zero actual e cannot absorb an absolute numerical floor.

## 4. Why ordinary retained mixed fields do not select the opposite word

This obstruction applies before any mixture-to-expectation step. Let F(W) be either the original E, the ancestor input mu_h, or any original-root VALUE source ignoring the fresh V. Pad it by zero derivative columns in V. Use one literal whole-input Gaussian covariance node with positive width v and root coefficient c:

    W=cR_W+vX,    V=cR_V+vY,

where X,Y are independent standard fine records. Define the actual conditional Gaussian Jacobians J_F and J_N, with all original nonlinear queries in N rebuilt on the same complete W fine record.

Since N(W,-V)=-N(W,V), its W derivative is odd in V and its V derivative is even. F has no V derivative. Consequently

    J_F(R_W,R_V)=J_F(R_W),
    J_(N,W)(R_W,-R_V)=-J_(N,W)(R_W,R_V),

and the ordinary retained mixed field is

    S_FN(R_W,R_V)=Sym[J_F(R_W) J_(N,W)(R_W,R_V)*].

It satisfies the exact parity identities

    S_FN(R_W,-R_V)=-S_FN(R_W,R_V),
    S_FN(R_W,0)=0.                                    (7)

The same holds for every finite positive sum of such clocks, keeping all their actual coarse roots. If each clock has its own auxiliary root, set every such root to zero. It also holds for signed linear combinations with known coefficients and for known deterministic output adapters. None of these statements replaces the field by its expectation.

As an additional check, the complete ordinary covariance vanishes:

    Cov(F(W),N(W,V))=0,                               (8)

because E_V[N|W]=0. Equation (8) alone would not prove the retained-field claim; (7) is the stronger same-root statement.

The missing opposite word in (2) ignores the newly introduced V and is generally nonzero. Thus a sum of the ordinary mixed fields in (7) cannot equal that word on the same retained record. In particular, merely granting N a small full-gradient source and feeding (F,N) to the closed mixed reserve would target the wrong field.

The V derivative D_V N, which contains the terminal SUM, can contribute only if the partner has a V input column or if a new root-dependent observer/current supplies that port. The original F has neither. Setting the baseline incoming port to mu_h is exactly such a coupling, and returns the nonlinear constraint (1). Multiplying an odd field by a root observer and integrating by parts is another possible route, but it is a new matched same-endpoint current with its actual firsts, not an ordinary mixed-covariance substitution. The present argument does not rule that route out.

## 5. Explicit test on the terminal-flat native fixture

Use the previous exact fixture with M=I,

    T=[[1/2,1/20],[1/20,1/2]], delta=1/40,
    g(y)=T y+delta q0 psi(y1/q0)e1,
    q0=a epsilon,
    psi(z)=(z+1/2) chi(10(z+1/2)),

where chi is the smooth compact bump from the pinned terminal-flat note. At S=U=0,Z=e1,

    T0=T1=T, A0=T, A1=T+delta N, N=e1 e1*,
    p0-p1=q0 T^2 e1,
    tbar=-(a q0/2)T^2 e1,
    mu_h=(a q0/(2h))T^2 e1.

Both constrained terminal points lie in the exactly linear terminal region. Thus (1) evaluates to the original E=ra^2 epsilon T^3e1. The common term (3) is exactly zero, while (2) is

    -kappa0 delta T N.                                (9)

Its nonsymmetric part is nonzero; its second-row/first-column entry is -kappa0 delta/20. This does not shrink when h changes. The leading physical orientation also remains nonzero: with B=-kappa0 delta TN and D_S E=cB, it is BB*-c^2 Sym(B^2), whose (2,2) entry is (kappa0 delta)^2/400. Thus the parity mismatch concerns a nonzero physical covariance word as well as the unsquared selector.

For a literal positive Gaussian covariance clock, the conditional average of (9)'s actual source word stays nonzero for sufficiently small positive fine width by continuity of this fixed smooth fixture and bounded domination. Such widths are finite and can be charged. Meanwhile (7) is exactly zero at R_V=0 for EVERY fine width. Hence the retained-field mismatch persists beyond the zero-width algebraic label.

This is a selected-record/field-identity test, not an assertion that a whole integrated law has an order-kappa0 lower bound after a particular clock cutoff. The separate observer lower bound explains why a retained-field mismatch cannot be hidden by a marginal covariance claim.

### The named rotation-lambda cross is also zero at this selected record

This is a separate exact check, not part of the parity proof. For the existing variance-preserving terminal rotation let S=D Sprime-a sigma Vprime in the M=I case, and H_Z(s)=g(s)-g(s+epsilon Z). Its physical lambda offspring is

    L_lambda=ra[H_Z(Sprime)-H_Z(S)].

At Sprime=Vprime=0,Z=e1 in the terminal-flat fixture, H_Z has S derivative H(0)-H(epsilon e1)=0. Both copies have the same Z derivative. Therefore the COMPLETE physical first D L_lambda is zero at this record. Every zero-width mixed coefficient J_F J_L_lambda* consequently vanishes there, whereas (9) and its physical orientation coefficient do not. Finite positive fine widths converge to these different limits for the fixed smooth source.

This does not discard the full G_lambda companions or the known Gaussian constraint rows, whose joint currents still need their real readouts. It says only that the already named physical lambda mixed reserve, by itself, is not the missing selected word. Other projections, Gaussian companion observers, or new coupled-input currents are not ruled out by this test.

## 6. The explicit restricted-generator self-return

Consider a proposed endpoint completion using only these stated operations:

1. exact same-record midpoint/baseline rewrites and known linear rescalings, with all final readouts charged;
2. common-mode sources and ordinary retained mixed-covariance reserves between original-root sources and independently introduced opposite baselines;
3. signed linear combinations or independent-bank variants which preserve the declared unit Gaussian selected port;
4. affine ambient genuine-gradient lifts followed by restriction to the SAME original nonlinear feature constraints;
5. deletion only within the actual energy-plus-floor allowance.

These operations do not include a new root-dependent cross-port current, a new nonlinear constraint-law kernel, or an unproved full-gradient completion on the original tape.

Under this generator list, the opposite-word obligation has a self-return:

- The exact rewrite gives (N_h,mu_h,2sigma).
- The common mode is absent on the terminal-flat test.
- The closed ordinary mixed fields have (7), not the required field (9).
- Retaining the selected baseline port does not make its unmarked body small; the companion baseline note proves its unit-port energy floor.
- Restricting an affine ambient lift again requires the same nonlinear ancestor input and its actual derivative. Equations (1)-(3) recover the same normalized coefficient, not a smaller one.

Thus no strict positive energy/word-grade gain is supplied for this edge by these operations. A type-weight proof cannot assign a positive gain merely because the coordinate mu_h became small or a source was lifted to more dimensions: the emitted constrained readout remains of the same physical grade. This is a restricted-generator obstruction, not a lower bound against arbitrary finite VALUE programs.

The marked mixed-residual branch proved earlier remains closed and useful. Its energy contraction alone does not remove this self-return because its required partner covariance field has not been identified with (2). In particular, the parity obstruction does NOT cover using a newly generated physical residual as a partner after that residual acquires an auxiliary V input column, nor a finite decomposition into different projected-gradient offspring. Such an attempted escape needs its exact selected-field identity, genuine source/ownership contract, and all amplification/caller/energy rows; these are not supplied merely by the existing residual energy estimate. No assertion is made that a law-comparison discrepancy is an executable residual.

## 7. Cost, precision, and the next actual gate

Every displayed source is finite and uses a fixed number of original VALUES per call. At fixed rank, whole-source repetitions and polynomial-logarithmic admitted clock/filter lists have extra query/sweep exponent b=0. Every changed source argument rebuilds its original record. Distinct fine banks required by conditional products remain distinct.

The widths h and sigma, the input multiplier a/(2h), outer multiplier 2h/a, original-root first kappa0/h, actual caller paths, retained root identities and all known-row normalizations are explicit. If h=A^delta, these are real heat powers in amplitudes and precision requirements; they are not hidden in a rank constant. They do not by themselves introduce an empirical replica count, but they also do not prove a rank gain.

An admitted child with physical VALUE error epsilon_N contributes at most 2sigma epsilon_N through (1); an incoming approximation error epsilon_mu contributes at most 2r h epsilon_mu. Any root-dependent current has its own derivative/readout bill. All floors are absolute and chosen from a positive declared budget at these actual amplifications. A finite clock's retained-field approximation, cardinality and inverse output shares remain input gates.

A successful next rule must genuinely join the baseline's selected Gaussian port to mu_h while retaining the original roots, or provide a matched observer current that creates the missing V column and controls all descendants. A constant covariance or an ordinary same-root mixed field with a partner ignoring V is insufficient. The present note returns the exact finite obligation and a reusable falsification test; it does not claim the general all-rank compiler.

## 8. Diagnostics

check_midpoint_retained_parity.py passes 16,832 checks of the literal midpoint VALUE identity, its complete W Jacobian chain rule, the normalized baseline's odd/even derivative blocks, exact zero retained mixed fields at zero auxiliary coarse root, and the nonzero terminal-flat physical orientation at sixteen strictly positive fine-clock choices. The numerical identities use the complete original source at every changed query. The results are in midpoint_retained_parity_checks.json. The exact parity proof supplies the theorem; the diagnostics are not a replacement for a retained continuum-clock or all-order law proof.
