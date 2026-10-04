# Opposite baseline: full-tape admission and an exact linear-return obstruction

Publication copy: nonmathematical provenance wording and/or local paths were sanitized. Original/public SHA-256 values are recorded in `INVENTORY.json`; historical source/audit pins refer to the original versions.

New bounded result, 2026-10-04. This note investigates the opposite-baseline source, not the separate ambient endpoint-current calculation. It supplies two rigorous obstructions for explicit operations, rather than an impossibility theorem for every VALUE-only repair.

## Result first

1. An opposite baseline is a gradient in its fresh auxiliary variable. That fact does **not** admit it as a bounded-full-first gradient on the complete original-plus-auxiliary tape. At every fixed positive heat width, a smooth SAME-potential fixture with original Hessians in [.45,.55] makes **every** completion having the stated auxiliary gradient block have arbitrarily large full first. The proof allows an arbitrary root-only correction to the potential and unread padding. It is therefore stronger than failure of one proposed completion.
2. The new terminal-flat fixture has an exact globally linear second output. On that fixture, normalized opposite-baseline channels at arbitrary widths carry the same nonzero Gaussian row. Signed same-bank linear combinations retaining its selected word cannot reduce that row's energy. For independent banks, an honest unit-coisometry input gives the same lower bound, with no replica gain. The apparent 1/sqrt(N) from deterministic coefficient averaging also shrinks the selected unit-input coefficient by 1/sqrt(N); recovering the original port costs sqrt(N).
3. The affine ambient lift is executable and does not have this higher-jet failure. Its exact physical readout still lies on p=g0(cS+sU). Freezing that constraint loses the active ancestor difference; substituting it restores the same constrained-source state. None of the displayed baseline operations yields the missing admitted small-full-gradient offspring on the same tape.

The physical first/heat tradeoff alone is **not** a general rank obstruction. Nor is the baseline energy floor an obstruction to its absorption through a marked parent if another genuine full-gradient lift and exact endpoint join are supplied.

## 1. Exact baseline current and normalizations

Write kappa0=ra (this note's literal native mass, before conservative logarithmic envelopes), h=a sigma, and initially M=I for normalization formulas. At the actual original root W=(S,U,Z),

    x=cS+sU,
    x1=x+a[gt(S)-gt(S+epsilon Z)],
    t0=S+a g0(x),  t1=S+a g0(x1),
    Fminus_h=r[gt(t0+hV)-gt(t1-hV)],
    B_h=r[gt(t0+hV)-gt(t0-hV)],
    C_h=r[gt(t0-hV)-gt(t1-hV)].

There is the exact VALUE identity Fminus_h=B_h+C_h. Define the charged normalizer

    N_h=B_h/(2 sigma)=kappa0/(2h)[gt(t0+hV)-gt(t0-hV)].       (1)

At V=0,

    D_V N_h=kappa0 T0,
    (1/(2 sigma))D_V Fminus_h=(kappa0/2)(T0+T1).

Therefore, with DeltaA=A0-A1,

    (1/(2 sigma))D_V Fminus_h DeltaA
         =(kappa0/2)(T0+T1)DeltaA,
    D_V N_h DeltaA=kappa0 T0 DeltaA.                         (2)

Combining (2) with C_h yields the asymmetric exact split

    kappa0(K0-K1)=kappa0(T0-T1)A1+kappa0 T0(A0-A1).

For a general known M, insert M immediately after every terminal Hessian. In the physical coisometric V readout the actual source is M*B_h; a free inverse of M is never implied. The M=I subclasses below suffice to disprove the proposed generic admission/deletion rules.

Equation (1)'s inverse width is part of the actual readout, not a derivative oracle or a free coefficient extraction. Selected derivatives in this note are analytical tests of finite VALUE programs. No original HVP is executed as a VALUE.

## 2. Physical full first has a real width tradeoff

Although the fresh V row is O(kappa0) for N_h, the root row is

    D_W N_h=(kappa0/(2h))[Ht(t0+hV)-Ht(t0-hV)] D_W t0.       (3)

Bounded primitive Hessians give O(kappa0/h), not O(kappa0). The same-potential fixture in section 3 attains a fixed numerical multiple of kappa0/h. This is a necessary bound for **any** realization N_h=B G with BB*=I, since Lip(N_h)<=Lip(G).

Scope is important. If kappa0=A^(1+g) and h=A^delta with 0<delta<g, then kappa0/h=A^(1+g-delta) is still a small full physical radius relative to A. A repair allowed to spend this smaller grade could conceivably work. Equation (3) only rejects assigning the original kappa0 radius uniformly, or forgetting the inverse-width bill. It does not reject all fixed-rank graded repairs.

Multiplying a source by eta reduces its first and energy by eta, but recovering its same selected current multiplies the relevant readout or mixed weight by eta^(-1). No strict grade gain follows from that algebra alone. A nonlinear marked-parent construction may have a real gain; that requires its actual source and endpoint contracts.

## 3. No bounded-full-first completion of the actual V-gradient block

### Precise admitted-frame statement

Consider d=1, M=1, 0<a<=1/16, r>0, and any fixed 0<h<=1/4. Let g0=gt=g be a smooth anchored original gradient with .4<=g'<=.6. Let

    t(W)=S+a g(cS+sU),  c=3/5, s=4/5,
    N_h(W,V)=(ra/(2h))[g(t(W)+hV)-g(t(W)-hV)].

Suppose a genuine gradient G=grad Phi on the complete tape (W,V,padding) has its actual V block equal to N_h. Padding may be arbitrary and Phi may include any root/padding-only correction. There is **no** bound on Lip(G) depending only on r,a,h, dimension and the original Hessian bounds, uniformly over these g.

This is a statement about the literal V selected frame and its known orthogonal renamings. It does not exclude an as-yet-unconstructed full-gradient realization through a different coisometry mixing original roots, nor a weak joint constrained-law comparison that avoids such a realization.

### Every completion has the same V-dependent potential difference

Let P'=g analytically; potential VALUES are not queried. Integrating the prescribed V derivative forces

    Phi(W,V)-Phi(W,0)
      =(ra/(2h^2))[P(t+hV)+P(t-hV)-2P(t)].                 (4)

Any unknown root-only correction cancels in (4). If grad Phi is L-Lipschitz, the root gradient of the left side is 2L-Lipschitz. Hence every root-direction second derivative of the smooth right side has absolute value at most 2L. This reasoning does not require the root-only correction itself to have a classical Hessian at the chosen point.

### Exact smooth SAME-potential frequency fixture

Use

    chi(u)=exp(1-1/(1-u^2)) for |u|<1, and 0 otherwise,
    eta=1/80,  zeta=h/480,
    g_k(y)=y/2+(eta/k)(1-cos(ky)) chi(8y)
                       +zeta chi(4(y-1-h)/h),   k>=48.  (5)

The two perturbations have disjoint supports. Also g_k(0)=0. The elementary envelope |chi'|<=3 follows by bounding 2b^2 exp(1-b)<=8/e<3, where b=(1-u^2)^(-1). Consequently

    |d/dy[(eta/k)(1-cos(ky))chi(8y)]|
         <=eta+48eta/k<=2eta=1/40,
    |d/dy[zeta chi(4(y-1-h)/h)]|<=12zeta/h=1/40.

In particular .45<=g_k'<=.55 globally. There is no growth of any original Hessian bound with k. At the ancestor origin and the three terminal points,

    g_k(0)=0,  g_k'(0)=1/2,  g_k''(0)=eta k,
    g_k(1+h)+g_k(1-h)-2g_k(1)=h/480,
    g_k'(1+h)+g_k'(1-h)-2g_k'(1)=0.                       (6)

Choose the **actual native root** S=1, U=-c/s, Z=0. Then x=x1=0, t0=t1=1, and the original K3 output E is exactly zero. At this root,

    t_U=a s/2,  t_UU=a s^2 eta k.

Taking the U/U derivative of (4) at V=1 gives exactly

    d_UU[Phi(W,1)-Phi(W,0)]
      =(ra/(2h^2))[(a s^2 eta k)(h/480)+(a s/2)^2 0]
      =ra^2 s^2 eta k/(960h).

Thus every completion obeys

    Lip(G)>=ra^2 s^2 eta k/(1920h).                       (7)

At fixed positive r,a,h, the right side is unbounded as k increases. Small a or a small prefactor chosen independently of k cannot cure this uniform-source failure.

For completeness, at V=1+1/8 the terminal bump is at normalized coordinate 1/2. Equation (3) gives exactly

    |D_S N_h|=(ra/h) |chi'(1/2)|/240 (1+a c/2),           (8)

which proves section 2's physical-first lower bound.

The natural full pullback computes root companions involving g_k'(x) times a terminal VALUE second difference. It would already use original HVP results as VALUES. Equation (7) proves the stronger fact that allowing those VALUES would still not give the required bounded full first. Differentiating them is neither needed by the proof nor authorized as execution.

### Finite VALUE circuits and copies do not evade this specific statement

Any fixed finite composition/linear combination of original gradient VALUES with supplied bounded known linear maps has a first bound computed from its finite graph and primitive first bounds. Fixed finite difference widths and fixed-rank copies change that bound and their cost; they do not make it depend on the unbounded unknown k. Such a circuit therefore cannot be an exact V-block completion uniformly over (5). Approximate/non-gradient replacements require separate quantitative error and source-type certificates; they are not admitted by this exact argument.

## 4. The executable affine ambient lift returns the constraint

With independent ambient coordinates Xi=(S,p,V), set

    q_plus=S+aMp+hMV,  q_minus=S+aMp-hMV,
    C_plus=[I,aM,hM],  C_minus=[I,aM,-hM],
    G_amb=r[C_plus* gt(q_plus)-C_minus* gt(q_minus)].     (9)

This is a genuine full-gradient original-VALUE source. It uses two terminal VALUES, known maps, and has full first at most r(||C_plus||^2+||C_minus||^2). Its blocks are

    S: r[gt(q_plus)-gt(q_minus)],
    p: ra M*[gt(q_plus)-gt(q_minus)],
    V: rh M*[gt(q_plus)+gt(q_minus)].

At p=g0(cS+sU), the S readout is the actual baseline B_h. That equation is a nonlinear graph constraint, not an ambient Gaussian covariance identity.

Dividing (9) by 2 sigma gives an ambient gradient lift of N_h with full radius O(kappa0/h), consistent with section 2's possible reduced-grade choice. Thus the affine ambient source itself is not being rejected on size grounds. The missing operation is its Gaussian-compatible restriction with the original active ancestor tangent retained.

The following precise restricted rewrite system is what the present obstruction covers:

- known linear changes/padding and finite copies, with their true input/output norms;
- opposite/common baseline identities and linear combinations at literal roots;
- affine primitive-gradient lifts such as (9), followed either by freezing their nonlinear coefficients as retained labels or by substituting their literal feature graph;
- completing the literal auxiliary V-gradient block on the complete tape;
- the already admitted mixed-reserve service only when its exact complete-tape full-gradient source and SAME-root selected covariance field premises have actually been supplied.

Under these operations the baseline's active ancestor constraint has no new closing rule:

- Freezing p makes (9) an admissible conditional affine gradient, but removes dp/dx=A0. The analogous p-minus label removes (A0-A1)/sqrt(2). Keeping a label is not the same as selecting its active derivative.
- Substitution restores dp/dx and the missing same-root word, while the active source is again the nonlinearly constrained baseline. Its generic full-gradient claim does not survive substitution.
- A literal full V-gradient completion is excluded uniformly by (7), even after root-only corrections.
- Baseline subtraction returns exactly C_h, whose selected terminal difference vanishes on the terminal-flat record below.
- The mixed-reserve rule does not itself establish a premise for this constrained baseline. Already admitted lambda offspring can be processed, but the premise for this different source remains absent.

This is a self-return/type obstruction for the explicit rewrite rules, not a syntactic impossibility theorem for arbitrary original-VALUE algorithms. In particular a new known-coisometry full-gradient lift, a joint constrained-law generator, or a proved exact field identity can escape it. Simply relabeling the same constrained source or its physical row cannot.

## 5. Terminal-flat SAME-potential test and the exact bank invariant

Use exactly section 4 of RETAINED-MIXED-KERNEL-AND-TERMINAL-FLAT-GATE.md:

    T=[[1/2,1/20],[1/20,1/2]],  N=e1e1*, delta=1/40,
    q=a epsilon, epsilon=a^(9/10), r=a,
    g(y)=Ty+delta q psi(y1/q)e1,
    psi(z)=(z+1/2)chi(10(z+1/2)).

At S=U=0,Z=e1,

    T0=T1=T,  A0=T, A1=T+delta N,
    common selected word=0,
    missing opposite word= -ra delta T N.               (10)

Its skew norm divided by ra is sqrt(2)/800, independent of a. The same potential's second output is **globally**, not just locally, linear:

    g_2(y)=(1/20)y1+(1/2)y2.

Since t0=0, for every positive width sigma and every V,

    [N_(a sigma)(V)]_2=ra[(1/20)V1+(1/2)V2].             (11)

The common channel's second output is independent of V, so conditional centering removes its entire V row. It cannot cancel the row in (11). Formula (11) also holds at the same S=U=0 root with Z=0, where the original E is exactly zero. Thus this is an actual baseline body, not an E-marked residual.

### Same bank, arbitrary signed scalar combinations

Let Y=sum_i alpha_i N_(h_i)(V), allowing arbitrary finite widths and signed scalar coefficients. Preserving (10)'s nonzero second-row selected ancestor word forces sum_i alpha_i=1. Therefore

    Y_2=ra[(1/20)V1+(1/2)V2],
    ||Y_2||_2=ra sqrt(101)/20.                           (12)

No width selection, signed cancellation or fixed number of same-bank copies improves this body energy. An overall readout scale is included in alpha_i; restoring the current restores (12).

### Independent banks, with the actual selected input kept unit-normalized

Let V_i be independent standard 2-vectors and define the physical incoming port

    V_in=sum_i beta_i V_i + an independent Gaussian fill,
    sum_i beta_i^2<=1.

For Y=sum_i alpha_i N_(h_i)(V_i), the second-row covariance/selected coefficient against V_in is the original row multiplied by sum_i alpha_i beta_i. Preserving the same selected port requires

    sum_i alpha_i beta_i=1.

Bessel/Cauchy-Schwarz gives sum_i alpha_i^2>=1. Hence exactly

    ||Y_2||_2=ra sqrt(101)/20 sqrt(sum_i alpha_i^2)
              >=ra sqrt(101)/20.                       (13)

There is no improvement even when the finite bank count N is large. The statement extends to correlated banks by expressing them in an independent basis and charging the actual unit input row. It is a same-port result; no fresh independent original root is substituted for the observed root.

### Port-only Bessel bound beyond scalar combinations

If an alternative construction preserves only the actual missing second-output/first-input Gaussian covariance port, write its unit standard incoming coordinate as V_in,1. The exact requirement at a retained root is

    E[Y_2 V_in,1 | W]=ra/20.

Conditional Cauchy-Schwarz gives ||Y_2||_(L2|W)>=ra/20. This weaker numerical bound permits arbitrary matrix adapters, unread padding, bank counts, and even nonlinear processing: the preserved Gaussian covariance coefficient itself costs that much energy. It applies to an actual covariance/selected Gaussian port, not merely a pointwise derivative at V=0. It also allows the coefficient to be carried by a separately admitted known Gaussian carrier, provided that carrier and its joint current are accounted for. It does not prohibit the marked-parent absorption discussed below.

### Different interface: deterministic coefficient averaging

If the requested target is instead the arithmetic mean of separately reported deterministic derivatives, one may impose sum_i alpha_i=1. Equal weights give body energy ra sqrt(101)/(20 sqrt(N)). But against the genuine unit input V_in=sum_i V_i/sqrt(N), the selected row is also reduced by 1/sqrt(N). Restoring the old port requires sqrt(N) readout and returns (13).

If only deterministic coefficient averaging were needed, obtaining a factor A^eta in body energy would require N>=A^(-2eta) up to numerical constants, a real original-query/first-sweep replication exponent. Fixed-rank finite N gives at most a numerical factor. It is invalid to use the deterministic-mean normalization as if it preserved a unit Gaussian input.

### What this invariant does not rule out

Equations (12)-(13) cover the displayed linear baseline-copy algebra, with its actual selected port. They do not rule out a nonlinear stationary kernel or marked-parent absorption. Nor do they prohibit splitting off a truly supplied known affine carrier and executing its joint current by a separate admitted adapter. Removing that carrier without retaining its current would remove the nonzero second row of (10). Such a supplied carrier/current executor would be new information beyond the linear deletion rule tested here.

## 6. VALUE, first, caller and floor ledger

All operations are finite and all amplitudes below are actual readouts.

- A standalone B_h or N_h evaluation needs g0(x) and its two changed terminal VALUES: three original queries. It ignores Z and does not need to execute the unused feedback branch.
- Fminus_h alone uses the complete six-query native graph with its two final terminal arguments changed. Adding B_h shares the q_plus VALUE and needs one additional terminal VALUE, giving seven distinct generic queries. If the original E is also evaluated, two unheated terminal VALUES are retained, giving nine generic queries in the joint E/Fminus/B evaluation. Exact aliases may reduce these only when semantic keys agree.
- N_h's literal original-source numerical error and direct caller are multiplied by |1/(2 sigma)|. Its actual root first includes the kappa0/h term. New conditional innovation widths add their own 1/v and incoming-first factors; they are not canceled by merely describing the old roots as retained.
- An ambient (9) query has two primitive terminal VALUES. Supplying p=g0(x) adds the actual ancestor query and constraint; it does not turn the ambient input into a standard Gaussian. Full first/adjoint replays use the corresponding original HVP sites only for verification, with every changed nonlinear source graph rebuilt.
- A finite bank list multiplies the actual query/replay count. Unit input adapters, output coefficients, retained identity rows, known fills, and any deterministic-coefficient-to-Gaussian-port sqrt(N) recovery are all charged.
- A numerical source floor epsilon_i in an output with coefficient alpha_i contributes at most the propagated sum of |alpha_i| epsilon_i times subsequent actual incoming firsts. It remains an absolute floor. At Z=0, E=0 while (11) is nonzero, so neither this floor nor the opposite body may be absorbed into a putative original-E factor.
- A comparison error is not an executable source. The existing mixed reserve's useful residual energy gain remains valid for genuinely admitted offspring; these statements do not authorize feeding a coupling discrepancy back as a VALUE or integrating away the original root observer.

## 7. Diagnostic result and remaining positive obligation

check_opposite_baseline_return.py passes 102,406 finite checks over 45 scalar completion fixtures, four terminal-flat native scales, arbitrary signed widths, and independent-bank normalization examples. It checks the original-Hessian sandwich, exact terminal VALUE/first cancellations, the full-potential difference's root Hessian, the physical kappa0/h row, the same-potential terminal-flat word, and the Gaussian bank coefficients. Its result is opposite_baseline_return_checks.json. Analytic equations (4)-(13), not numerical experiments, prove the stated obstructions.

No finite VALUE identity reducing the complete terminal-sum/ancestor-difference current to the existing small-full-gradient offspring has been obtained here. The substantive missing source gate is a full-tape/joint-law realization that preserves the nonlinear ancestor constraint and its covariance field. A genuine new realization can still be absorbed by the retained mixed kernel. The current baseline subtraction, conditional V-gradient observation, width rescaling, and finite linear bank manipulation do not supply it.
