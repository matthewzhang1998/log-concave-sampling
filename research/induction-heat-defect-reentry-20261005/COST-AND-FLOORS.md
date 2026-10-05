# Complete cost, clocks, precision, and replay for the marked shell

2026-10-05. Companion to HEAT-SHELL-REENTRY.md. The finite native implementation is imported, not executed here. A raw VALUE census is not a native wrapper census.

## 1. Exact raw source census

For a generated n-force history h, let v_h=sum_v 2^(k_v+1) be its ordinary active-source original-VALUE count. Its n marked telescopes cost

    sum_m [v_h+2^(k_m+1)]=(n+1)v_h.

This includes BOTH sets of same-center anchors in the marked occurrence. If the original tau family is rebuilt alongside the shell, the total raw count is (n+2)v_h. It is an exact unmerged call count before wrapper replication and exact same-key caching.

With V_n from the sealed generator,

    V_1=1,
    V_n=sum_(i=1)^(n-1) [(n-i)(i+1)V_i T_(n-i)
                       +i(n-i+1)T_i V_(n-i)],

the shell has raw census (n+1)V_n per one matched original clock tuple at each history. The common old-plus-shell construction has (n+2)V_n. Each term has exactly 2n-1 original clocks; there is no new heat-interpolation time grid.

At n=8, V_8=27,020,800:

    shell only: 243,187,200 raw original VALUES,
    rebuilt old plus shell: 270,208,000 raw original VALUES.

These large numbers are fixed-rank constants. A raw history's maximum ordinary count is 142, its shell count is 1,278, and its common old-plus-shell count is 1,420. Actual histories can be cheaper.

## 2. Complete finite native bill

Choose one positive original-clock rule fine enough for both heat levels. For M_n positive nodes per clock, a deliberately unmerged bound is

    (n+2)V_n M_n^(2n-1)

RAW VALUES before native wrappers. The actual bill is the versioned sum

    Q_new = Q_outside-old,complete(required version)
       + sum_(old or shell history,node,force occurrence v)
              2^(k_v+1) * mu_v * M_native(k_v,b_v,epsilon_v)
                    * Q_actual_original-caller/replay(v)
       + Q_coefficient-captures + Q_known-Gaussian-rows
       + Q_scalar-clocks-and-roots + Q_physical-readouts
       + Q_all-affected-ancestor-and-anchor-replays,

where mu_v=2 at a marked occurrence and mu_v=1 elsewhere. The expanded versioned DAG counts each operation once: if a caller term already includes a capture or ancestor replay, do not charge that same operation again in the separate capture/replay terms. Native M includes the complete finite source/filter/pair/calibration/mean/gate program. No tensor expectation is supplied as a query. Fresh heat increments and the coarse/private/source/filter/reset/cut roots are all charged; the same Z can be stored and reused only under the exact common-dictionary key prescribed by the proof.

The shell needs at most n additional D-dimensional heat-increment rows per common history/node key, before exact dictionary sharing. Cold/warm/marked calls reuse a stored row only if every caller, clock, source order, shield, source version, precision, and active input matches. A warm field from a different native dynamic tape is a different VALUE call even when its analytic coefficient label agrees.

All readout/covariance shares are frozen after the complete old-plus-shell group count is known. Their inverse products enter the actual scalar root coefficient and the native moment/first constants. They may be public-log losses for the fixed-rank construction; they are not assumed constant before computing the count. There is exactly one reserved independent keep; no variance is counted twice.

The earlier tau sampler is rebuilt at its newly required finite order, allocation and source/clock floors. The operation is not an incremental mutation of a previously sampled output. Reusing its weaker recorded LAW error or source version is insufficient. Unaffected exact cached keys may be reused; all changed ancestors and anchors are replayed.

## 3. Clock and source floors

Let C_n=S_n 2^(n-2)(n-2)! and D_n=(2n-1)2^(2n-1)C_n, as in the sealed generator. Multiply these by the ACTUAL telescope, normalized-mark factor, current-normalization, readout, frame and common-keep constants. A safe worst-floor coefficient envelope is

    C alpha^n sigma^(-(n-2)) sqrt(D).

The same positive dyadic original-resolvent rule, with scalar multiplier error delta_0, has finite-to-continuous allowance C delta_0 alpha^n sigma^(-(n-2)) sqrt(D) over the original standard Y. At sigma=alpha^(1/(n-1)) and target P_n=n+1/(n-1), choose delta_0<=c alpha after all displayed constants and consumer amplification are divided out.

Likewise choose each normalized local response floor at most c_v alpha, with its finite tensor telescope and downstream amplification divided out. At a marked source, normalized sign weights including its outer one-half have total absolute weight compatible with the SAME bound

    source VALUE error <=2 delta_g/(A t).

Thus a sufficient original VALUE floor is delta_g<=epsilon_value A t/2, after splitting all local preparation floors. At the smallest t this has scale const*A*alpha^(1+1/(n-1)), before any stricter wrapper/Sobolev floor. At n=8 it is const*A*alpha^(8/7). Captured center floors similarly pay their ACTUAL 1/t; an increment-root or delta-arithmetic error also pays its actual row derivative.

Use one radial pullback per field with R=sqrt(D)+s and the sealed Gaussian tail certificate; the normalized difference inherits the same half-sum response floor. Pick s from the required epsilon via the displayed logarithmic bound. Its radial Jacobian action remains O(D); a dense D-by-D matrix is never formed.

The exact curl and first statements concern the mathematical source field. Absolute VALUE rounding alone does not certify numerical Sobolev or curl accuracy. The existing finite native/source/filter Sobolev, arithmetic-first, mean, covariance-root and calibration floors remain independently required and propagated. Every source error travels through its actual strict-ancestor amplitude product. A signed root never makes a numerical error allowance negative.

## 4. Prior order and alpha query exponent

For a fixed requested final P and computed complete downstream inverse-alpha amplification J, a conservative starting native order is

    b > 4(n-1)(P+J),

because the smallest starting source exponent is beta_shell=1/[4(n-1)]. Then evaluate the literal prior

    sum_v |rho_v|^b product_(strict ancestors a) |rho_a|

at the final finite program, with all paths, old versions and amplified floors. It is a starting rule, not a substitute for that check.

The fixed-k, fixed-b serial native theorem located by the input audit gives M_native public-log dependence on fixed-power radii and numerical tolerances, including nested original-source calls. The original resolvent rule is also public-log. The marked shell adds no new quadrature dimension and only a finite twofold local source expansion. Therefore its new local inverse-alpha query exponent is zero RELATIVE TO its actual original caller/replay program, under those imported implementation hypotheses.

The complete exponent is still the max-plus cost of the fully expanded versioned DAG:

    c_shell <= max(c_outside-old,complete(required version),
                   all actual caller/ancestor replay exponents),

with any non-public-log replication, readout normalization cost or native guard correction inserted rather than hidden. The cost of the underlying requested-precision g oracle, if more than one unit per VALUE, remains a separate bit/operation cost. No uniform-in-P zero exponent, terminal arbitrary-order queue, or c(P)=o(P) follows.

## 5. Numerical geometry near alpha=1

At alpha=1 the two heats coincide and the shell is identically zero; omit it. For 0<alpha<1 compute delta=sqrt(tau^2-sigma^2) by certified scalar arithmetic, for example using expm1 on the difference of logarithms to avoid cancellation. No derivative with respect to alpha or a sampled clock is included in the fixed-caller source-first theorem. Scalar labels are frozen and their numerical errors have the separately charged interval-propagation floor. Near equality one may also omit a shell only if its proved coefficient error fits the explicitly assigned tolerance; it is not silently zeroed.
