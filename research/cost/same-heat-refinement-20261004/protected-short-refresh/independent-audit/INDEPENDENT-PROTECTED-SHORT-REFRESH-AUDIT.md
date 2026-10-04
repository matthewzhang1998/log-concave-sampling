# Independent audit of a same-heat short protected refresh

Date: 2026-10-04. Bounded independent source audit; new work rather than recovery of a historical theorem.

## Result and exact boundary

**PASS, conditional on the separately supplied full-graph normalization:** one appended short Hamiltonian refresh supplies the required late protected row at a heat exponent independent of the requested posterior order R. For weak J=29/10, take h=a^(19/30), so the protected variance is gamma=sin²(h) comparable to a^(19/15). Projecting its fresh Gaussian L out of nonterminal short-force queries, only in the auxiliary retained/proxy graph, has force price O(r a^(29/10)). The emitted finite sampler retains its complete L dependence.

The exact short kernel preserves the SAME posterior Q_(a,y). Its finite implementation does not preserve Q exactly; its quadrature and Picard VALUE errors must be added to the law bill. The short quadrature count is max(1,a^(-(R-17/5))), below the main quarter-flow count max(1,a^(-(R-3/2))). The short duration, projection grade, and beta exponent do not depend on R.

**Important unresolved premise belonging to the parent audit:** an unweighted stack of N short query rows generally has restricted Gaussian norm O(sqrt(N) h), not O(h). The pointwise O(h) rows and the kernel-weighted projection error proved below are uniform in N. Importing LOW30's signed-gradient embedding additionally needs a weighted global Hilbert-space normalization controlling the complete graph, its transpose, its full affine-row norm, and its physical terminal row. This audit does not silently identify row-sum stability with that stronger certificate.

This is not a proof of the whole hidden family, a sublinear-cost compiler, a complete numerical/source batch, or an all-caller host admission.

## Source contracts inspected

1. `SAME-HEAT-POSTERIOR-REFINEMENT-AND-FLOW-GATE.md`, in the parent directory tree: exact same-heat quarter flow, finite product-integration/Picard graph, original-query count, and remaining source gates.
2. `no-copy-rank-20261004/host-observer-audit/ACTUAL-CW7-WEAK-E-HOST-SEPARATOR-AND-COMPARISON-ORDER.md`: law comparison before the physical T-only host, but retained/proxy reconstruction from the actual graph and actual private tapes.
3. LOW30, `research-source/High Acc Ideas/ai-bucket/30_low_acc.tex`, lines 741–865 and 7240–7445: normalized return, moving anchors, retained graph, actual protected projection, signed twin, and full center port.
4. LOW30 lines 420–456 and 867–904: finite physical columns/rows proved from literal iterates rather than differentiation of a small VALUE error.
5. `exact-slack/cw7/01_retained_hidden_force_certificate.md`, sections 2–6: source-level retained-state errors, original full record, actual caller first, and late protection.
6. `p-weak-mean/WEAK-29-ACTUAL-REMAINDER-AND-RETAINED-MEAN-RECONSTRUCTION.md`, sections 1–3: same-version weak-J batch, common origin, actual-versus-retained error, and separate finite numerical restoration.

All quoted source paths are local research inputs. No external action or push was taken.

## 1. Exact short kernel at the original caller

Fix y before drawing fresh L~N(0,I), independent of the entire incumbent record conditionally on the captured caller. Let X be the completed main quarter-flow output at this SAME heat a and caller y. Evolve

    x''(t)=-(x(t)-y)-a grad V(x(t)),
    x(0)=X, x'(0)=sqrt(a)L, 0<=t<=h.

The invariant phase density is proportional to

    exp[-V(x)-|x-y|²/(2a)-|v|²/(2a)].

The vector field is divergence-free and conserves the displayed Hamiltonian. Under the inspected C2 bounded-Hessian assumptions its flow exists and has the regularity needed for this statement. Thus an exact Q_(a,y) initial position remains Q_(a,y) after discarding the final velocity. The refresh must be centered at y, not at the random X.

Variation of constants gives

    x(t)=y+cos(t)(X-y)+sqrt(a)sin(t)L
         -a integral_0^t sin(t-s) grad V(x(s)) ds.       (1)

If two inputs share L, let K_h=1-cos(h) and alpha=a K_h. For h<=1 and alpha<1, the integral inequality gives

    sup_(0<=t<=h)|x(t)-x'(t)| <= |X-X'|/(1-alpha).

Therefore this extra exact step multiplies an inherited W2 error by at most 1/(1-alpha)=1+O(a h²). It neither loses the quarter-flow posterior grade nor provides an additional order improvement asserted here.

The exact stationarity statement is not transferred to the finite integrator without a VALUE error comparison.

## 2. Explicit finite short graph and uniform projection bound

Use deterministic t_i=i h/N and nonnegative product-integration weights

    w_ij=cos(t_i-t_(j+1))-cos(t_i-t_j),  j<i,
    sum_(j<i) w_ij=1-cos(t_i)<=K_h.

Execute M Picard layers using one stored main endpoint X. All old calls, modes, anchors, and zero trajectories precede L and ignore it. An uncentered equivalent recurrence is

    x_i^[0]=y+cos(t_i)(X-y)+sqrt(a)sin(t_i)L,
    x_i^[k+1]=x_i^[0]-a sum_(j<i) w_ij grad V(x_j^[k]).  (2)

Centering around a finite mode or the actual zero trajectory changes the bookkeeping but not the following Lipschitz comparison. Finite mode residuals and mode VALUE errors must be retained or charged exactly as in the parent quarter-flow construction.

Construct the auxiliary graph by deleting the direct L row from every force-query path point in (2), while retaining sqrt(a)sin(h)L at the final physical endpoint. This is a graph modification; it is not a modification of the emitted sampler. All other records and rows remain identical.

Let D_k be the largest pathwise state difference at the nonterminal queries after layer k. Bounded Hessian and positivity give

    D_0<=sqrt(a)sin(h)|L|,
    D_(k+1)<=sqrt(a)sin(h)|L|+alpha D_k,
    D_k<=sqrt(a)sin(h)|L|/(1-alpha).

The two terminal endpoints have the same free harmonic row, so only the force sum differs:

    |X_full-X_prot|
      <=sqrt(a) a K_h sin(h)|L|/(1-alpha).             (3)

This exact finite bound is uniform in N and M. In particular

    ||X_full-X_prot||_p <=C_p sqrt(Dim) sqrt(a) a h³,
    ||F_full-F_prot||_p <=C_p sqrt(Dim) r a h³,        (4)

for an original-gradient terminal force with readout r/sqrt(a). Projection changes no value whenever L=0, even at nonzero old private tape or varying y. This gives a genuine zero-origin property for this new deletion; it does not repair a previously nonzero main retained-origin debt.

No differentiability of the Hessian is used in (3). The bound uses only original-gradient Lipschitz continuity.

## 3. Actual protected row identities

Suppose the main retained endpoint has known leading private row P_old with P_old P_old*=I. The appended harmonic endpoint has

    P=(cos(h)P_old, sin(h)I),   P P*=I.

Let Pi be the orthogonal projection onto the new L coordinates, zero on all old private and exposed external-center coordinates. Then

    P Pi P*=gamma I,  gamma=sin²(h).

After the auxiliary projection all nonterminal direct query rows C_N obey C_N Pi=0. The terminal direct row retains P. Consequently

    R0=P Pi/gamma,
    R0 C_N*=0,  R0 P*=I.                             (5)

Earlier queries had zero L columns before any deletion. Later short queries have original-coordinate direct L row sin(t_i)I, with norm at most h. Their actual total terminal force coefficient is a K_h=O(a h²). All recursive responses to changed force values were already included by the finite geometric bound in (3). No old row acquires an L column just because later computations use it.

For h<=1, sin(h)>=5h/6, hence

    (25/36) a^v <=gamma<=a^v,  v=19/15.              (6)

The original terminal-force readout is allowed to read L. Projecting that terminal row too would destroy (5) and would instead cost O(r h), which is not the intended construction.

## 4. Twin width and its conditional gradient import

Assume the separate weighted graph certificate supplies the exact setting of LOW30: a full original-gradient graph with bounded affine row, sufficiently small symmetric-twin edge norm, and physical terminal predecessor norm O(a). The old main predecessor paths still cost O(a); the appended short paths add only O(a h²).

With epsilon=a^(J-2)=a^(9/10), define

    R=[P Pi/gamma, -I/epsilon],
    beta=(gamma^(-1)+epsilon^(-2))^(-1/2), B=beta R.

For C_tw=[[C,0],[C,epsilon e_t]], equations (5) give exactly

    R C_tw*=[e_t*,0],
    C_tw B*=beta(e_t,0)*,
    B B*=I.

The protected VALUE price is r a h³=r a^(29/10). The imported signed-twin VALUE price is O(r epsilon a²)=O(r a^(29/10)). Moreover

    (5/sqrt(61))a^(9/10)<=beta<=a^(9/10),
    rho=O(r/beta), and rho=O(a^(1/10)) at r=a.       (7)

Thus the width power is independent of posterior target R. The exact signed fixed-point lift is genuinely gradient only under the full weighted symmetric-graph hypotheses, not merely under the componentwise comparison (3).

The short duration has effective protected heat a h²=a^(34/15), agreeing with the weak-J exponent. Unlike the serial-history construction, this new final block explicitly creates the fresh increment; no old-history cutoff is needed merely to prove this increment exists.

## 5. Law accuracy, finite query cost, and no double counting

After analytic recentering at the posterior mode, assume the inherited fixed-moment position and velocity profiles are O(sqrt(Dim)sqrt(a)). On the short interval the exact path has that same speed profile. Replacing each force by its left grid value has residual at most

    C_p sqrt(Dim) sqrt(a) a K_h h/N
      <=C_p sqrt(Dim) a^(3/2) h³/N.

Resolvent stability contributes 1/(1-alpha). Starting Picard at the harmonic path yields the tail

    C_p sqrt(Dim) sqrt(a) alpha^(M+1).

Consequently the short implementation has the safe local strong bound

    C_p sqrt(Dim)[a^(1/2)(a h²)^(M+1)
                         +a^(3/2)h³/N]
      +mode error+numerical error.                  (8)

For each fixed R choose fixed M large enough and

    N_short>=C_R max(1,a^(-(R-17/5))).               (9)

The exponent 17/5=3/2+3v/2 is correct. Every layer evaluates and stores the N interior original gradients once. Endpoint force and necessary origins are additionally charged; no nested old posterior is called by a short-force query. Original directional first/adjoint sweeps follow these same saved points with original HVPs.

The old main endpoint is one captured input. The new count is additive, bounded by a fixed-order multiple of N_short plus actual zero/mode work. The main finite quarter-flow count already has exponent R-3/2, which dominates R-17/5. Dense known kernel multiplication may still cost O(M N²) arithmetic and storage; this audit concerns the stated original-query exponent, not a claim that dense arithmetic is free.

The finite discretization error (8) is a posterior-law comparison with the exact short flow. It is not automatically an actual-versus-retained force deletion. For the new proxy, retaining the literal same finite graph and then applying (3) avoids charging the quadrature error a second time as a fictitious same-record source error. Any separate old-state retention and numerical-proxy restoration errors remain real.

## 6. Actual firsts, caller, and origin

These conclusions use literal finite recurrences and original HVPs, not derivative convergence inferred from (8).

If the main endpoint satisfies on its full actual record

    D_W X=sqrt(a)P_old+O(sqrt(a)a),
    D_y X=I+O(a),

then differentiation of (2), positivity, bounded Hessian, and alpha<1 give

    D_(W,L) X_new=sqrt(a)P+O(sqrt(a)a),
    D_y X_new=I+O(a),                                (10)

uniformly in N at fixed finite construction depth. In the center recurrence the exact direct term is I+cos(t_i)(D_y X-I); it is incorrect to drop the derivative of y or of an old moving zero computation.

For the terminal original-gradient force, (10) gives physical first O(r), caller first O(r/sqrt(a)), and full recorded P-square curl O(r a). The leading square lift r P* Hess(V)P is symmetric on the COMPLETE private record. Caller labels have their separate declared rows.

Adjoin U=(y-y_ref)/sqrt(a) as a full label BEFORE forming the gradient embedding. This label was fixed before L and contains no L. Equations (5) remain true with zero padding on U. Same-version zero paths and moving anchors are actual original-force computations with full caller derivatives. Neither an approximate mode nor a stored zero is declared constant under y differentiation.

For the finite signed twin, LOW30's literal recurrence

    J_(k+1)=H_(k+1)(C_tw+S J_k), J_0=0

and C_tw B*=beta e_terminal^+ give

    ||J_k B*||<=beta/(1-||S||).

The terminal weak-row identity separately gives the physical row bound. Thus both physical Jacobian sides hold at actual finite iterates, even though the finite lift need not itself be an exact gradient. No small derivative bound for the difference between finite and exact lifts is asserted.

If the prior main actual-to-retained state error is e_ret, the short recurrence transports it by at most 1/(1-alpha). Therefore the new total retained-source error includes that transported e_ret PLUS (3) and any numerical proxy/origin debt. The new refresh cannot erase an inadequate old retained grade. The parent's quarter step supplies any claimed extra attenuation of its incumbent error.

## 7. Diagnostic tests and falsification boundaries

`check_protected_short_refresh.py` writes `protected_short_refresh_checks.json` and checks:

- exact rational exponents v=19/15, J=29/10, short quadrature base 17/5, beta exponent 9/10, and final r=a radius exponent 1/10;
- 12 quadratic stationary covariance identities for the exact short refresh;
- exact quadratic full-versus-protected pathwise coefficients against (3);
- 27 nonlinear finite DAG cases with Hessian in [1/2,3/4], positive exact kernel weights, the N-independent projection bound, and zero-L identity;
- the coisometry and all twin/protected row identities with a full center port;
- 20 finite-column recurrences with varying Hessians, so the checks do not quietly assume one constant derivative matrix;
- the actual sqrt(N) growth of a naive unweighted restricted row stack.

The quadratic check also verifies that replacing the emitted sampler by its protected auxiliary map generally changes its variance. For V(x)=lambda x²/2, exact momentum coefficient is sin(sqrt(1+a lambda)h)/sqrt(1+a lambda); protection replaces it by sin(h). The two coefficients differ. The intended proof never makes this sampler substitution.

These are bounded diagnostics, not substitutes for the analytic proof or the independently supplied weighted graph normalization.

## 8. Stopping conclusion

The new late refresh closes the protected-increment existence and exponent issue without tying v or beta to R. Exact same-heat invariance, explicit source chronology, zero-preserving auxiliary projection, actual finite firsts, and the short quadrature exponent pass this audit.

Full source admission still needs the separately assigned polynomial-node weighted normalization, all old retained source obligations, actual caller-domain/profile closure, finite mode/origin and source-version bookkeeping, and every downstream source occurrence at its actual scale. No claim about all of those gates, or about asymptotically sublinear original-query cost, follows merely from this successful late-protection construction.
