# Same-heat posterior refinement: a concrete flow telescope and its remaining gate

Date: 2026-10-04. New mathematical research.

Publication copy: nonmathematical provenance was sanitized and the independent audit's scope clarifications are incorporated below. Original/public hashes are recorded in `INVENTORY.json`. The independent audit identifies its original source pin, not these edited bytes.

## Result first

There is a concrete same-heat refinement mechanism that really targets the posterior, rather than merely its old force mean: a **quarter-period exact Hamiltonian refresh**. It uses the incumbent once, contracts its law error by O(a), and introduces only original-gradient calls afterward. Its exact endpoint has the right new carrier, actual center first I+O(a), and attenuates old retained-state errors by an additional a.

This is not yet a zero-power or sublinear-power source compiler. The exact force integrals are the missing executable operation. A fully explicit finite product-integration/Picard DAG implements the mechanism, but its safe C2 law-level/query certificate for the displayed approximation is

    Q_R(a) <= Lambda_R [Q_P(a) + a^{-(R-3/2)}],
    R=P+1/2,
    c_R <= max(c_P,P-1).

For that displayed finite approximation, the old-call cost is additive rather than multiplicative/reheated, but the new forcing still grows linearly in order. The protected/retained-host join remains open, so this is not an admitted recursive-family improvement. Fixed Picard depth does not remove the displayed quadrature bill. The scalar C2 Hessian-sandwich bump fixture is a deterministic oracle-information lower bound only for the first-Picard integral at a fixed input, including adaptive queries whose baseline transcript is fixed. It is not a lower bound for the exact nonlinear flow, Gaussian-input Lp strong approximation, randomized methods, or posterior-law compilers, and it does not prove the upper exponent P−1 necessary.

The useful new target is therefore sharply defined: replace the displayed finite flow DAG's low-regularity quadrature by a **complete positive same-endpoint law compiler**, with actual old-call reuse and sublinear heat exponent, while returning the protected/retained/caller source ports. This is a separate candidate from Gaussianizing an old force source more accurately.

## 1. Inspected source boundary

Inputs inspected:

- `cost/recurrence-audit-20261004/FIXED-SEED-RANK-CORRECTIONS-AND-SHARING-AUDIT.md`, especially Candidate B.
- `exact-slack/cw7/01_retained_hidden_force_certificate.md`.
- `exact-slack/cw7/02_same_version_batch_and_numerical_census.md`.
- `p-weak-mean/WEAK-29-ACTUAL-REMAINDER-AND-RETAINED-MEAN-RECONSTRUCTION.md`.
- `cost/COST-OPTIMIZATION-PROOFS-RECONSTRUCTED.md` and `cost/GROWTH-AND-SUBCRITICAL-TARGETS-RECONSTRUCTED.md`.
- `exact-slack/exact32-proof-v2.tex`, for the C2 Hessian-sandwich/original-gradient/directional-HVP model and terminal-versus-source distinction.

The CW7 and weak-E identities legitimately preserve the conditional mean of their literal finite old force: E F_actual = E H_fin + E E_fin. They do not identify that mean with the true posterior force to arbitrarily high order. Their retained deletion and origin certificates likewise do not provide a posterior-bias correction. The missing historical Stein note is not promoted to an available theorem.

Normalize the upper Hessian bound to one. Fix the entering caller y before all new private sampling, and let

    pi = Q_(a,y),  pi(dx) proportional to exp(-V(x)-|x-y|^2/(2a)) dx.

Suppose the incumbent X0=R_P(a,y;W0) is one complete finite source call, with W2 error at most Lambda sqrt(Dim) a^P and the inspected source interfaces. Every W0 alias and captured outer label remains unchanged.

## 2. Exact same-heat posterior kernel

Sample one new Z~N(0,I), independent of W0 conditional on the exposed caller. Solve, for T=pi/2,

    x''(t) = -(x(t)-y)-a grad V(x(t)),
    x(0)=X0,   x'(0)=sqrt(a) Z.

The potential and kinetic energy are

    H(x,v)=V(x)+|x-y|^2/(2a)+|v|^2/(2a).

The flow conserves this energy and phase volume. Thus if X0 has law pi, independent v0~N(0,aI), then x(T) again has law pi. No Metropolis step or acceptance variable is involved in this exact identity. C2 with bounded Hessian supplies a globally Lipschitz vector field and its C1 finite-time flow.

The variation-of-constants identity is

    x(t)=y+cos(t)(X0-y)+sqrt(a) sin(t) Z
         -a integral_0^t sin(t-s) grad V(x(s)) ds.

At the endpoint it becomes the concrete source identity

    x(T)=y+sqrt(a) Z
         -a integral_0^T cos(s) grad V(x(s)) ds.             (H)

All path points are at the original physical scale sqrt(a). There is no smaller-heat old posterior call in (H).

### Old-law error is genuinely attenuated

Run two paths from x0 and x0' with the same new Z. The Hessian sandwich gives a matrix H_s with 0<=H_s<=I and

    d(t)=cos(t)(x0-x0')-a integral_0^t sin(t-s) H_s d(s) ds.

For a<1, because the positive kernel has row mass at most one on [0,T],

    sup_t |d(t)| <= |x0-x0'|/(1-a),
    |d(T)| <= a |x0-x0'|/(1-a).                            (C)

Couple the incumbent to an exact posterior initial state, independently of Z. Invariance and (C) prove

    W2(Law x(T),pi) <= [a/(1-a)] W2(Law X0,pi).           (L)

An approximate endpoint x_hat(T) with actual coupled error eta therefore has

    W2(Law x_hat(T),pi) <= Lambda sqrt(Dim) a^(P+1)+eta.

This corrects posterior bias and covariance together. It is stronger than a claim about the old force's mean.

## 3. Caller, retained, carrier, and chronology consequences

These are exact-flow statements, and also finite-graph first bounds for the positive-weight DAG in Section 4. They do not supply all native source ports by themselves.

### Actual caller improvement

The endpoint derivative in the starting position is O(a) by the same integral inequality. Consequently the incumbent is evaluated and differentiated once. A finite forward/adjoint sweep through its stored record receives the accumulated outgoing direction once; discarded records need paid replay.

If the incumbent's actual center derivative is I+O(sqrt(a)), as for the CW7 seed, differentiating the physical integral equation gives

    D_y x(T)=I-a integral_0^T cos(s) Hess V(x(s)) D_y x(s) ds,
    D_y x(T)=I+O(a).

The old actual-versus-kept g distinction is respected: this is a new endpoint calculation, not replacement of an actual derivative by a kept one.

### Old retained-state error gains one a

Replace X0 by its same-record retained state X0_ret and use the same Z. The endpoint error is at most a/(1-a) times the entering retained-state error. Thus an entering state deletion sqrt(a) a^K returns sqrt(a) a^(K+1), before independently assigned flow/numerical deletions. Origins are computed by executing this same flow graph on the old source's actual zero trajectory and Z=0.

This does not permit retaining an incorrect center elsewhere in a later hidden decoder. Every attached host still needs its actual kept-center transport and source-domain certificate.

### The complete leading carrier is explicit

If the entering retained state has leading row sqrt(a) P_old W0 with P_old P_old*=I, the harmonic row at an intermediate t is

    C_t=(cos(t) P_old, sin(t) I),   C_t C_t*=I.

At T the old row vanishes and the leading row is P_new=(0,I). All the old randomness that survives enters through the a-weighted force integral. The full endpoint derivative is

    D_(W0,Z) x(T)=sqrt(a) P_new + O(sqrt(a) a).

The old original-force readout at this retained endpoint consequently has the usual symmetric leading square lift and O(r a) full recorded curl. This uses the entire record, including old columns, and does not identify P_old, P_new, or any later proxy B.

The positive standard Gaussian input and source-zero carrier are genuine. They are **not** an independent unread additive keep: Z is read by the force path. A protected late increment/source theorem, or a justified protection operation with its actual state-error price, remains required. A lower bound on the output covariance alone cannot substitute for that chronology.

### Conditional mean and covariance use one shared endpoint record

Writing x(T)=y+sqrt(a)Z+D, one must keep

    E x(T)=y+E D,
    Cov x(T)=aI+sqrt(a) E[Z(D-E D)*+(D-E D)Z*]+Cov D.   (M)

The cross terms in (M) are not optional. Independent copies inside products are separate complete records; deterministic repeated queries inside a single flow record may be cached only with the complete matching semantic key.

## 4. A literal finite DAG and its honest query count

For clean moment bounds recenter at the posterior mode m, satisfying m-y+a grad V(m)=0. The mode can be approximated to the required fixed-target VALUE budget by the usual contracting original-gradient iteration; all finite mode and zero versions are frozen in the batch. Mode approximation error is separately charged.

Put q=x-m and f(q)=grad V(m+q)-grad V(m). Then f(0)=0 and Lip(f)<=1. Choose fixed uniform times t_i=iT/N. Define the nonnegative exact kernel weights

    w_ij=cos(t_i-t_(j+1))-cos(t_i-t_j),  0<=j<i,
    sum_(j<i) w_ij=1-cos(t_i)<=1.

Start with the harmonic paths

    q_i^[0]=cos(t_i)(X0-m)+sqrt(a) sin(t_i)Z.

Use M finite Picard layers

    q_i^[k+1]=q_i^[0]-a sum_(j<i) w_ij f(q_j^[k]).      (DAG)

At each layer evaluate the N original gradients f(q_j^[k]) once and store them. All successor readouts of that layer reuse those actual values. The preceding layer is not rebuilt for each readout. The same one incumbent X0 remains a captured input throughout. Layers and numerical schedules are fixed before caller differentiation.

The exact residual telescope is

    q_i^[k+1]-q_i^[k]
      =-a sum_(j<i) w_ij [f(q_j^[k])-f(q_j^[k-1])].     (T)

Its VALUE contraction is a. It does not assert a correspondingly high-order derivative bound on each bracket; under C2, close query points need not have close Hessians at any specified power rate.

At fixed target M, the literal VALUE count is

    Q_DAG <= Q_P(a)+M N+Q_mode+Q_zero.

Each source directional first/adjoint sweep is a constant multiple of the expanded original count, using the saved actual points. Known kernel multiplication costs O(M N^2) arithmetic if performed densely. That cost is real, even though the present exponent ledger counts original queries. A new late protected construction, if needed, must add its own executed queries, dimensions, and widths.

For an entering sqrt(a) sqrt(Dim) moment profile, bounded Hessian gives a path-speed profile of that same order. Piecewise-constant product integration therefore has the uniform strong VALUE estimate

    ||x_hat(T)-x(T)||_Lp
       <= Lambda_p sqrt(Dim) [a^(M+3/2)+a^(3/2)/N],    (E)

with numerical and mode errors added separately. One way to see the quadrature term is to compare the exact path force within each interval to its left value, use Lip(f)<=1 and the path-speed bound, and then apply the row-mass stability factor 1/(1-a). The Picard tail follows from (T).

Choose fixed M sufficiently large for target R and N at least Lambda_R a^{-(R-3/2)}. Combining (L) and (E), R=P+1/2 is valid at the law level, and the original-query certificate is

    c_R <= max(c_P,R-3/2)=max(c_P,P-1).                (B)

This proves an additive rather than multiplicative old-call ledger for this explicit approximation. It does **not** establish a complete new old-source family: retained controlled-block enumeration, protected increments, uniform actual caller domains, and numerical tape encodings remain source obligations. Even granting those ports, (B) is asymptotically linear rather than sublinear in order.

## 5. Two minimal falsification fixtures

### 5.1 Independent additive correction changes covariance already for a quadratic

Take V(x)=lambda x^2/2, y=0, omega=sqrt(1+a lambda), c=cos(omega T), s=sin(omega T)/omega. The exact map is

    x(T)=c X0+sqrt(a) s Z.

Its target variance a/(1+a lambda) is preserved exactly because

    c^2 a/(1+a lambda)+a s^2=a/(1+a lambda).

The correct additive decomposition around the fresh carrier is

    sqrt(a)Z + [cX0+sqrt(a)(s-1)Z].

Replacing the bracket's Z by a fresh independent Z' preserves its marginal mean but changes the variance to

    c^2 a/(1+a lambda)+a[1+(s-1)^2].

The covariance error is 2a(1-s)=lambda a^2+O(a^3), and the W2 error from the true Gaussian is (lambda/2) a^(3/2)+O(a^(5/2)). Thus independent-bank additive residuals do not preserve the source identity even in the easiest model. Shared same-record Z is lawful and necessary here; it is not permission to alias records that a different product contract requires to be independent.

### 5.2 Fixed-depth Picard is not arbitrary-order finite quadrature under C2

In normalized one-dimensional coordinates choose z0=0, Z=1 and g_a(z)=a h(z). At the first Picard endpoint,

    z^[1](T)=1-a integral_0^T cos(t) h(sin(t))dt
             =1-a integral_0^1 h(u)du.               (Q)

Suppose a deterministic quadrature has q queried points in [0,1], with gradients and Hessian actions available there. Add the endpoints and let ell_i be the intervening gap lengths. In each gap put

    psi(u)=ell_i r^2(1-r)^2,  r=(u-u_i)/ell_i,
    h_+(u)=u/2+psi(u),   h_-(u)=u/2-psi(u).

Set psi=0 outside [0,1]. Both psi and psi' vanish at every query point. Therefore h_+ and h_- give identical original gradient and HVP data at those points. Their derivatives remain between approximately .3075 and .6925; in particular both satisfy the strict [1/4,3/4] Hessian sandwich. A globally C2 physical potential is obtained from V_a(x)=a H(x/sqrt(a)), H'=h.

But

    integral_0^1 psi(u)du=sum_i ell_i^2/30 >=1/[30(q+1)].

The two first-Picard physical endpoints in (Q) differ by at least

    a^(3/2)/[15(q+1)].                              (F)

This is a worst-case deterministic information bound for the first-Picard integral at the fixed normalized input (z0,Z)=(0,1). Adaptive deterministic queries are covered by choosing the two admissible potentials after the baseline oracle transcript is fixed; the potentials may depend on the algorithm and a. A common deterministic answer has worst-case error at least half the displayed two-endpoint separation. No lower bound for the exact nonlinear flow, Gaussian-input Lp strong approximation, randomized methods, or posterior-law compilation follows. The upper query exponent P−1 is not shown necessary. Random law-corrected quadrature, native endpoint counterpackets, and special posterior cancellations remain possible.

## 6. A source-level Stein identity that states the alternative gate exactly

There is also a useful exact finite-source defect identity. Center at the exact mode for notation, standardize S=(X-m)/sqrt(a), and write its actual finite map as

    S(W)=P W+r(W),   P P*=I,
    b(z)=sqrt(a)[grad V(m+sqrt(a)z)-grad V(m)].

The target density is proportional to exp(-|z|^2/2-U(z)), with grad U=b and Lip(b)<=a. Gaussian integration by parts on the **complete owned W** gives, for smooth vector tests phi,

    E[(S+b(S)) dot phi(S)-div phi(S)]
      =E[(r+b(S)) dot phi(S)+(P Dr*) : D phi(S)].     (S)

Thus the actual posterior defect is a paired rank-one/rank-two source, not just the mean of r+b(S). It is computable as a formal coefficient from the finite VALUE map and its actual first/adjoint actions. The identity does not authorize treating P Dr* as a new differentiable VALUE source: differentiating its saved HVP children would request unavailable third derivatives unless a native frozen-child/descendant contract supplies the missing rows. Full matrix formation also pays all directions; (S) alone is not a cheap matrix oracle.

For S=tau W and b(S)=a lambda S, the two terms in (S) combine to

    [(1+a lambda) tau^2-1] E phi'(tau W).

The mean drift is zero for every tau, while the true covariance condition is tau^2=(1+a lambda)^(-1). This is a minimal test for any proposed mean-only posterior correction.

A same-heat Stein compiler would close the original target if it took the paired defect (S), preserved its complete source ownership, and emitted an actual positive corrected map with a strictly higher **grouped posterior-defect grade** and the inherited caller/retained/protected ports, at Q_P(a)+a^(-b_P) cost with a sublinear running envelope of b_P. Ordinary L2 contraction of a Poisson insertion is not used or asserted here.

## 7. Exact theorem that would close the flow version

For each fixed target R, prove a finite original-gradient/HVP VALUE compiler K_hat_(a,y,R)(x,Z;U) for (H), with the following complete return:

1. Its same-caller joint law is close enough to the exact quarter-flow endpoint to contribute at most Lambda_R sqrt(Dim) a^R uniformly over the admitted incumbent input profiles. An actual coupled bound would suffice; a weaker law bound needs the full observed-input contract.
2. Every old source call is the single R_P(a,y) plus at most fixed/polylog same-heat complete copies explicitly required by the new compiler. Original path queries and fixed-target independent duplication are enumerated. No old suffix at a^eta is hidden in U or numerical precision.
3. It returns actual center first I+O(a), complete physical fresh first, exact carrier/origin, retained graph and state-level deletion, a truly protected late Gaussian port with positive gap, and all hidden-decoder/captured-label caller contracts.
4. Its total original query count beyond the incumbent is Lambda_R a^(-b_R), with max_(r<=R) b_r=o(R), including every clock, source call, covariance/root correction, actual tape dimension, precision floor, and discarded-tape replay.

Then the actual recurrence is

    c_R<=max(c_P,b_R),

and the old positive cost can be diluted as intended. The exact posterior invariance and old-error contraction are already supplied by (H)--(L); the finite law/source compiler, particularly same-endpoint covariance and protected chronology, is the open part.

## 8. Checks and status

`check_same_heat_refinement.py` generates `same_heat_refinement_checks.json`. It verifies the quadratic invariance and independent-bank discrepancy, the blind-bump Hessian/integral bounds, and the nonnegative product-integration row sums. All checks passed.

Proved here: exact-flow posterior invariance; O(a) incumbent contraction; the source/carrier/caller identities within their stated scope; explicit finite DAG and its query/strong-error certificate; the quadratic covariance obstruction; the deterministic fixed-input first-Picard quadrature obstruction; exact finite-source Stein defect identity.

Not proved: a sublinear-power complete same-heat posterior source; arbitrary-order positive endpoint law compilation; a new protected late increment theorem; or a closed high-order hidden family. The tests deliberately separate those open claims from the exact algebra.


## Publication audit clarification

The independent audit (original SHA-256 ff7414e849925eeea4b204758ab43661052175144ccc68d2eddafc76a41dbecd) gives scoped PASS for the exact posterior kernel, contraction, actual first/carrier statements, legal finite VALUE/HVP DAG, and its query upper certificate. Protected/retained-host closure remains open. Two separable prefix sums reduce the same DAG's vector arithmetic to O(MN), plus original queries and the cost per vector operation; this does not improve its original-query heat exponent. See the accompanying independent audit for the exact identity and complete scope.
