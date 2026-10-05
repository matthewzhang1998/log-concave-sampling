# A positive VALUE repair of the first nonlinear ancestry debt

2026-10-05. A new order-three canonical-m3 component for a wider nonlinear class, with the existing matrix-free resummed shift kept in the executed program. This is not general-C2 order-four closure.

## 0. Scope and result

Assume the original anchored convex-gradient source satisfies

    g(0)=0, 0 <= Dg <= A I, 0<A<=1/2.

For the improved bias certificate, impose the following additional structural qualification. In an arbitrary fixed orthogonal decomposition into blocks of dimension d_j<=b, g acts blockwise, and each block satisfies

    Lip(Dg_j) <= B A.

Here b and B are declared constants, independent of A and ambient dimension D. The decomposition is analytical: the executing algorithm never reads it. This includes arbitrary anisotropic quadratics, all rotated coordinatewise smooth nonlinearities with this curvature bound, and fixed-size genuinely noncommuting nonlinear blocks. It excludes the uncontrolled high-dimensional radial examples. The source ports and finite program below are valid under the original C2 assumptions even when the improved bias qualification fails.

Let

    k=sqrt(3/8), d0=1+1/sqrt(2),
    u_A=1+A+1/sqrt(2)+A k,
    C_(b,B,A)=1+k+(B/2)(u_A^2+d0^2)sqrt(b+2).

There is a finite original-VALUE source with 4D private Gaussian coordinates and (5+3J)N_out VALUES per residual occurrence, where J=O(log^2(1/delta_br)), whose own-mean completion satisfies

    integrated W2 error to N(m3(Z),I)
      <= [C_(b,B,A) A^3 + delta_br A^2
          +delta_out A C_F + Lambda A^4] sqrt(D)
          + restored absolute floors,

    C_F=1+(1+sqrt(2))A+k A^2.

Here m3 is the genuine-history canonical target from the pinned packet, and Lambda denotes the collected guarded completion constant, including the fixed numerical factors below rather than an unchanged sharp native prefactor. Taking delta_br<=A and delta_out<=A^3 gives order-three target accuracy, while the complete own-mean consumer costs order four. Fixed-grade original VALUE inverse-A exponent remains zero. All native compiler guards remain imposed; the safe normalization below gives private first <=4A and curl <=12A^2. This is an actual A^2-to-A^3 nonlinear accuracy step on the sine obstruction, not merely a smaller residual qualification around a matrix.

The original matrix-free resummed shift is executed unchanged:

    v=x/2+N/2, w=x/4+N/2+M/4,
    S=x-g(v)+g(g(w)).

The correction uses centered finite differences around S, not a producer derivative or an unshifted covariance expansion.

## 1. The missing nonlinear ancestry and its exact finite Gaussian bridge

On a stationary OU history, put

    H1=int_0^inf e^-s X_s ds,
    H2=int_0^inf s e^-s X_s ds.

Given X0=x, the pair (H1,H2) has exactly the law (v,w) above, with the SAME N in both rows. For s>=0 define

    r_s=e^-s, a_s=2s e^-s, b_s=2s(s-1)e^-s,
    c_s=sqrt(1-e^-2s[1+4s^2+4s^2(s-1)^2]),
    U_s=r_s x+a_s N+b_s M+c_s L,

where L is another fresh standard D-Gaussian. Then

    (X0,H1,H2,X_s) =_law (x,v,w,U_s)

when x is standard Gaussian independent of N,M,L. In particular this reproduces ancestry conditional on the already retained x,v,w, not merely the marginal law of X_s.

Indeed Cov(X_s,H1)=e^-s(s+1/2), Cov(X_s,H2)=e^-s(s^2/2+s/2+1/4). Subtracting the X0 regression gives s e^-s and (s^2+s)e^-s/2, respectively. The displayed a_s,b_s are the unique coefficients in the existing N,M rows. The remainder variance is positive for every s>0 and zero at s=0. The independent quadrature audit provides an exact rational certificate stronger than real positivity. No source-dependent covariance root or D-by-D matrix calculation is used; c_s is a known scalar.

Define the Gaussian operator

    T_s h(x,N,M)=E_L h(U_s).

For F1=int e^-s g(X_s)ds the exact conditional ancestry identity is

    E[F1 | X0,H1,H2] = int e^-s T_s g(x,N,M) ds.

This is precisely the nonlinear quantity that g(v) failed to reproduce.

## 2. A positive logarithmic-node bridge rule

There is a source-independent positive interior rule (p_j,s_j), j=1,...,J, with

    p_j>0, s_j>0, sum p_j=1, sum p_j e^-s_j=1/2,
    ||sum p_j T_s_j - int e^-s T_s ds||_(L2(gamma_D)->L2(gamma_3D))
         <= delta_br,
    J=O(log^2(1/delta_br)).

The lemma and exact arithmetic certificates are in bridge-quadrature-proof.md. Its mechanism is important. The row q(z)=e^-z(1,2z,2z(z-1)) has complex Euclidean norm at most one on every disk |z-s|<=s/64, s>0. Gaussian-chaos second quantization therefore gives a bounded holomorphic operator T_z on these disks, without analyticity of g. Exponentially weighted positive Gauss rules on geometric intervals, with small positive head/tail atoms, approximate the integral exponentially in the node count per interval. There are O(log(1/delta_br)) intervals. A positive convex admixture at e^-s=1/4 or 3/4 repairs the first exponential moment exactly, with only a constant-factor increase in error and one extra node.

Only the LINEAR operator T_s g is approximated in operator norm. No claim is made that the entire nonlinear correction is an operator of this form. Its Taylor remainder is bounded separately at each actual positive node in Section 4.

## 3. Literal finite source

Using the bridge rule, set

    D_j=g(v)-g(U_s_j),
    C_j(S,D_j)=[g(S+D_j)-g(S-D_j)]/2,
    F_br(x,N,M,L)=g(S)+sum_j p_j C_j(S,D_j),
    psi_br(x)=E_(N,M,L) F_br(x,N,M,L).

N,M are shared with the original v,w; L is shared across every bridge node. The old nonlinear response g(S) is retained. Signed vector arithmetic in the centered difference is not a signed probability measure.

The shared L does NOT make the collection (U_s_j)_j a genuine multi-time OU path. The proof needs only each one-time joint law with x,v,w and linearity of expectation, followed by nodewise remainder bounds. It never uses a false multi-time covariance identity.

With the pinned positive outer OU rule (omega_i,t_i), mass one and first moment 1/2, let c_i=sqrt(1-t_i^2) and use the SAME G,N,M,L across all outer nodes of one source occurrence:

    x_i=t_i Z+c_i G,
    F=sum_i omega_i F_br(x_i,N,M,L),
    B0_raw=sum_i omega_i g(x_i),
    E=F-B0_raw.

The exact own mean of this literal finite program is Q_out psi_br(Z). The notation B0_raw here names the random gradient baseline, not its caller origin; below it is denoted simply B.

At a single outer node there are four shared base VALUES: g(v),g(w),g(g(w)),g(S). Each bridge node adds three: g(U_j),g(S+D_j),g(S-D_j). Thus F costs (4+3J)N_out and E costs (5+3J)N_out before complete-key sharing. F/E own 4D raw private Gaussian coordinates; B owns D.

## 4. Actual bias theorem: shifted, ancestry-matched proof

Couple S to the genuine history by v=H1,w=H2. Write

    F1=int e^-t g(X_t)dt,
    F2=int e^-t g(X_t-F1,t)dt,
    B_shift=g(v)-g(g(w)), S=x-B_shift.

For any one block of dimension d, let sigma_d=[d(d+2)]^(1/4), the L4 norm of a standard d-Gaussian. Minkowski and Lipschitzness give

    ||F1||4<=A sigma_d,
    ||F2||4<=A(1+A) sigma_d,
    ||B_shift||4<=A[1/sqrt(2)+A k] sigma_d,
    ||B_shift-F2||4<=A u_A sigma_d,
    ||g(v)-g(U_s)||4<=A d0 sigma_d.

The last bound holds at EVERY s, since U_s is marginally standard. Taylor's theorem at the ACTUAL shifted S gives

    g(x-F2)=g(S)+Dg(S)(B_shift-F2)+R_true,
    ||R_true||2 <= (B A^3/2)u_A^2 sqrt(d(d+2)).

Each centered difference obeys, at its actual bridge node,

    C_j(S,D_j)=Dg(S)D_j+R_j,
    ||R_j||2 <= (B A^3/2)d0^2 sqrt(d(d+2)).

No derivative is executed: these are comparison identities only. With the exact bridge integral, the conditional true-minus-repaired linear discrepancy is exactly

    E[Dg(S)(F1-F2-g(g(w))) | x,v,w].

Its L2 norm is at most A^3(1+k)sqrt(d), because ||F1-F2||2<=A^2 sqrt(d). The bridge quadrature creates one additional global error at most A times delta_br ||g||2, hence delta_br A^2 sqrt(D).

Conditional Jensen, positive weights, sqrt(d(d+2))<=sqrt(b+2)sqrt(d), and summing squared block norms prove

    ||psi_br-psi2||_(L2 gamma)
       <=[C_(b,B,A)A^3+delta_br A^2]sqrt(D).

The outer R1 is an L2 contraction. Also ||F_br||2<=A C_F sqrt(D), so the pinned outer Hermite rule adds delta_out A C_F sqrt(D).

This proof does not invoke third derivatives, commute Hessians, or assume A sqrt(D) is small. The extra block curvature qualification is essential: replacing sqrt(d(d+2)) by a dimension-independent multiple of sqrt(d) for an unrestricted d=D block would be false.

### 4.1 Fixed-shape general-C2 qualitative improvement

There is also a weaker general-C2 conclusion without the block or curvature-Lipschitz qualification. Fix a dimension D and a single anchored C1 convex gradient h with 0<=Dh<=I, and set g_A=A h. With delta_br<=A and delta_out<=A^3,

    ||Q_out psi_br- m3||2 = o(A^2) as A -> 0.

This has NO uniform rate over h or dimension, and does not permit h to change with A. It is not a substitute for the quantitative block theorem.

For the bridge remainder let R=sqrt(|x|^2+|N|^2+|M|^2+|L|^2). The unit row norm gives |U_s|<=R for every s; for A<=1/2, S and S+/-D_s are within 3AR of x. Continuity of Dh implies that the centered Taylor remainder divided by A^2 tends to zero uniformly in s, pointwise in the roots. It is dominated by a fixed multiple of R in L2. Therefore dominated convergence applies even when all quadrature nodes depend on A, since their weights are positive with mass one.

For the genuine-history remainder, write J0=int e^-t|X_t|dt and J1=int t e^-t|X_t|dt. The displacement divided by A is bounded by |v|+|w|/2+J0+J1/2 for A<=1/2, an L2 random variable. Both endpoints of the shifted Taylor segment tend to x, so continuity of Dh and the same dominated-convergence argument give o(A^2). The remaining linear and bridge errors are O(A^3 sqrt(D)), and the outer error is O(A^4 sqrt(D)). Guarded own-mean completion adds its order-four allowance and restored floors.

## 5. Exact matrix quadratics and genuinely nonlinear examples

If g(x)=Kx, the centered difference is exactly K D_j. Consequently

    F_br=Kx-K^2 sum_j p_j U_j+K^3 w,
    E[F_br|x]=(K-K^2/2+K^3/4)x.

The outer first moment 1/2 gives the exact canonical m3. This holds for every 0<=K<=AI; its eigenbasis is never read. Define rbar=1/2, abar=sum p_j a_j, bbar=sum p_j b_j, cbar=sum p_j c_j, L_K=I-K/2+K^2/4 and beta=sum_i omega_i c_i. The actual shared-root raw output is

    F=(1/2)K L_K Z+beta K L_K G
      +(-abar K^2+K^3/2)N+(-bbar K^2+K^3/4)M-cbar K^2 L.

Its conditional covariance is the sum of the four displayed coefficient squares. In particular the beta^2 and cbar^2 terms retain cross-node root products. This is the raw own-source covariance, not a claimed target-history covariance.

The one-dimensional obstruction g_A(x)=A(x+sin x)/2 has b=1 and B=1/2. The new bias is O(A^3), whereas the old graph has the pinned lower bound 0.009194688258588945 A^2-3.04 A^3. Thus the new construction cancels the identified first nonlinear ancestry error on the actual obstructing family.

A genuinely noncommuting two-dimensional block family is

    g(z)=A[z/2+(sin(a.z)a+sin(b.z)b)/8],

where a,b are nonparallel nonorthogonal unit vectors. Its Hessian lies in [A/4,3A/4], Lip(Dg)<=A/4, and Hessians at different points generally do not commute. Arbitrary orthogonal sums and rotations of these blocks satisfy the theorem with block size two. No near-linear residual of order A^3 is assumed or needed.

## 6. Firsts, curl, caller origins, and literal zeros

At one outer node let V=Dg(v), J=Dg(g(w)), W=Dg(w), H=Dg(S), H0=Dg(x), H_j^U=Dg(U_j), and H_j^+=Dg(S+D_j), H_j^-=Dg(S-D_j). All original Hessians are symmetric in [0,AI]. Define

    P=I-V/2+JW/4, Q=-V/2+JW/2, R=JW/4,
    C_j=(H_j^+-H_j^-)/2, K_j=(H_j^++H_j^-)/2,
    Hbar=H+sum p_j C_j.

Then ||Hbar||<=3A/2, ||Hbar-H0||<=3A/2; both are symmetric. For y in {x,N,M,L},

    D_y E_node=Hbar S_y+sum p_j K_j (D_j)_y-H0 x_y,

    (S_x,S_N,S_M,S_L)=(P,Q,R,0),
    ((D_j)_x,(D_j)_N,(D_j)_M,(D_j)_L)
       =(V/2-r_j H_j^U, V/2-a_j H_j^U,
         -b_j H_j^U,-c_j H_j^U).

These are the exact noncommuting product rules. Put

    eta=A/2+A^2/4, q=A/2+A^2/2, r0=A^2/4,
    Cx=3/2+(3/2)eta+(3/2)A,
    Cn=(3/2)q+(3/2)A, Cm=(3/2)r0+A.

Since |r_j|,|a_j|,|b_j|,|c_j|<=1, the actual raw private first is at most

    A sqrt(beta^2 Cx^2+Cn^2+Cm^2+A^2),

and raw caller-Z first is at most A Cx/2. For the square lift P_G^*E=(E,0,0,0),

    curl <= beta[3A eta+3A^2]
            +A sqrt(Cn^2+Cm^2+A^2).

The leading G block Hbar-H0 is symmetric; only the explicitly small product terms contribute skew. No Hessian commutation or curvature qualification is used here. Under half-variance normalization, A<=1/2 and beta<=sqrt(3)/2 give private first <=4A and curl<=12A^2. Safe declarations are ell_E=4A,a_seed=3A,padding mu=A. The ell_E<=1/4 native guard therefore holds when A<=1/16; all OTHER native guards remain mandatory.

For every fixed finite p>=2, with ||N(0,I_D)||p<=kappa_p sqrt(D),

    ||E||_(Lp|Z=z)
       <=A^2[(3/4+A/8)|z|
          +kappa_p(sqrt(2)+A k+1)sqrt(D)].

The exact first exponential moment of the bridge rule was used here. The caller origin is the ACTUAL program at G=N=M=L=0; it has x_i=t_i z, v_i=t_i z/2, w_i=t_i z/4, U_ij=r_j t_i z, all subsequent descendants evaluated normally. It costs (5+3J)N_out before exact-key sharing and satisfies |E_origin(z)|<=A^2(3/4+A/8)|z|. Subtract and restore it. Anchoring leaves private first/curl unchanged and doubles the raw caller bound to A Cx. The baseline caller origin costs N_out and has magnitude <=A|z|/2, raw anchored caller first <=A. The half-variance-normalized compiler inputs therefore have caller bounds sqrt(2)A Cx and sqrt(2)A, respectively; native caller guards must use these actual normalized values.

At z=G=N=M=L=0 every descendant is literally zero; at nonzero z the captured origin is not declared zero. At g=0 the entire source vanishes for every root and the imported completed consumer retains its literal known Gaussian row. No normalization divides by a residual or vanishing energy.

## 7. Complete positive own-mean reentry and executable bill

Use the SAME pinned gradient/near-gradient finite compilers and coisometry adapter as the original matrix-free packet, with independent COMPLETE banks conditional on Z and every actual exterior label. The square lift now has dimension 4D. Use half variance shares, gradient order at least four, ell_E=4A,a_seed=3A,mu=A, and impose every native radius/curl/active-dimension/caller/finite-clock/filter/precision/mode guard at these actual dimensions and radii.

The imported near-gradient error is

    Lambda e_E[ell_E(a_seed+mu)+ell_E^3(1+mu^-1/2)]
       +absolute floors.

Its bracket is O(A^2) and anchored source energy is O(A^2)(|z|+sqrt(D)), yielding the O(Lambda A^4 sqrt(D)) integrated allowance. The gradient branch has its separate order-four allowance. Independent completed targets add to N(Q_out psi_br(Z),I); only then is product coupling used. No raw root or private source record is appended as an observer, and no completed output is read as a strong mean/covariance statistic.

All outputs are ordinary finite Gaussian pushforwards. The positive node weights, variance shares, fills, and buffers are actual executed rows. Centered finite-difference subtraction is deterministic arithmetic and creates no negative probability weight. Full compiler positivity is inherited only through the same guarded complete construction, not inferred from the scalar quadrature alone.

The safe original-VALUE bill is

    Q_rule/certificate_setup+Q_captured
      +N_out N_B+(5+3J)N_out N_E
      +Q_known/numerical/replay.

N_B,N_E must include every native replay, bank, filter, clock, pair, mark, finite-mode anchor and later occurrence at the ACTUAL normalized radii and dimension 4D. Captured origins and repeats are charged. Rules and their scalar square roots are source-independent; rigorous numerical construction and guard checking remain charged in setup/absolute floors. No matrix reconstruction, Hessian VALUE producer, expected-action leaf, or hidden continuous path is executed.

Each changed caller, input, node, coefficient or numerical version replays the whole graph. Each requested first/adjoint sweep uses one original HVP at every recorded VALUE site, with the above actual chains. A discarded primal pays complete replay; no HVP is differentiated. Exterior source parameters and original caller anchors remain live. The block decomposition, curvature certificate, and matrix witnesses are analytical only, not stopped fitted parameters in the graph.

## 8. Absolute numerical floors

At a node let errors in g(v),g(w),g(g(w)),g(S) be nu_v,nu_w,nu_b,nu_S. Let nu_Uj,nu_+j,nu_-j be errors at the actual bridge and terminal sites. A valid absolute F allowance is

    nu_S+sum_j p_j(nu_+j+nu_-j)/2
       +3A nu_v+2A nu_b+2A^2 nu_w+A sum_j p_j nu_Uj.

For a uniform original VALUE error nu this is (2+6A+2A^2)nu. E adds one baseline nu; anchoring doubles the resulting allowance before exact-key sharing. Root/node/weight/row/arithmetic/compiler/finite-mode/replay errors retain their absolute Lipschitz and native floors. Never divide a floor by measured energy or an A-dependent vanishing norm. Interior bridge nodes avoid an executed zero-variance endpoint square root; coherent exact-key reuse is needed for numerical literal zeros.

More explicitly, changing the bridge probability weights by Delta p gives an integrated raw floor at most A^2 d0 sqrt(D)||Delta p||_1. An absolute Euclidean error epsilon_j in the four scalar coefficients of row U_j gives an additional floor at most A^2 sqrt(D)sum_j |p_j'|epsilon_j, under the standard integrated Gaussian input. These are absolute row-coefficient tolerances, requiring no division by the small endpoint variance. Positive normalized weights and the first/curl bounds must be checked on the actual numerical coefficient version. A residual first-exponential-moment error epsilon_m contributes at most A^2 epsilon_m sqrt(D) to quadratic raw mean accuracy, rather than being silently treated as exactness.

## 9. Limitations and provenance

This proves the next nonlinear accuracy grade, not full order four. General C2 only gives executable ports and positivity; the displayed order-three bias additionally uses bounded block dimension and Lip(Dg_j)<=BA. Conditional normalization must re-certify those properties and its actual translated source; no uniform translated claim follows automatically. This does not solve endpoint joining, growing-order recurrence, strong conditioning on private roots, or arbitrary high-dimensional radial nonlinearities.

The repaired graph is not claimed to preserve the old five-Gaussian-residual order-four certificate: its additional shifted VALUE sites are new, and extrapolating the old certificate to them would be unjustified. A reusable service may retain the old graph on its original certified near-matrix class and use this repaired graph on the new block-curvature class, with the branch fixed by the declared structural certificate before execution and differentiation. That union preserves the old result without asserting that the new bound dominates it everywhere.

Pinned parent theorem SHA256: 772a0b758f7a4294ec9127629bec5bdd78b1f2ba167f32912ed5e7495e3ba06b.

Pinned parent manifest SHA256: 50af4b0dc59da95e5c656ffc5415515b61ff8882652f3de6538e3af46c8fbcb0.

The source packet is /workspace/shared/recovery-20261004/cost/endpoint-stein-quadrature-20261004/positive-law-repair/higher-cumulant-gate/same-carrier-feedback/matrix-free-polynomial-backbone. The parent finite compiler services are imported unchanged, with the new dimensions, counts and normalized bounds stated above; they are not reimplemented by the numerical diagnostics.
