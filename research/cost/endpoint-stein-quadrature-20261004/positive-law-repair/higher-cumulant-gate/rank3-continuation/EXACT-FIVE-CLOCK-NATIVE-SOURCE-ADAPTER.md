# Exact five-clock native-source adapter for the conditional cubic tensor

2026-10-04. Literal specialization of the native three-marked path. This uses the exact compressed target and clock proof in the sibling order-reentry notes; it does not form a finite path cube. Feedback is handled by the separately pinned positive-rank3 feedback lemma, pending independent audit.

## 1. Five clocks and the original query records

The target is

    kappa_3(H|Z)=24 R3 Sym_3[B . R3((DB).B)],
    B=R2(Dg), DB=D R2(Dg),
    Rk h=integral_0^1 r^(k-1) P_r h dr.

The inner tensor is raw J_jka=sum_b partial_a B_jb B_kb. It must not be symmetrized across the derivative index a before contracting with the outside B. Final physical symmetrization is legal.

At a positive five-clock node (q,s,r1,r2,r3), with its complete known scalar weight w already containing q^2 s^2 r1 r2^2 r3 and all quadrature weights, draw

    Y=qZ+sqrt(1-q^2)G0,
    X=sY+sqrt(1-s^2)G1,
    Q1=r1Y+sqrt(1-r1^2)G2,
    Q2=r2X+sqrt(1-r2^2)G3,
    Q3=r3X+sqrt(1-r3^2)G4.

These are the exact original g query sites for the Dg--D^2g--Dg tree. G2,G3,G4 are independent; the Y,X ancestry remains shared.

Put ci=sqrt(1-ri^2), sigma_i=ci/sqrt(2), and split each ci Gi into a coarse sigma_i Wi and an independent sigma_i Zi. The coarse bank Q consists of G0,G1,W1,W2,W3. Its injections into the three actual query centers have norm bounded by a numerical constant. Conditional on this full coarse bank, the original forces have independent sigma_i shields. No min-clock artificial common floor is needed.

## 2. Literal original-gradient VALUE sources

At each captured query center qi define

    Fi(u)=[g(qi+sigma_i u)-g(qi)]/(A sigma_i).

Every Fi is a genuine D-dimensional gradient source with fresh radius<=1, zero at u=0, and coarse/retained-Z first<=C/sigma_i. The pointwise PSD Hessian sandwich is preserved up to the positive normalization. The caller origin g(qi) is a real original VALUE record, reused only at its identical complete key.

Run the native C0 selected source at each leaf with scalar normalization sqrt(2)/s0. Run the finite separate-first-chaos C1 selected source at the center with normalization sqrt(3)/s0, including its private OU/Sobolev preparation. Their analytical selected matrix targets are

    H1=E Dg(Q1)/A,
    H2(P)=(sigma_2/A) E D^2g(Q2)[P],
    H3=E Dg(Q3)/A.

The original D^2g expression means a derivative of the privately smoothed g coefficient; it is never a primal query. All actually evaluated objects are the original g VALUES in Fi and the finite native pair/filter graph.

Use all-A leaf and center amplitudes. For a packet with readout coefficients c,a,b and positive untouched share eta, use

    rho_1=4w A/(cab sigma_2),
    rho_2=rho_3=A,
    q_packet=rho_1 rho_2 rho_3=4w A^3/(cab sigma_2).

Then the native reverse symbol has leading logarithmic third coefficient

    cab q_packet * (sigma_2/A^3) tree =4w tree.

This is precisely kappa_3/6 at this node, because the compressed target has prefactor 24. The source sign chooses the desired sign of the cubic current; probabilities remain positive.

## 3. Positive shares, guards and actual first paths

For N total positive history/node packets, reserve fixed total variance v3 and assign nu=v3/N to each. One example is c=a=b=eta=sqrt(nu)/2. Allocate the remaining caller/mean/covariance/untouched variance separately before execution. All inverse c,a,b factors are known public-log numbers.

Check the admitted pair/source guards at

    Lambda A,
    Lambda w A/(nu^(3/2) sigma_2),
    q_packet * (bounded native field envelope) <=1/4.

Positive dyadic cell weights obey bounded w/sigma_2. Thus the numerical small-A/public-log window makes every actual source admissible. No tiny realized source energy is used in a denominator.

The root fresh/private radius is Lambda rho_1. A root coarse-Z path costs Lambda rho_1/sigma_1. A center coarse path costs Lambda rho_1 A/sigma_2; the final leaf costs Lambda rho_1 A^2/sigma_3. Known public readouts are part of the Gaussian carrier. Consequently a conservative complete residual/retained-Z first is bounded by

    Lambda A sum_nodes w [1/sigma_2
                      +1/(sigma_2 sigma_1)
                      +1/sigma_2^2
                      +1/(sigma_2 sigma_3)] <=Lambda A.

The sigma_2^-2 sum is only logarithmic, and the others are integrable positive products. Strict-ancestor factors make several terms smaller; this display deliberately keeps a simple unchanged-scale return.

A fixed pair order b>=4 puts the root prior-kernel floor at order Lambda A^4 sqrt(D), with descendant errors attenuated further. Use larger fixed b and corresponding finite filter precision if an arbitrarily small separate absolute floor is required. Source/clock versions, shares, gates and numerical tolerances are frozen before first or adjoint sweeps.

## 4. Feedback, one energy and finite work

For the normalized three-tensor A_Q of this node, the native heat-cut facts are

    proper cuts <=C,
    proper cuts of D_Q A_Q <=C(1/sigma_1+1/sigma_2+1/sigma_3),
    ||A_Q||_(L2;HS)<=C sqrt(D).

The positive-rank3 feedback lemma therefore gives node error

    Lambda w^2 A^6 sqrt(D)/sigma_2^2
               *[1+1/sigma_1+1/sigma_2+1/sigma_3],

besides the explicitly propagated finite calibration/prior/clock/numerical floors. All squared-weight sums are integrable or at worst public-log at the endpoint cutoffs. Thus its substantive sum is Lambda A^6 sqrt(D), below the desired A^4 sqrt(D) law allowance.

There is also an integrated actual-H-energy refinement when g has PSD Hessian. From C=2R2(B^2),

    ||e_H(Z)||_(L2(Z))=||B||_(L2(Z);HS).

The marked leaf before owned-coarse averaging is not P_r Dg: its half-shield split gives P_rprime Dg at a standard analytical root, where rprime=sqrt((1+r^2)/2). Therefore a PSD inequality applied directly to the B-clock terms would not suffice.

The positive dyadic B-clock has the dilation moment envelope

    sum_j w_j r_j (rprime_j)^m <= C/(m+2), m>=0.

To prove it, note that 1-rprime is comparable to 1-r on every endpoint dyadic panel, and sum the positive panel masses against exp(-c m(1-r)); the final midpoint has the same bound. Away from the endpoint the multiplier decays exponentially. For reference, the exact integral is

    integral_0^1 r ((1+r^2)/2)^(m/2) dr
       =2[1-2^(-(m+2)/2)]/(m+2).

Hermite isometry now gives

    ||sum_j w_j r_j P_rprime_j(Dg)||_(L2 HS)
        <=C ||R2(Dg)||_(L2 HS).

Every P_rprime_j(Dg) is PSD. Reindexing these fields onto one analytical standard Gaussian root (only for this energy calculation), their nonnegative HS cross inner products imply

    sum_j w_j r_j ||P_rprime_j(Dg)||_(L2 HS)
      <=sqrt(N_B) ||sum_j w_j r_j P_rprime_j(Dg)||_(L2 HS)
      <=C sqrt(N_B) ||B||_(L2 HS).

This is the actual half-shield marked energy, summed with its literal positive leaf-clock weights. It does not identify or replace the different executed coarse roots. Thus the integrated centered H energy is retained up to public logs; finite coefficient calibration keeps its separate absolute floor. This is one physical energy, not a product of two dimension-sized vector energies. The PSD/dilation argument is source-qualified and is not asserted for arbitrary signed gradient sources or arbitrary deterministic retained Z.

Each raw native source occurrence consists of original g VALUES at the displayed affine query and its same-center anchor. C0 uses at most two such Fi calls; C1 plus its finite first-chaos filter uses its actual logarithmic number of complete Fi calls. A safe total bill is

    O(N_q N_s N_r1 N_r2 N_r3
        * N_pair * N_filter)

complete original-source occurrences, with fixed three-vertex and anchor factors, plus captured caller/mode work and numerical operations. The complete Gaussian bank includes all original coarse roots, private shields, pair/filter replicas, marked publics, and Gaussian fills. Original HVPs occur only in requested first/adjoint sweeps at the recorded original-g VALUE sites. Discarded primal records require their complete replay.

## Scope

The exact five-clock target and its quadrature certification, this literal native adapter, and the positive-rank3 feedback lemma together specify a positive cubic-current packet with known carrier and conservative residual first Lambda A. The independent law-consumer comparison is in the positive-quadratic-reference note. All should be audited together before the packet is imported as closed.

This packet does not itself upgrade a merely cubic covariance service to fourth order and does not turn its completed output into a genuine gradient. Those remain separate completion and re-entry questions.
