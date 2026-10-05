# Native ports and the old radial family for the new delayed graph

This accompanies DELAYED-CONDITIONAL-NOISE-THEOREM.md. The raw finite-source ports are displayed first. The pinned half-variance convention in nonlinear-bridge-repair/INDEPENDENT-REPAIRED-SOURCE-AUDIT.md, Section 3, divides the source by sqrt(1/2) without rescaling its inputs, so it multiplies both first and curl by sqrt(2).

## 1. Exact moment and derivatives

Impose the exact first moment sum_j a_j r_j=1/2 using the positive correction in Section 5 of the theorem. Put w=sqrt(A), q=1-w, sigma=sqrt(1-q^2), and

    b=V_Q=q sum_j a_j g(r_j(qx+sigma Z)+sqrt(1-r_j^2)N),
    u=x-b,
    E(x,Z,N)=g(u-wg(u))-g(x).

The same Z,N are used through every occurrence. Write H_f=Dg(u-wg(u)), H_i=Dg(u), H_0=Dg(x). Then

    E_x=(H_f-H_0)-w H_f H_i-H_f(I-wH_i)b_x.                 (1.1)

The first term is symmetric and has operator norm <=A. The remaining term has operator norm <=A^2 K, where

    ||b_x||<=q^2 A/2,
    K=w+q^2/2=(1+w^2)/2=(1+A)/2.                           (1.2)

This estimate preserves matrix factor order. The norm ||H_f(I-wH_i)||<=A follows from the contraction I-wH_i, without assuming that H_f and H_i commute.

For the complete private bank,

    ||b_(Z,N)|| <= q A sqrt(1-q^2/4),
    ||E_(Z,N)|| <= q A^2 sqrt(1-q^2/4).                     (1.3)

Indeed ||b_Z||<=q sigma A/2, while ||b_N||<=q A sum a_j sqrt(1-r_j^2)<=q A sqrt(3)/2. Combining the two blocks in Euclidean operator norm gives (1.3). This is a full-bank bound, independent of the number of quadrature nodes.

## 2. Outer roots and curl

Let the finite positive outer rule have exact first moment sum_l omega_l t_l=1/2. Introduce a D-dimensional Gaussian root G, shared at all outer nodes, set x_l=t_l z+sqrt(1-t_l^2)G, and run the whole preceding graph with common Z,N. Define

    beta_out=sum_l omega_l sqrt(1-t_l^2)<=sqrt(3)/2,
    E_raw=sum_l omega_l E(x_l,Z,N).

The raw complete private first and retained-z first obey

    private first <=sqrt[beta_out^2(A+A^2 K)^2
                           + A^4 q^2(1-q^2/4)],
    retained-z first <=(A+A^2 K)/2.                          (2.1)

Lift E_raw into the G output block with zero output in the Z,N blocks. The symmetric H_f-H_0 term causes no G-block curl. The remaining diagonal defect and off-diagonal blocks give

    raw curl <= 2 beta_out A^2 K + A^2 q sqrt(1-q^2/4).       (2.2)

These formulas, rather than a guessed normalization constant, are the exported ports. Uniformly for A<=1/12 they are bounded by 1.2A for raw private first and 2A^2 for raw curl. Under the pinned half-variance normalization, both quantities are multiplied by sqrt(2), so safe normalized ports are first<=2A and curl<=3A^2. This uses the actual source-only scaling and does not silently rescale Gaussian inputs.

Subtract and restore the actual caller origin with G=Z=N=0; preserve every descendant. Private first/curl remain the same and the retained-z first at most doubles. The original caller/scale/anchor dependencies remain live and retain their own recorded derivative/adjoint chains.

There are 3D private roots after the outer rule. A raw terminal occurrence has J+2 original VALUES; a residual occurrence adds the original baseline and therefore has J+3 VALUES before exact-key sharing. The full source bill includes every occurrence in the completed native program, all captured origins, full replays, independent complete banks, fills, modes, clocks, buffers and numerical versions. Each derivative/adjoint sweep acts on every recorded original VALUE. No compiler completion guarantee follows merely from these ports.

## 3. Parameter and precision guards

The delay is represented by q=1-sqrt(A); the executed graph never requires delta=-log q or an infinite OU trajectory. In the intended A<=1/12 regime, q>=1-1/sqrt(12)>0.71, and sigma=sqrt(2sqrt(A)-A). There is no near-zero denominator in any executed original-g input. The factor sigma^(-1) occurs only in the proof's weak-error allowance. Source derivatives hold A,q,sigma and all quadrature nodes fixed; any exterior dependency of A must retain its original scalar derivative chain.

Quadrature uses O(log^2(1/A)) nodes. Its known scalar setup, error certification and finite encoding must be paid. The exact-first-moment correction is positive and of size at most twice the initial operator tolerance. Stable evaluations of 1-r^2, sigma, normalization weights, and exact prescribed moments need their own arithmetic allowances. Neither tiny positive variances nor scalar moment discrepancies are silently clipped or erased.

To preserve the A^(11/4)sqrt(D) target allowance, operator quadrature tolerance is epsilon=O(A^(3/4)) internally and epsilon_out=O(A^(7/4)) externally. All original numerical VALUE floors stay absolute. The full compiler complexity also depends on actual caller/radius/filter/mode/tolerance and precision guards, and is not reduced to the displayed polylogarithmic source-node count.

## 4. Explicit old radial-family check

Use exactly the anchored C2 source from RADIAL-STAGE-TWO-OBSTRUCTION.md:

    g(x)=A sqrt(D) phi_(A,A^2)(|x|/sqrt(D)) x/|x|,
    phi(q)=(1/2)H_(A^2)(q-1/2)+(1/2)H_(A^2)(q-(1-A tau)),
    tau=(1-1/sqrt(2))/4.

Its Hessian interval is the required [0,A], so the theorem applies directly, including all descendants. The old unchanged-graph lower bound cannot be transferred to this new graph.

There is also an explicit leading-row check. Each original U_j is marginal standard. In the Gaussian radial row reduction, g(U_j)=A c_U U_j, where c_U=phi(1)=1/4+O(A). Since the rule has exact first moment, the averaged U row has x coefficient q/2. Thus b has leading x coefficient A q^2/8. The inner u=x-b and final u-wg(u) produce total first-order radial contraction

    c_1 [q^2/2+w] = (c_1/2)(1+w^2),     c_1=1/4.

With w^2=A this differs from the genuine history's leading contraction c_1/2 only by O(A). Therefore the old order-A^2 first-chaos mismatch is absent; this contribution is O(A^3). For sufficiently small A both the final graph row and genuine row lie below the upper hinge shell, so the stated leading expansion is valid. Positivity bounds the coherent Gaussian-row-reduction error by O(A) independently of J; the C2 hinge smoothing transfer contributes O(A^3sqrt(D)). In particular the old witness divided by A^2sqrt(D) tends to zero along A=1/n,D=n^6.

The executable check_delayed_radial.py uses the exact displayed smooth phi, the certified positive quadrature construction, exact scalar first moment, and 70-digit arithmetic for the final cancellation. It records quadrature masses, first moments, sampled spectral errors, theorem constants and the Gaussian-row first-chaos witness in delayed_radial_checks.json. The sampled spectral check is an arithmetic sanity check; the all-Hermite-degree operator certificate is the theorem's analytic bound. The radial row check is likewise an independent sanity check, not a substitute for the uniform theorem.

The raw radial witness divided by A^2 is approximately -4.20e-6 at A=1e-4 and -4.22e-8 at A=1e-6, consistent with an O(A^3) leading defect. This check makes no claim that the complete bias norm equals that witness.

## 5. Optional exact-quadratic correction, separate from the main graph

For g(x)=Kx, the delayed source's conditional mean has coefficients

    Kx - [(1+w^2)/2]K^2x + (w q^2/2)K^3x.

The genuine target before outer R1 is Kx-(1/2)K^2x+(1/4)K^3x. Therefore adding the actual source

    C(x)=(w^2/2)g(g(x))+(1/4-w q^2/2)g(g(g(x)))             (5.1)

restores exact anisotropic quadratic behavior. For w=sqrt(A), both coefficients are nonnegative in the intended range, and ||C||2 <=(3/4)A^3sqrt(D). It adds three original VALUES, or two if g(x) was already recorded as a baseline. All nested descendants in (5.1) must be executed and paid. The main theorem does not need this correction; using it adds its stated O(A^3sqrt(D)) general-source allowance and requires correspondingly updated first/curl/numerical ports.
