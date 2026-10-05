# A finite conditional resummation of the genuine history pair covariance

2026-10-05. Bounded-block, general-C2 result. The source is the actual coherent history, not a local atom silently relabeled as history. This does not close unrestricted-D m4, the singleton/base history sources, or the marked whole-joint/native port.

## 1. Result and exact boundary

Let X be stationary standard OU, g=grad U, g(0)=0, 0<=Dg<=AI, 0<A<=1/2. If a bounded-block conclusion is claimed, assume g separates in a supplied orthogonal physical decomposition with block dimensions n_j<=b. The unrestricted choice is b=D, and the displayed losses must then be retained.

Use the actual coherent W=(F2,Delta3), Delta3=F3-F2, on one future. For T>=1, let S_k consist of Y=X0 and the endpoints of the N_k=2^k cells of length d_k=T/2^k in [0,T]. Given S_k, the new midpoint roots are independent. Write P_{k,2} for the exact Hoeffding projection onto subsets of exactly two midpoint blocks and define

    C_pair,T(Y) = sum_{k>=0} E[(P_{k,2} E[W|S_{k+1}])
                              (P_{k,2} E[W|S_{k+1}])* |Y].

This is the genuine positive 2D-by-2D pair part of the history block Gram, including its cross block. It is not the entire history covariance.

There is an explicit finite, positive, original-VALUE law approximating a Gaussian with this covariance and an allocated positive buffer. Its pair covariance target error is bounded, in L2(Y;HS), by

    C sqrt(D) [
       (1+b) A^(38/9) + b^(3/2) A^(40/9)
       + b A^(9/2)(1+log_+ T)
       + b^(3/2) A^5(1+log_+ T)
       + sqrt(b) A^5 K + b A^7 + b^(3/2) A^10 ],       (1)

where d_K is within a factor two of A^(10/9) and K=O(log(T/A)). For fixed b and logarithmic T the leading error is O_b(A^(38/9) sqrt(D)), up to the explicitly smaller logarithmic terms. The exponent 38/9 is strictly greater than four. The source-class assumption is only C2 and the Hessian interval: no modulus of continuity for Dg is used.

For the active F2 slot, the positive Gaussian realization has additional conditional integrated W2 error

    C A^8 b^(3/2) (1+L)^2 sqrt(D)/v^(3/2)
      + C A^4 sqrt(D)/sqrt(v) (1+sqrt(b/L)) exp(-L),     (2)

under the explicit guard R_j^2/v<=1/8 in every physical block (sufficiently, C A^4(b+L)/v<=1/8); fixed numerical constants can be reduced further. Here L>=1 is the cap logarithm. The history target covariance error (1) adds at most its value divided by 2sqrt(v). A fresh independent buffer supplies the inactive Delta slot; a [I I] readout is charged in the usual way. No mark is retained through this marginal LAW comparison.

At most fourteen original full-vector g VALUES are needed per complete law draw when the known physical blocks are executed in parallel inside those full-vector calls. The Gaussian root count is at most 17D for the active slot, plus the separately declared inactive block buffer and O(D) scalar clock roots. The shared clock-table setup is polylogarithmic; selection is polylogarithmic per physical block. Total vector and clock arithmetic is O(D polylog(1/A,D,1/epsilon)) when the supplied block coordinates are the oracle coordinates; otherwise add the cost of 28 applications of the supplied orthogonal basis transforms. Dense transforms can cost O(D^2) arithmetic. This is not dimension-free scalar work. Every draw includes all source ancestry it uses. There is no path grid, derivative oracle, conditional-mean oracle, conditional covariance oracle, or executed analytical field.

The existing full coherent extension, COHERENT-BLOCK-ORDER-TWO-EXTENSION.md, supplies the order-greater-than-two tail and adds at most C[A^(22/5)b^(6/5)+A^5b^(3/2)(1+log_+T)]sqrt(D). The cited nontrivial tail bound is stated under A^3 b^(3/2)<=1; in the complementary regime the trivial CA^2sqrt(D) cap is already bounded by its displayed leading expression. Therefore the same producer also approximates the ENTIRE interaction-order-at-least-two part, with that additional paid error. For fixed b it is smaller than the leading A^(38/9) term. Alternatively, apply the Delta3 omission argument below directly to P_{>=2}, which still annihilates F1, paying CA^5 K sqrt(bD) on retained levels and using the full coherent sensitivity cap below the cutoff. This independently covers the Delta3 higher-order cross/diagonal slots. This does not include singleton components.

This is a bounded-block covariance-component producer. It does NOT provide the coarse base covariance, the exact singleton history fields, a retained-Delta current law, cubic/quartic marked cuts, a bounded full-first native source, or an unrestricted-D grade improvement. Formula (1) with b=D does not give a dimension-uniform gain.

## 2. Uniform comparison to the actual nonlinear pair component

Work in one physical block of dimension n. Put R_g=F2-F1. Then

    ||R_g||_L4 <= C A^2 sqrt(n),
    ||Delta3||_L4 <= C A^3 sqrt(n).                     (3)

The existing uniform coherent mollification argument applies at each dyadic level. Define the anchored Gaussian mollification

    h_e(x)=E[g(x+eZ)-g(eZ)],
    B_h^T=integral_0^T e^-t Dh(X_t) F1,h,t dt.

Its source and every ancestor are changed coherently in the comparison. With

    J_d = C0(1+d^(-1/2)),  e_d=sqrt(A/J_d),

where C0 is a fixed sufficient bridge-variance constant, that proof gives an actual random comparison V_d satisfying

    P_{k,2} E[R_g+B_h^T|S_{k+1}] = P_{k,2} V_d,
    ||V_d||_L4 <= C A^(5/2) n sqrt(J_d).               (4)

The future t>=T has zero midpoint dependence given S_k. The future beyond T inside the genuine F1 ancestors is still present in (4); it was not replaced by a mean. The two inputs to the proof are conditional Gaussian weak transport for the entire ancestral shift and the explicitly paid bound Lip(Dh_e)<=CA/e.

Set X_d=P_{k,2}E[R_g|S_{k+1}] and E_d=P_{k,2}V_d. Orthogonality is used only in conditional L2. In particular

    E[|X_d|^2|S_k] <= E[|R_g|^2|S_k],
    E[|E_d|^2|S_k] <= E[|V_d|^2|S_k].

Average onto Y, apply conditional Cauchy-Schwarz to the cross covariance, and then ordinary Holder in Y. This gives

    ||Cov_pair(F2)-Cov_pair(-B_h^T)||_(L2(Y;HS))
      <= 2||R_g||_4 ||V_d||_4 + ||V_d||_4^2
      <= C[A^(9/2)n^(3/2)sqrt(J_d)+A^5 n^2 J_d].     (5)

There is no L4-contraction assertion for a high-order or exact-order Hoeffding projector. P_{k,2}F1=0 exactly, which is why the small A^2 residual, rather than the A-sized F2 field, appears in (5).

To retain the full coherent W target, approximate its pair block by (-P_{k,2}B_h^T,0). The omitted Delta3 pair cross term costs at most C A^5 n per level and its own covariance at most C A^6 n. This follows from (3) and the same conditional-L2 contraction, again pairing Delta3 with R_g because the F1 pair projection vanishes. It is a bounded-block energy argument, not a falsely small Delta3 first.

For a block-separable source, all these covariance fields are block diagonal. The inequalities sum n_j=D, sum n_j^3<=b^2D, sum n_j^4<=b^3D, and sum n_j^2<=bD give the b factors in (1).

At finer levels discard the TRUE pair terms. The existing complete-history sensitivity bound gives their positive sum the one-HS estimate

    C A^2 d_K^2 sqrt(D).

For d<=1, sum sqrt(J_d)<=C d_K^(-1/4) and sum J_d<=C d_K^(-1/2); the levels with d>1 cost O(1+log_+T). Choosing d_K comparable to A^(10/9) balances A^2 d_K^2 with A^(9/2)d_K^(-1/4), producing A^(38/9). This proves the first five terms of (1). It does not approximate a singleton by a linearized field.

## 3. The complete coarse environment disappears for this comparison field

For two ordered coarse cells I=[a,a+d], J=[c,c+d], a+d<=c, the exact smoothed ordered response is

    B_h = integral_{t<s} e^-s Dh(X_t)h(X_s) dt ds.

Its pair Hoeffding projection factors as

    (B_h)_{I,J}=A_I(xi) B_J(zeta),
    A_I=integral_I {E[Dh(X_t)|endpoints,xi]
                           -E[Dh(X_t)|endpoints]}dt,
    B_J=integral_J e^-s {E[h(X_s)|endpoints,zeta]
                           -E[h(X_s)|endpoints]}ds.              (6)

All other coarse cells and the entire external future integrate out exactly. Ordered supports force t in I and s in J. This is an exact Markov/locality statement for B_h, following the uniform actual-history comparison (5); it is not a locality assertion about the full nonlinear F2.

Let P_t be ordinary physical OU. For the canonical cell [0,d], define the later-cell field

    B_d(x)=E_{X_d, zeta|X0=x}[B_local B_local*],

where B_local uses weight e^-s inside that cell and is centered in the midpoint given both endpoints. Define the completely positive local sandwich

    (Gamma_d M)(x)=E_{X_d,xi|X0=x}[A_local M(X_d) A_local*].

Then the pair covariance after averaging the full coarse skeleton, retaining Y, is exactly

    C_{a,c}(Y)=e^(-2c) (P_a Gamma_d P_{c-a-d} B_d)(Y).  (7)

The matrix order is literal. Since ||A_local||op<=Ad, Gamma_d has L2(matrix-HS) operator norm at most A^2d^2. Invariance of the OU law and Jensen prove this bound. P_z is a contraction for Re z>=0. Thus the two clock positions enter through genuine Markov semigroups with fixed local sources; the generic relative-gap obstruction for arbitrary interacting physical fields does not apply to (7).

On an infinite equally spaced mesh, the exact sum is

    e^(-2d) R_d Gamma_d R_d B_d,
    R_d=(I-e^(-2d)P_d)^(-1).                           (8)

On the actual finite N-cell mesh use i>=0, r>=0, i+r<=N-2, where a=id and c=(i+r+1)d. The positive total weight is explicit:

    M_{N,d}=e^(-2d)[1-Nq^(N-1)+(N-1)q^N]/(1-q)^2,
    q=e^(-2d).                                        (9)

Equations (7)-(9) are target identities. Execution below uses Gaussian rows and discrete positive random clocks, not P_t or Gamma_d as oracle leaves.

## 4. Seven original VALUES for one conditional response sample

Fix a level d with N>=2 and one ordered cell pair. Generate their four coherent endpoints from Y by the actual OU transitions; shared endpoints are literally identified when the gap is zero. Generate independent standardized midpoint roots xi,xi',zeta,zeta'. They are shared by both replicas in Section 5.

Choose local times t uniformly on [0,d] and s with density e^-s/(1-e^-d). Conditional on each cell's two endpoints and its midpoint, generate the selected bridge point using its exact scalar Gaussian regression rows. Use the SAME residual root for x,x', and the SAME residual root for z,z':

    x=a_t+sigma_t xi+kappa_t N0,
    x'=a_t+sigma_t xi'+kappa_t N0,
    z=b_s+tau_s zeta+lambda_s M0,
    z'=b_s+tau_s zeta'+lambda_s M0.

Here sigma_t,tau_s<=sqrt(tanh(d/2)). The means a_t,b_s depend only on the corresponding two endpoints. Add two independent smoothing roots N,M and execute

    u0=g(e M),
    u=g(z+e M)-u0,  u'=g(z'+e M)-u0,
    Psi_delta = [g(x+e N-delta u)-g(x'+e N-delta u)
                 -g(x+e N-delta u')+g(x'+e N-delta u')]/(2delta).
                                                               (10)

Exactly seven VALUES are used. The full nonlinear shifted g is evaluated; no expectation has passed through g. All branches retain their saved same-bank ancestors and anchors.

For every fixed outer midpoint pair and every nuisance root/time/endpoint, exchange of zeta,zeta' gives conditional mean zero. Differentiating only these two inner midpoint roots proves

    ||D_(zeta,zeta') Psi_delta||op <= A^2 tau_s/sqrt(2).           (11)

This bound is independent of delta. The complete first can be as large as CA/delta and is NOT declared small.

Average (10) over the nuisance roots and local times, keeping the endpoints and all four midpoint roots. Its delta->0 limit equals MINUS the rectangular difference of the local product A_I(xi) [e^c B_J(zeta)], with the correct factor 1/2, divided by

    L_d=d(1-e^-d).

The absolute e^-c factor from (6) is supplied separately through the e^-2c geometry weights in its covariance. Its sign does not affect the Gram.

The outer smoothing yields Dh_e. The inner anchored difference yields h_e. Independence of the two smoothing roots produces their product. Uniform Taylor control is legitimate only after the outer smoothing; it gives an integrated L4 error at most C delta A^3 n/e for the unweighted sample. Since |u|<=A|z|, no growing smoothing-root moment is hidden here.

## 5. Positive sampling of every level and cell pair

For retained k<K with N_k>=2 define

    w_k=M_{N_k,d_k} L_{d_k}^2 tanh(d_k/2)/2,
    S=sum_k w_k, p_k=w_k/S.

If S=0 the target is zero and only the buffer is emitted. Otherwise S<=C uniformly in T,K: for small d, w_k<=Cd^3; for large d, w_k<=Cd^2e^(-2d).

Sample the level with positive probabilities p_k. Sample i,r with probability e^(-2(i+r+1)d)/M_{N,d} on i+r<=N-2. These choices, all four endpoints, and all four midpoint roots are SHARED between the two replicas. Within each replica resample local t,s, bridge residuals N0,M0, and smoothing roots N,M independently.

Scale (10) by

    a_k=L_d sqrt(M_{N,d}/p_k),
    delta_k=e_k A^2 sqrt(tanh(d_k/2)).                 (12)

The scaled sample has the uniform conditional-inner radius

    ell=A^2 sqrt(S)<=CA^2.

Its conditional nuisance mean targets the normalized rectangular pair field. Consequently the Gram of this conditional mean, averaged over the shared geometry/endpoints/midpoints, is exactly the finite positive sum of the smoothed pair Grams, up to the explicit secant bias.

The factor a_k delta_k/e_k is uniformly O(A^2). Hence the secant mean-field error is at most CA^5 n in integrated L4. The target mean-field has L4 norm at most CA^2sqrt(n) by (11). Its covariance error is therefore

    C[A^7 n^(3/2)+A^10 n^2],                          (13)

which gives the final two terms of (1). This is a paid, uniform general-C2 secant error, not an unquantified derivative limit.

## 6. Two-replica covariance resummation inside a positive law

Let Psi and Psi' denote the two scaled samples with the shared roots just specified. Put

    R=ell(sqrt(n)+sqrt(2L)),
    F=projection of Psi onto the radius-R ball,
    F'=projection of Psi' onto the same ball,
    H=(F F'^*+F' F*)/2,
    Z=sqrt(v)(I+H/(2v))G,                             (14)

where G is fresh. The product H G is evaluated as [F(F' dot G)+F'(F dot G)]/2 in O(n) arithmetic; no matrix is materialized. This is one positive pushforward law, even though the random matrix H need not be positive. No signed probability or negative Gaussian variance is used. Its mean covariance coefficient is positive:

    E[H|Y]=E_outer[(E_nuisance F)(E_nuisance F)*|Y]=C_cap(Y)>=0.

Its exact covariance is

    Cov(Z|Y)=vI+C_cap(Y)+E[H^2|Y]/(4v).               (15)

The independent audit proves, using the exact noncommuting Gaussian likelihood overlap and rank(H)<=2, that for h=R^2/v sufficiently small,

    W2(Law(Z|Y), N(0,vI+C_cap(Y))) <= C sqrt(v) h^2.  (16)

Unlike independent positive covariance increments or direct nuisance averaging, the leading coefficient in (15) is the conditional-mean Gram itself. The unwanted variance appears only quadratically in that small covariance coefficient. This removes the old A^4/R_replica nuisance floor without increasing replica count.

The norm concentration following from (11) is uniform in every fixed non-inner root. The independent audit bounds the cap covariance bias and its Gaussian-law cost by

    ell^2/sqrt(v) [1+sqrt(n/(2L))] e^-L.              (17)

Combining (16)-(17), then coupling independently across the supplied physical blocks, proves (2). For L comparable to log(1/A) this has an explicit logarithmic-squared factor. No literal uncapped quartic rate without that factor is claimed.

The comparison consumes geometry, endpoints, all midpoint/nuisance roots and G. Y and declared unread labels remain. It is a conditional MARGINAL law. It neither preserves an actual history mark nor makes the sampled H readable after replacement by its Gaussian target.

## 7. Literal finite execution and numerical floors

One active-block draw uses two copies of (10): fourteen original VALUES. Shared roots are at most four endpoint innovations and four midpoint roots per physical coordinate. Each replica owns four Gaussian nuisance roots. The keep contributes one. Thus at most 17n Gaussian scalar roots are required, plus scalar clock randomness. The inactive Delta slot, if emitted, has its own fresh n-dimensional buffer. Source block parallelization can pack all calls into fourteen original D-vector queries, provided the separating physical basis is supplied and used consistently. If it is a dense rotated basis, separately charge 14 forward and 14 inverse basis applications; the original-VALUE count is unchanged but arithmetic need not be O(D polylog).

For frozen clock labels the source in (10), including all endpoint and nuisance rows, has a valid complete/caller first bound C a_k(A/delta_k+A^2). With the chosen parameters this is at most C A^(-3/2)d_K^(-5/4), up to harmless coarse-scale constants. Radial cap projection is nonexpansive. Differentiating the transport then gives a source/caller first profile C R L_source |G|/sqrt(v), while its G-first is at most sqrt(v)(1+R^2/(2v)). These are honest moment profiles, not globally bounded native firsts; the presence of |G| alone prevents that promotion. Changed Y, endpoints, geometry, smoothing scale or secant scale replays every affected original ancestor and anchor.

There are K=O(log(T/A)) positive level weights. For cell-pair selection, s=i+r has mass proportional to (s+1)q^s, 0<=s<=N-2. Its partial sums have the same closed form as (9). Inversion by binary search takes O(log N) arithmetic operations, followed by one uniform i in {0,...,s}. It does not enumerate N or N^2 cells. The local s-time is obtained by the explicit inverse exponential CDF; t is uniform. All weights are positive and finite-clock discretization is paid.

For precision, cap projection is nonexpansive, ||H||op<=R^2, and the map from H to Z has conditional L2 amplification at most sqrt(n)/(2sqrt(v)) in operator norm. Perturbing a source vector by epsilon therefore costs at most CR sqrt(n)epsilon/sqrt(v), before any later readout. Each original VALUE error enters the final scaled (10) with an explicit factor at most C a_k/delta_k; the anchored inner errors have the additional original A factor. At each fixed cell length, Gaussian regression rows are one-half Holder in physical time, with explicit constants bounded by a fixed polynomial in 1/d_K on the retained range; those scale factors must be included in the precision allocation. Cap the Gaussian root magnitudes for numerical analysis only and pay their Gaussian tail; then choose time precision, root precision, weight total-variation error, and VALUE floors so that these finitely many displayed amplifications sum to the allocated absolute error. Since delta_k, e_k and d_K are fixed powers of A and all tail/clock caps are logarithmic, the required bit precisions and per-block clock arithmetic are polylogarithmic in 1/A, D and the absolute floor. Total arithmetic includes the explicit O(D) vector/root work and the number of physical blocks. No numerical floor absorbs the substantive errors (1), (2), or (13).

This is a mathematical finite sampler specification, not an assertion that an implementation was executed against an unspecified g oracle. It does not invoke or certify a native compiler. Products with G, cap nonsmoothness and 1/delta complete firsts must not be promoted to a bounded-first source contract.

## 8. Remaining ports, stated minimally

The construction resolves the full coarse-environment supplier for the grade-four DYADIC PAIR component at fixed physical block size, including a paid general-C2 approximation of actual history and a finite positive marginal realization. It goes beyond a standalone local secant atom.

To obtain the complete next-history block, a new supplier is still needed for the base conditional covariance and singleton coefficients. Those retain the whole nonlinear ancestor and are not recovered by (6). General-D closure additionally needs a dimension-safe replacement for the bounded-block conditional transport/remainder estimates and the capped rank-two law bounds. The native marked response, higher marked cuts, same-carrier cubic/quartic/Wick feedback, prefix/endpoint/completion rows, and a whole-joint mark-preserving comparison remain separate.

The minimal new singleton port is a finite original-VALUE conditional field for E[g(X_t-F1,t)|retained endpoint/midpoint labels] (or a proved different whole-law representation), keeping the full ancestor nonlinearity, with the required covariance precision and qualified full first/cuts. The present proof never moves that expectation through g and never claims it is supplied.
