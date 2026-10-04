# Nonlinear shear: conditional terminal heat and a scalar weak restoration

New candidate, 2026-10-04. Sections1-3 give an exact common-input decomposition
and a proved coefficient comparison. Sections4-6 specify an executable source
and the quantitative join to be checked against the positive-fork/old-core
interfaces. This is not a general nonlinear-DAG theorem.

## 1. Source and exact common-input curl

Let u,v be fixed orthogonal unit vectors in R^n, Pi=I-vv*, and a in(0,1/2).
Let g0:R->R have |g0'|<=1, and let g=grad V:R^n->R^n have ||Dg||<=1. These are
original VALUE primitives (or expressly admitted original-gradient charts),
with their actual calls counted. A literal single-potential instance is

    g0(z)=v* [g(x0+z v)-g(x0)]-lambda z,

where x0 and lambda are known fixed data and the declared bound |g0'|<=1 is
checked from the original Hessian sandwich. This is still a scalar gradient,
uses original g VALUES and known linear subtraction, and vanishes at zero.
The proof permits this shared-potential specialization; it does not assume
independence of the two primitive functions. Let

    q(W)=W+a u g0(v*W),
    f(W)=r[g(q(W))-g(W)],    kappa0=ra.

The second terminal cancels from curl. Write H=Dg, so

    Curl f(W)=kappa0 g0'(v*W)
          [H(q(W))u v* - v u*H(q(W))].                  (1)

The component of H(q)u parallel v cancels from the bracket. Since Pi u=u,
replace H in that bracket by Pi H Pi EXACTLY. This cancellation is important:
only the Hessian on the directions perpendicular to the ancestor is needed.
The whole f remains the marked source; no energy of its two terminals is
substituted. Its actual centered energy is denoted e. A declaration A>=3r
bounds its first, and a numerical multiple of kappa0 bounds its curl.

At a whole-input OU clock c^2+s^2=1, write

    W=cX+s(Z_perp+Z_v v),    z=c v*X+s Z_v,
    q_z(Z_perp)=c Pi X+s Z_perp+z v+a u g0(z).

Here Z_perp is standard on v-perp and Z_v is an independent scalar standard.
Define the symmetric terminal matrix on v-perp

    Hbar_X(z)=E_Zperp Pi H(q_z(Z_perp)) Pi,
    b_X(z)=Hbar_X(z)u.

Let J_f(X)=E_Z Df(cX+sZ), using a separate complete whole-f conditional
averaging bank. Then the exact conditional curl is

    C_f(X)=kappa0 E_Zv g0'(z)[b_X(z)v* - v b_X(z)*].       (2)

J_f depends on X, but does NOT depend on the newly exposed Z_v in (2).
The orientation node is O_X=-Sym(J_f(X) C_f(X)). Both copies preserve their
original complete f genealogy at the SAME X.

## 2. The conditional terminal matrix has a dimension-free z first

For every z, ||Hbar_X(z)||<=1. It is Lipschitz in z with

    ||Hbar_X(z')-Hbar_X(z)|| <= C |z'-z|/s.              (3)

No third derivative bound is used. For unit alpha,beta in v-perp, Gaussian
integration by parts in alpha gives

    alpha* Hbar_X(z) beta
      =s^(-1) E[(alpha*Z_perp) beta* g(q_z(Z_perp))].

The same Gaussian is used at z,z'. The query difference is
(z'-z)v+a u[g0(z')-g0(z)], of norm at most(1+a)|z'-z|. Since g is unit-first,
Cauchy-Schwarz bounds the difference by(1+a)|z'-z|/s, uniformly in the unit
test vectors. This proves (3), also for C2 potentials without Hessian modulus.
In particular ||b_X||<=1 and Lip_z b_X<=C/s, uniformly in X and dimension.

## 3. Weak smoothing of the scalar ancestor derivative

For eta>0 let m_eta(z)=E_x g0'(z+eta x), where x is scalar standard. At fixed X,
z has density p_mu,s of N(mu,s^2), mu=c v*X. For every unit vector t, set
B(z)=t*b_X(z). Then |B|<=1 and Lip B<=C/s. Translating the scalar density gives

 |E_z[(m_eta(z)-g0'(z))B(z)]|
 <=E_x integral |B(w-eta x)p_mu,s(w-eta x)-B(w)p_mu,s(w)|dw
 <=C eta/s.                                             (4)

Indeed the total variation of B p_mu,s is bounded by
Lip B+||B|| integral|p_mu,s'|<=C/s. The bound is uniform in mu; g0' need only
be bounded. The vector version follows by taking the supremum over t. Thus
replacing g0'(z) by m_eta(z) in(2) changes C_f(X) by operator norm at most
C kappa0 eta/s, uniformly in X.

Gaussian first-chaos Bessel on the WHOLE marked f gives

    ||J_f||_(L2_X;HS)<=e/s.

Consequently the orientation-node comparison, after Z_v is averaged but with
X still retained, has the one-energy bound

    ||O_X-O_X,eta||_(L2_X;HS)<=C kappa0 eta e/s^2.        (5)

The argument never bounds a pointwise derivative of g0 by a finite difference.
It smooths the derivative only against its actual correlated terminal consumer
and scalar Gaussian density. It also never differentiates the marked factor
in Z_v, since that factor is independent of Z_v conditional X.

For a positive orientation clock with weights w_j, choose

    eta_j=b s_j^2,   0<b<1.

Then sum w_j times (5) is at most C b kappa0 e, since sum w_j is numerical.
This is an energy-relative coefficient calibration for the exact common-input
query geometry. It is different from independently heating a frozen primitive
and asserting that its pointwise matrix is unchanged.

## 4. Literal gradient sources for the smoothed coefficient

At fixed(X,z) the terminal source is the genuine gradient on v-perp

    h_t(x)=(rho/s) Pi[g(q_z(0)+s Pi x)-g(q_z(0))].

Its selected matrix is rho Hbar_X(z) on that subspace. The scalar source

    h_0(y)=(rho/eta)[g0(z+eta y)-g0(z)]

is a genuine gradient on R and has selected matrix rho m_eta(z).
Both have active radius O(rho). Their old-label paths are respectively
O(rho/s) and O(rho/eta); the SAME z occurs in the terminal center and scalar
source. Each anchor is an actual original-gradient VALUE call and retains
its caller path.

Use the block source diag(h_t,h_0). The known skew contraction on
(v-perp) direct-sum R is

    J=[0 u; -u* 0].

Between the two complete gradient-pair banks use J and the known fill
sqrt(I-JJ*)Z_mid. This preserves standard input at the second pair; J need
not be orthogonal. The physical map L(x,y)=x+yv is an isometry from that
block space to R^n. The selected two-pair channel is a known numerical
multiple of

    rho^2 m_eta(z)[Hbar_X(z)u v* - v u*Hbar_X(z)].

The exact sign and c_sel factors are those of the general two-pair skew
adapter, with this fixed middle contraction and its additional known fill.
All dimensions/known maps can be represented using the supplied u,v rows;
no basis of unknown vectors or matrix-Hessian oracle is required.

The marked pair uses [f(cX+s w)-f(cX)]/s and is independent of Z_v conditional
X. Its complete private bank remains independent of the two side pairs,
conditional only on the shared passive and captured (X,Z_v). The selected
coefficient is exactly J_f(X). Thus the positive fork targets the smoothed
node O_X,eta after these owned roots are integrated.

## 5. Candidate positive-clock bill

Run a positive orientation-only clock at relative tolerance delta=c0 b, with
a fixed sufficiently small numerical c0 (and only the declared logarithmic
allocations). Its minimum width is then comparable to b up to those factors.
If a smaller unrelated delta is chosen, retain the actual s_min comparable to
delta and replace the last bill below by kappa0/(b^2 rho delta). At
each node use the complete fork readout sqrt(w_j), marked amplitude b and
side amplitude of order kappa0/(b rho^2). The marked deletion/law rows are
b e and b[kappa+A^(7/3)]e, up to their established positive sums.

The worst actual side-root path is rho/eta_j. Therefore its complete absolute
weighted root/caller bill is

    (kappa0/(b rho)) sum_j sqrt(w_j)/eta_j
      <=Lambda kappa0/(b^2 rho s_min)
      <=Lambda kappa0/(b^3 rho).                       (6)

The dyadic sum bound uses w_j=O(s_j^2), positive panel mass, and
sum sqrt(w_j)/s_j^2<=Lambda/s_min. It is an actual path bill, not the first of
an averaged coefficient. With kappa=A^(1+g), b=A^zeta and rho=A^gamma, the
ordinary A first/physical sqrt(A) caller guard has room if

    3 zeta+gamma<g.

Finite side priors and all eta_j floors must be priced through the real fork
readout; fixed positive zeta,gamma give finite fixed orders for each requested
nominal allowance. The variance shares remain feasible since the side
coefficient has positive exponent 1+g-zeta-2gamma under an appropriate choice.
There is no extra inverse-A replication, but all finite source/clock work and
higher compiler constants remain real.

## 6. Owned-root and old-endpoint gates to verify

Because m_eta'(z) is bounded by C/eta and the terminal derivative is O(1/s),
the complete smoothed O-node has a root-to-matrix-HS derivative bounded by
C A kappa(1/s+1/eta). After weights and eta=b s^2, this is Lambda A kappa/b.
The matrix derivative here remains dimension-safe: the skew field
b_X(z)v*-v b_X(z)* has rank at most two. A derivative of m_eta multiplies
J_f times that rank-two matrix, with HS norm at most C A, rather than a full
ambient HS norm of J_f. A derivative of J_f has the uniform input-to-HS
frame A/s; a derivative of the skew field has rank at most two and vector
norm C/s by the transverse heat cut. These give the displayed complete
A kappa(1/s+1/eta) frame. Its centered HS energy remains Lambda kappa e by
the w/s sum. The new-root
covariance Price current therefore has the expected one-energy bill
Lambda A kappa^2 e/b when the old reference is independent of these new roots;
with an old derivative frame A^2 there is the additional A^2 kappa e row.
This needs the exact conditional Gaussian reference/keep chronology of the
positive-fork consumer, rather than just the matrix norms.

In particular, the favorable s_min~b in(6) is available if the negative fork
is compared at fresh independent roots against an old constant covariance
reference containing +q O_f. If the old forward/Cov clock is regenerated at
these same finite nodes, its ordinary A e quadrature error instead requires
delta_clock<=b kappa/A and changes the minimum width and path bill. The
independent-constant-reference join is therefore a substantive premise to
check, not a free choice of labels.

The scalar weak restoration (5) is proved above. The complete positive source,
original-potential admission of g0, all theta profiles and retained host, and
the old constant-reference consumer are the remaining finite program checks
before claiming a nonlinear shear orientation repair.
