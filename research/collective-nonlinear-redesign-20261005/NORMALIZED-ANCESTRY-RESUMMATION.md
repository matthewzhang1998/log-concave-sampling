# A cheap normalized-ancestry escape, and its general-class boundary

2026-10-05. New original-VALUE graph. This is not a Taylor correction to the sealed stage-two graph, not an any-order general-class theorem, and not a modification of any sealed packet.

## Result

Let g=grad U, U in C2, g(0)=0 and 0<=Dg<=A I, 0<A<=1/2. Retain the genuine Gaussian linear histories

    v=x/2+N/2, w=x/4+N/2+M/4,
    r=1/sqrt(2), k=sqrt(3/8).

Execute the three original-VALUE graph

    h=k g(w/k),
    Ftilde=r g((v-h)/r),
    Ttilde=g(x-Ftilde).                                  (N)

Its own mean is psi_tilde(x)=E[Ttilde|x]. The new graph retains the full terminal nonlinear response. It queries the inner forces at variance-normalized radii rather than applying g directly to compressed Gaussian histories.

For every radial convex A-Lipschitz gradient on the FULL D-dimensional space,

    ||R1 psi_tilde - m3||_(L2 gamma_D)
       <= 4 A^2 + (A^3/sqrt(2)) sqrt(D).                 (1)

All anisotropic quadratics have exactly the genuine canonical-m3 mean. There is no radial oracle, radial basis input, Hessian, matrix reconstruction, or expectation query in (N). The same code runs on every member of the unrestricted class.

For the two-shell obstruction sequence A=1/n, D=(3000 n)^2, its bias divided by A^2 sqrt(D) is at most A/sqrt(2)+4A/3000, hence tends to zero. Thus the prior radial obstruction is escaped at THREE terminal original VALUES rather than at inverse-A replication cost.

The graph does NOT improve the uniform class. The explicitly smooth one-dimensional source g_A(x)=A(x+sin x)/2 has

    ||R1 psi_tilde - m3||_2 >= A^2/30000

for 0<A<=1/100000. Section 5 proves this with a first-chaos witness. Thus normalization fixes a real collective radial pathology, but a general-class port must retain additional history information. The D=1 counterexample itself is radial.

## 1. Genuine target and a radial linearization lemma

Let X_t be the stationary OU process with Cov(X_s,X_t)=exp(-|s-t|)I. Conditional on X_0=x, its genuine exponentially weighted first and second linear histories are v and w above. Write

    H1,a=integral exp(-h) X_(a+h) dh,
    F1,a=integral exp(-h) g(X_(a+h)) dh,
    F2=integral exp(-a) g(X_a-F1,a) da,
    psi(x)=E[g(x-F2)|x], m3=R1 psi.

These are analytical targets, never executing leaves. Also

    integral exp(-a)H1,a da=w.

For radial g(y)=h(|y|)y/|y|, let

    lambda(q)=h(q sqrt(D))/(q sqrt(D)), q>0.

Then 0<=lambda(q)<=A. If Y is a centered Gaussian with covariance q^2 I,

    ||g(Y)-lambda(q)Y||_2 <= sqrt(2) A q.              (2)

Proof: h'(s) and lambda(q) lie in [0,A], so s -> h(s)-lambda(q)s is A-Lipschitz and zero at q sqrt(D). Also E(||G|-sqrt(D))^2<=2. This applies to arbitrary Gaussian correlations with every retained row. No independent-time replacement is made.

For p,q>0,

    |lambda(p)-lambda(q)| <= A |p-q|/min(p,q).         (3)

This follows from |s h'(s)-h(s)|<=A s. Only the original Hessian interval is used.

## 2. Proof of the radial escape

Put lambda=lambda(1),

    q=sqrt(1-lambda+lambda^2/2),
    Y0=v-lambda w,
    s=sqrt(1/2-3lambda/4+3lambda^2/8),
    qhat=s/r=sqrt(1-3lambda/2+3lambda^2/4).

All are scalar analysis parameters. In particular q<=1, s<=r, qhat>=sqrt(7)/4>1/2 for lambda<=1/2.

By (2), Minkowski, and Lipschitz propagation,

    ||F1,a-lambda H1,a||_2 <=sqrt(2)A,
    ||F2-lambda(q)Y0||_2 <=sqrt(2)A(q+A).             (4)

For the executed graph,

    ||h-lambda w||_2 <=k sqrt(2)A,
    ||Ftilde-lambda(qhat)Y0||_2
       <=sqrt(2)A(s+kA).                              (5)

The first term in (5) is the propagated change of Y=v-h; the second applies (2) to Y0/r, whose covariance is qhat^2 I. Thus no purported Gaussian identity is applied to the nonlinear actual Y.

Since q^2-qhat^2=lambda/2-lambda^2/4 and q+qhat>1,

    0<=q-qhat<=lambda/2,
    |lambda(q)-lambda(qhat)|<=A lambda<=A^2.

Combining (4),(5) and ||Y0||_2=s sqrt(D) gives

    ||F2-Ftilde||_2
       <=sqrt(2)A[q+s+A(1+k)]+A^2 s sqrt(D).

The entire terminal g is retained, so multiply this by its A-Lipschitz constant. Conditional Jensen and the L2 contraction of R1 give (1), since sqrt(2)[1+r+(1+k)/2]<4.

This proof applies directly to the genuine smooth two-shell sources. There is no hinge approximation or neglected descendant.

## 3. Optional fourth probe: dimension-free radial error

For the SAME arbitrary full-dimensional source, fix a public unit vector e and query

    lambda_* = <e,g(sqrt(D)e)>/sqrt(D) in [0,A].

Set q_*^2=1-lambda_*+lambda_*^2/2,
 s_*^2=1/2-3lambda_*/4+3lambda_*^2/8, rho=s_*/q_*.
Replace only Ftilde in (N) by

    Ftilde_*=rho g((v-h)/rho).

For radial g, lambda_*=lambda(1) exactly, and the comparison coefficient in BOTH (4) and (5) is lambda(q_*). The mismatch term A^3 sqrt(D)/sqrt(2) disappears:

    ||R1 E[g(x-Ftilde_*)|x]-m3||_2 <=4A^2.             (6)

It remains exact on every anisotropic quadratic because rescaling a linear map cancels. The probe does not certify radiality or improve the unrestricted class. Its source/caller dependence must remain live when the source itself varies. The frozen-source private ports in Section 4 extend by replacing r with rho; arbitrary exterior caller ports are NOT claimed here. This optional variant is therefore not used to import a completed service without a fresh exterior-caller audit.

## 4. Three-VALUE graph ports, readset, and positive-law interface

Write H1=Dg(w/k), H2=Dg((v-h)/r), H3=Dg(x-Ftilde), H0=Dg(x). These are analysis labels only. Exact firsts are

    (Ttilde)_y = H3[x_y-H2 v_y+H2 H1 w_y].            (7)

Factor order is preserved. No Hessian is queried to obtain a VALUE. A requested first/adjoint sweep uses one original HVP at each of the three retained source sites; no HVP is differentiated.

Let E=Ttilde-g(x). Then

    ||E_x|| <=A(1+A/2+A^2/4),
    ||E_N|| <=(A^2/2)(1+A),
    ||E_M|| <=A^3/4,                                  (8)
    |E| <=A^2(|v|+A|w|).                              (9)

The leading E_x term H3-H0 is symmetric and has norm <=A. The remaining ordered products supply the curl, not a commutation assumption.

For any finite positive outer rule with t_i in [0,1], sum omega_i=1, and sum omega_i t_i=1/2, set x_i=t_i Z+sqrt(1-t_i^2)G and reuse the SAME G,N,M at every outer node. Let E_Q=sum omega_i E(x_i,N,M), beta=sum omega_i sqrt(1-t_i^2)<=sqrt(3)/2. Its raw private first is at most

    A sqrt[beta^2(1+A/2+A^2/4)^2
               +(A^2/4)(1+A)^2+A^4/16] <1.2 A.

For the square lift (E_Q,0,0) on (G,N,M), the raw skew/curl norm is at most

    A^2 [beta(1+A/2)+(1/2)sqrt((1+A)^2+A^2/4)]<1.85A^2.

The half-variance normalized bounds first<=2A and curl<=3A^2 are therefore safe. Raw retained-Z first is <=A(1+A/2+A^2/4)/2; subtracting/restoring the actual caller origin at G=N=M=0 doubles that caller bound. Exterior source/anchor/scale ports require their own recorded chains and guards.

The conditional residual energy satisfies, with kappa_p=||N(0,I_D)||_p/sqrt(D),

    ||E_Q||_(Lp|Z=z)
      <= A^2[(1/4+A/8)|z|+kappa_p(r+A k)sqrt(D)].       (10)

The captured origin is bounded by the first |z| term in (10). It is not declared zero at a nonzero caller. At all roots and caller zero, every original VALUE site is literally zero. At g=0 the full source vanishes.

The original-VALUE bill per replayed residual occurrence is 4 N_out: the three sites of (N) plus g(x_i). The baseline costs N_out. Private Gaussian dimension is 3D. The complete positive-own-mean interface, IF all pinned native compiler guards are verified at these new ports, has the explicit bill

    Q_rule + Q_captured + N_out N_B + 4 N_out N_E
       + Q_known/numerical/replay.

N_B,N_E are complete expanded native occurrence counts at these new radii and dimension 3D, including all captures, later replays, banks, clocks, marks, fills, modes and precision. Independent COMPLETE baseline/residual banks and positive variance shares are mandatory. Using the existing near-gradient bracket with ell_E=2A, a_seed=2A, padding mu=A gives the bracket ell_E(a_seed+mu)+ell_E^3(1+mu^(-1/2))=6A^2+8A^3+8A^(5/2), hence an O(A^4 sqrt(D)) residual completion allowance. This is an interface substitution, not an independently re-certified or executed compiler. A<=1/8 only supplies ell_E<=1/4; all other native guards remain mandatory. Baseline gradient order and outer quadrature tolerance must also meet the chosen error target.

Positivity is elementary for the raw graph (finite Gaussian pushforward), but that raw law is NOT N(its mean,I). The imported completed law is positive only through the guarded full compiler and independent coisometry. We do not call the raw graph a completed mean sampler.

With uniform original-VALUE error nu, the raw terminal Ttilde error is at most

    (1+A r+A^2 k)nu;

E adds one nu, and a separately executed caller capture adds its own full allowance. The bound propagates every changed descendant. Positive outer weights add no node-count factor. Scalar coefficient, weight/moment, Gaussian-row, arithmetic, original-HVP, replay and native numerical floors remain absolute, separately priced errors.

## 5. Exact source-valid obstruction to a uniform general-class claim

Take D=1 and f(x)=(x+sin x)/2, g=A f. This is anchored C-infinity, convex-gradient, and 0<=g'<=A. Define

    C(u)=E[X f'(X)f(Y)]
        =[u(1+exp(-1/2))+exp(-1)(u cosh(u)-sinh(u))]/4,

where X,Y are standard Gaussians of correlation u. Gaussian integration by parts verifies the displayed closed form.

The normalized graph's order-A^2 first-chaos term uses r C(r), r=1/sqrt(2); the genuine target uses integral_0^1 C(u)du. Their difference is

    d=(exp(-1)/4){r[r cosh(r)-sinh(r)]
                              -[sinh(1)-2 cosh(1)+2]}
      =-0.000134719878344... .                         (11)

It is strictly negative without numerical reliance. Its bracket has the convergent series

    sum_(j>=1) [2j/(2j+1)!]
                  [2^(-j-1)-1/(2j+2)].

The j=1 term vanishes; every j>=2 term is negative, and j=2 is -1/720. Hence |d|>=exp(-1)/2880.

For completeness the order remainder is explicit. Because |f|<=|x|, |f'|<=1 and |f''|<=1/2, the genuine terminal first-chaos remainder is bounded by

    [1+(sqrt(3)/4)(1+A)^2] A^3,

and the normalized graph's by

    [k+(sqrt(3)/4)(r+A k)^2] A^3.

Their sum is <4A^3 for A<=1/2. This uses L4 Gaussian norms and the genuine history integrals, not pointwise bounded random forces. Outer R1 multiplies first chaos by exactly 1/2. Therefore

    ||R1 psi_tilde-m3||_2
      >= exp(-1) A^2/5760-2A^3 >=A^2/30000

for A<=1/100000. The last rational comparison can use exp(1)<11/4.

Tensorizing this same scalar sine source coordinatewise gives an unrestricted-class lower bound A^2 sqrt(D)/30000 for the three-VALUE graph. It is not a full-D radial source, but it confirms the dimension-uniform failure on the general class.

The same leading obstruction applies to the optional fourth-probe variant: as A->0, rho=r+O(A), and the change of rho enters Ftilde only at order A^2, hence the terminal only at order A^3.

## 6. Meaning and remaining route

Variance normalization is a real cheap escape for the high-dimensional radial pathology, and it retains exact quadratic ancestry. It does not keep enough history information for arbitrary sources, as (11) proves. This isolates a boundary for these two specific normalized graphs: a general-class upgrade must alter their retained information or conditional law. No impossibility claim for all finite variance-normalized nestings is made.

A separate true-joint-history stratification construction is being investigated in this packet. It retains the entire empirical collective force before terminal evaluation, unlike (N), and can trade explicit path/readset growth for general-class accuracy. No claim of zero inverse-A exponent is made for that construction.
