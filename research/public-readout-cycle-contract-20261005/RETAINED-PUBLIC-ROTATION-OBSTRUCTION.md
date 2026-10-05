# A source-qualified obstruction to removing an exposed bank/public by Gaussian rotation

2026-10-05. This rules out a specified transform: preserve the actual exposed Gaussian coordinate and physical readout jointly, replace its dependency by a zero-readout fresh bank or an independently taped centered correction, and claim a higher-order joint error. It is not an impossibility theorem for a new joint native producer with explicit additional drift/covariance sources.

## 1. Exact covariance invariant for reparameterizations

Let B be the standard coefficient bank and X0=L B+R U the source-zero physical row, with B,U independent standard Gaussian blocks. The exact joint Gaussian covariance is Cov(B,X0)=L*. Any invertible Gaussian reparameterization of (B,U) that retains the SAME B and X0 as exported variables retains this covariance. Adding independent fresh coordinates does not change it.

If a proposed new independent coefficient bank Bnew has zero row in X0, replacing the original query ancestry B by Bnew changes the joint law unless all relevant query coefficients are independent of the replaced directions. Exact disintegration is instead

    B = L* C^(-1) X0 + S V,
    S S* = I-L* C^(-1)L,   C=Cov(X0),

with V independent of X0 (use the supported spectral subspace if singular). This formula preserves the law but leaves the ORIGINAL query center dependent on the public X0. It relocates the public derivative; it does not remove it. No unknown covariance inverse is required for the impossibility proof: L,C are the known Gaussian source-zero rows.

For the scalar standardized pair B=r X0+sqrt(1-r^2)V, a source coefficient f(B) becomes f(r X0+sqrt(1-r^2)V). Keeping the exact retained conditional density gives score proportional to r, with no alpha gain when r is O(1). If r=0, the canonical zero-row boundary is genuinely preserved, and this obstruction does not apply.

## 2. One actual C2 original source and a positive native channel

Take the one-dimensional strongly convex potential

    V(x)=x^2/2-eta cos x,   0<eta<1/2,
    g(x)=V'(x)=x+eta sin x.

Then V is C-infinity, 1-eta <= V'' <= 1+eta, so it lies within the original C2 bounded-Hessian source class. Choose a fixed shield t>0. At the actual captured center B, the scalar VALUE source

    f_B(x)=delta [g(B+t x)-g(B)-t x]/[(1+eta)t]

uses two original gradient VALUES and a known affine subtraction. It is a genuine scalar gradient, its selected first radius is at most delta eta/(1+eta), and its B-first is at most 2 delta eta/[(1+eta)t]. The ordinary native pair interface has ideal matrix

    M(B)=s0 E f'_B(Z)=epsilon cos B,
    epsilon=s0 delta eta exp(-t^2/2)/(1+eta).

All constants and source widths are fixed. At sufficiently small delta the covariance gap is positive. Its actual finite original-VALUE pair program approximates the following ideal kernel to any fixed higher power, with the same retained supplied input and center B:

    Q_epsilon=M(B) P+sqrt(1-M(B)^2) Z,

where B,P,Z are initially independent standard Gaussians. This is a legal positive source-qualified conditional kernel, not a stipulated inaccessible tensor.

There is also a version requiring NO affine subtraction in the source grammar. Use the literal anchored C0 source

    fraw_B(x)=delta [g(B+t x)-g(B)]/[(1+eta)t].

Its selected first is at most delta, its captured-center first is the same bound as above, and its native matrix is Mraw(B)=m0+epsilon cos B, with known m0=s0 delta/(1+eta). Compare to the known constant-m0 Gaussian pair rather than the zero-matrix pair. Every formula below is unchanged after replacing the baseline physical coefficient aP by (a+c m0)P and the baseline readout variance v by v+2ac m0. In particular the conditional mean DIFFERENCE remains c epsilon cos(B)P, and the same-grade tilted logarithmic DIFFERENCE remains exactly the displayed expression. Thus the obstruction does not depend on allowing any source adapter beyond the literal anchored original-g C0 one.

Let a,c,kappa>0 and retain the actual public readout

    W_epsilon=a P+c Q_epsilon+sqrt(kappa) G,

with independent standard G. At epsilon=0 its physical row is aP+cZ+sqrt(kappa)G and B has ZERO physical readout.

## 3. Retaining an internal public defeats an independent centered repair

The exact conditional cross moment is

    E[P W_epsilon | B]=a+c epsilon cos B.

Suppose a new correction T_epsilon is built on fresh dynamic tapes, conditionally independently of (P,Z,G) given B, and has any finite conditional mean m(B). This includes every independently taped centered permanent-leaf counterpacket and even noncentered packets whose mean depends only on B. Then

    E[P T_epsilon | B]=0.

Consequently W_epsilon+T_epsilon retains the SAME order-epsilon change in its joint (B,P)-cross moment. Adding an opposite-sign packet with its own new P' can cancel an unconditional readout covariance/current, but cannot cancel the old retained P-cross current.

There is an immediate conditional Wasserstein lower bound for a coupling that retains B and P exactly. The target independent-baseline port has conditional mean aP+m(B), whereas the corrected channel has mean aP+c epsilon cos(B)P+m(B). Jensen gives

    inf_{couplings retaining B,P} ||W_corrected-W_target||_2
        >= c |epsilon| sqrt(E cos^2 B)
        = c |epsilon| sqrt((1+exp(-2))/2).

The target can have arbitrary independent centered Gaussian reserve or independent correction noise; those do not change its conditional mean. Actual native errors o(epsilon) perturb the bound by o(epsilon), so the obstruction survives the real finite VALUE implementation. A correction explicitly depending on old P with an appropriate nonzero conditional mean could change the result; producing that mean is an additional joint source problem, not fresh-bank relabeling.

The discrepancy also gives an ordinary joint-W2 lower bound when the corrected output has the bounded fourth moments required of an admitted port. Use the explicit joint test F(B,P,W)=cos(B) P W. Its expectation difference is c epsilon E cos^2(B). For any coupling to (Bprime,Pprime,Wtarget), splitting the product and using |cos B-cos Bprime|<=|B-Bprime| yields

    |E F_actual-E F_target|
       <= [E(P^2 W^2)+E W^2+E(Pprime)^2]^(1/2)
          W2_joint(actual,target).

The bracket is O(1) under the finite fourth-moment guards, and E(Pprime)^2=1. Thus even allowing B and P to move within a general joint coupling leaves an order-|epsilon| obstruction. The explicit conditional lower bound above is sharper and matches the unchanged-input native contract directly.

This counterexample is already canonical ZERO-B-readout. Thus zero coefficient-bank readout alone does not authorize individual exposure of old internal physical publics. It specifically refutes the stronger recursive port discussed in the parent task.

## 4. Nonzero bank readout has same-grade physical-tilt descendants

Now let the output additionally contain beta B, beta nonzero and independent of epsilon:

    X_epsilon=beta B+a P+c Q_epsilon+sqrt(kappa)G.

Conditional on B, the remainder is centered Gaussian with exact variance

    v+2ac epsilon cos B,   v=a^2+c^2+kappa.

The characteristic function is exactly

    E exp(i theta X_epsilon)
      = exp(-v theta^2/2)
        E[ exp(i beta theta B) exp(-ac epsilon theta^2 cos B) ].

Differentiate at epsilon=0. Since E exp(iuB) cos B=exp(-(u^2+1)/2) cosh u, the first logarithmic coefficient is

    -ac epsilon theta^2 exp(-1/2) cosh(beta theta).

The same-grade rank-four coefficient is -ac epsilon exp(-1/2) beta^2 theta^4/2; all even higher ranks occur at the SAME epsilon power. Subtracting only the constant mean coefficient exp(-1/2) leaves those terms. No exact joint Gaussian reparameterization can alter this characteristic function. In the exact conditional-density rewrite B=(beta/Var X0)X0+S V the same coefficients reappear from differentiating the public-dependent conditional coefficient.

All variables, source centers and positive covariance completions were explicit. This therefore strengthens a generic polynomial observer warning to an original-C2-source-valid obstruction for the proposed rotation/averaging shortcut. It does not assert that the theta-series itself defines a positive counterkernel, nor that no new source family could cancel the displayed currents.

## 5. Connection to the sealed permanent-leaf packet

The lower-bound example is a simpler admitted native source, not a replacement for the sealed C2 eight-force tensor. The analogous issue in that packet is the mixed moment between its root output and its retained marked physical leaf inputs: the intended rank-four coefficient is precisely such a cross-channel correlation before the physical rows are consumed. An independent counterpacket on new leaf publics has no matching cross moment with the old leaf publics.

For a source-valid nonvanishing coefficient witness in the same original class, at a zero affine center

    E g'(sigma Z)=1+eta exp(-sigma^2/2),
    E g'''(t Z)=-eta exp(-t^2/2).

The C0 leaf and C2 center coefficients are both nonzero. Four such center factors and four positive leaf factors produce a nonzero literal scalar cubic-core coefficient with the original finite positive clocks and native normalization. No nonexistent constant cubic derivative on all of R is assumed. Exposing those marked inputs changes the joint current target. The grouped port avoids the issue only by declaring its own new physical marks internal before execution; it cannot retroactively consume an old marked input needed by the recursive caller.

## 6. Exact conclusion

Allowed: retain the complete external B and old publics, and insert a newly grouped zero kernel with fresh internal marked/public tapes into an unused reserve. All zero-target coefficient cancellation is conditional on the retained record.

Not justified: take an existing packet whose internal physical roots are still read later, replace it by its marginal grouped Gaussian law, and assert that a fresh bank or an orthogonal rotation preserved those old correlations. The exact native example has a strict order-epsilon retained-public lower bound. A wider recursive interface must expose the relevant roots before comparison and supply new joint drift/covariance/current producers for that larger target.
