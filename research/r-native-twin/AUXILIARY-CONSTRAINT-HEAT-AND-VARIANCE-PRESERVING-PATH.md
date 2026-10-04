# Auxiliary constraint heat: exact effects and a variance-preserving path

New bounded source/current identities, 2026-10-04, for the six-VALUE ambient gradient lift64ed4369. The constructions below use only original VALUES and known Gaussian rows. Their covariance frames do not by themselves prove an ambient stationary law on the nonlinear constraint graph.

## 1. Common terminal heat and its two opposite constraint defects

At the original graph record R, keep S,Z,x0,x1,d=aDelta and lambda=0. Set p_i=g0(x_i)+sigma V with ONE new standard V for both i, independent of R. Then

    E_sigma=r[gt(t0+aM sigma V)-gt(t1+aM sigma V)],
    G_S=E_sigma,
    G_x0=-r sigma V,        G_x1=+r sigma V,
    G_lambda=0,
    G_p0+G_p1=aM*E_sigma+rM*d.                         (1)

The actual paired mark remains |E_sigma|<=r a^2 epsilon|Z| pointwise. This is a nominal envelope, not an automatic bound by the ACTUAL energy of unheated E.

The two defect coordinates are opposite known Gaussians, but they share V with E_sigma. Conditional Gaussian integration by parts gives the EXACT cross covariance

 Cov_V(E_sigma,G_x0)
    =-r^2 a sigma^2 E_V[Ht(t0+aM sigma V)-Ht(t1+aM sigma V)] M,

and the opposite sign for G_x1. Thus they cancel under an identical readout, not under arbitrary distinct observers. Declaring them independent Gaussian fills would discard a real coefficient.

On the original graph, the combined constraint contribution to a full differential pullback is

    r sigma(Dx1-Dx0)*V=r sigma Dd* M V.                (2)

Since ||Dd||<=C a, its conditional cross with E_sigma has a one-energy bound C ra sigma ||E_sigma-E_V E_sigma||_(L2V) by Bessel. This is the heated conditional mark. The derivative of the matrix Dd in(2) is not bounded by the original first-only contract.

For a fixed old record and a smooth test of an endpoint Y_sigma consisting of E_sigma, the defects with FIXED linear readouts, and an independent Gaussian keep, differentiation gives the ordinary width current E[DY_phi·partial_sigma Y]. If the only sigma dependence is the displayed common V scaling, then partial_sigma Y=sigma^-1(D_VY)V. Distributional Gaussian integration by parts yields

 d/dsigma E phi(Y_sigma)
    =sigma^-1 E[(Delta_VY_sigma)·Dphi(Y_sigma)
                    +(D_VY_sigma D_VY_sigma*):D^2phi(Y_sigma)].    (3)

The defect part has zero V-Laplacian, but the terminal part has the actual third-jet drift of gt. Formula(3) is a weak analytical identity, not an allowed higher-derivative VALUE oracle or a proved dimension-free one-mark norm bound. The opposite defects alone do not remove that drift. A changing nonlinear readout adds its own actual derivative terms.

## 2. Preserve the original Gaussian source law exactly

There is a useful alternative path which keeps the actual marked E law rather than asserting relative stability of terminal-only heat. Assume a sigma||M||<=1/2, and define the known matrices

    D=(I-a^2 sigma^2 M M*)^(1/2),
    D2=(I-a^2 sigma^2 M* M)^(1/2).

Set

    S'=D S+a sigma M V,
    V'=-a sigma M* S+D2 V.                            (4)

This is a known orthogonal Gaussian rotation because DM=MD2. Hence S',V' are independent standard, and both remain independent of U,Z. Known dense matrix/fill work is charged unless the native scalar-block geometry makes it cheap.

Rebuild the COMPLETE original paired graph at R'=(S',U,Z): x0'=cS'+sU, d'=a[gt(S')-gt(S'+epsilon Z)], x1'=x0'+M*d', and the corresponding t_i'. Do not reuse the old feedback at changed S'. Define the known Gaussian increment

    h=sigma V-a sigma^2 M*(I+D)^(-1)S.                (5)

No inverse of M is used. The identity (D-I)=-a^2 sigma^2 M M*(I+D)^(-1) gives aMh=S'-S. It is also useful that

    h=sigma V'+a sigma^2 M*(I+D)^(-1)S'.              (6)

In the AMBIENT source use S unchanged, x_i=x_i', p_i=g0(x_i')+h, ambient d=a[gt(S)-gt(S+epsilon Z)], and lambda=0. Its actual readout is now

    G_S=E(R'),
    G_x0=-rh,               G_x1=+rh,
    G_lambda=r M*(d'-d).                               (7)

The physical readout has EXACTLY the original E law and actual centered energy e, since R' is a standard original input. The price is the explicit new nonlinear constraint defect in G_lambda. The old zero blocks cannot all be claimed simultaneously along this law-preserving path.

## 3. Dimension-free one-energy covariance frames for the new defects

Because V' is independent of the complete R', equation(6) gives

 Cov(E(R'),G_x1)
    =ra sigma^2 E[(E(R')-EE) S'^*](I+D)^(-1)M,         (8)

and G_x0 has the opposite sign. Gaussian first-chaos Bessel bounds its Hilbert norm by C ra sigma^2 e. This is an actual one-energy statement about the preserved E mark, with no sqrt(d) loss.

For the lambda defect, fix the original Z and set H_Z(S)=gt(S)-gt(S+epsilon Z). Then d'-d=a[H_Z(S')-H_Z(S)]. The pair(S,S') has the reversible standard-Gaussian correlation D. For every fixed test u, the scalar function (Mu)*H_Z has Lipschitz constant at most2|u|. The Gaussian OU Dirichlet estimate gives

 E|u*G_lambda|^2
   <=8 r^2 a^2 ||I-D||op |u|^2
   <=C r^2 a^4 sigma^2 |u|^2.                         (9)

The conditional mean of G_lambda given Z is zero because S and S' have the same standard marginal. Thus (9) is a dimension-free ordinary row frame, and Hilbert/operator factorization yields

    ||Cov(E(R'),G_lambda)||HS <=C r a^2 sigma e.       (10)

This estimate is uniform in the captured Z before its supplied source profile is integrated. It uses variance-preserving common inputs, not a small Gaussian-norm displacement bound or independent primitive heat.

The actual unmarked Hilbert energies of the defect bodies can still include sqrt(d); equations(8)--(10) are mixed covariance frames with one intact E factor. Complete original-root firsts of G_lambda are O(ra), its fresh V first is O(ra^2 sigma), and its caller retains the literal four original-gt paths. A small covariance does not upgrade these actual firsts or its own Hilbert energy.

## 4. Precise remaining current/producer gate

The path is a finite VALUE construction of the original law together with explicitly priced defect fields. It does not make the full ambient feature vector Gaussian. The actual E and G_lambda share their original roots. To use their covariance frames inside a SAME-endpoint comparison with a nonlinear observer requires the complete corresponding Riesz/Price currents and proper cuts; multiplying by a root-dependent derivative or contracting a trace is not licensed by the fixed-test bounds above.

In particular the full graph pullback contains Dd and the paired original Hessian factors. Their next firsts are not controlled by differentiating a small error. A successful positive auxiliary-width generator must either realize those currents with original VALUES and an intact mark, or keep the associated larger body with a valid complete graph return. The present note supplies exact Gaussian-law preservation, exact constraint/cross identities and bounded covariance frames, not that final stationary producer.
