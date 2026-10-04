# Two exact joint-gradient lifts of the whole K3 composite

New bounded source analysis, 2026-10-04. This tests the actual paired nonlinear queries from c8b20bf4 and the unresolved matrix interface36d25b9e. The first lift fails the admitted VALUE/first/one-energy contract. The second is an executable ambient genuine-gradient source with small first, but its exact K3 readout lives on a nonlinear constrained input law. Neither failure is an impossibility theorem for another producer.

Throughout, M is a known contraction, r,a,epsilon are fixed, and g0,gt are anchored original gradients of first at most one. Put

    Delta=gt(S)-gt(S+epsilon Z),
    x1=x+aM*Delta,
    t0=S+aM g0(x), t1=S+aM g0(x1),
    E=r[gt(t0)-gt(t1)].

The actual source has x=cS+sU with independent standard S,U,Z. For differentiation of the proposed lifts, first treat S,x,Z as independent coordinates; the known linear pullback x=cS+sU is stated separately.

## 1. Natural paired-potential lift: exact extra channels

Consider the analytical scalar potential

    Psi(S,x,Z)=r[Vt(t0)-Vt(t1)].                       (1)

Write H0=Dg0, Ht=Dgt, D=Ht(S)-Ht(S+epsilon Z), and Hplus=Ht(S+epsilon Z). Direct chain differentiation gives

    G_S=E-r a^2 D M H0(x1) M*gt(t1),
    G_x=r a[H0(x)M*gt(t0)-H0(x1)M*gt(t1)],
    G_Z=r a^2 epsilon Hplus M H0(x1)M*gt(t1).           (2)

These are the full gradient of (1); no companion can be discarded while retaining that genuine-gradient claim. Under the actual linear coordinate map x=cS+sU, the gradient blocks become (G_S+cG_x,sG_x,G_Z). The known readout first block minus c/s times the second recovers G_S, still with its explicit feedback correction.

The ancestor companion splits exactly as

    G_x=ra[H0(x)-H0(x1)]M*gt(t0)
                      +a H0(x1)M*E.                  (3)

Only the second term has the intact E mark. The first is a Hessian difference times an unmarked terminal value. Bounded Hessians do not give it an extra epsilon or the actual E energy. Likewise the corrections in G_S,G_Z carry gt(t1), rather than E.

The gradient in (2) would require original HVPs as executed VALUES, beyond the present oracle model. Worse, differentiating G_x includes

    ra[(D H0(x))[h]M*gt(t0)
          -(D H0(x1))[h]M*gt(t1)]                     (4)

when S,Z are fixed. Such third jets have no uniform bound under the original C2 potential assumption. For merely C2 potentials they need not even exist classically. The next smooth same-potential fixture proves that pairing does not cancel this defect.

## 2. Exact same-potential unbounded-first fixture

Use d=1, scalar M=1, and the SAME smooth strongly convex primitive at both nodes:

    g(y)=m y+(d0/k)sin(ky), m=1/2, d0=1/10.

Fix an integer L>=8, a=r=1/(2L), epsilon=a^(9/10), and let k=2 pi N with N a positive integer that may grow arbitrarily. Every original Hessian remains in[.4,.6], independent of N.

At the actual source coordinates

    S=1, x=0, Z=-L/(N epsilon),

both S and S+epsilon Z have integer-multiple2pi phases. Hence

    Delta=L/(2N), x1=1/(4N)=pi/(2k),
    t0=1, t1=1+a C0/k, C0=pi/4+1/10.

Consequently

    g'(x)=.6, g'(x1)=.5,
    g''(x)=0, g''(x1)=-.1k,
    g'(t0)=.6, g'(t1)=.5+.1cos(aC0).

The EXACT second derivative of(1) in x is

    Psi_xx=ra[.1k g(t1)
                    +a{(.6)^3-(.5)^2(.5+.1cos(aC0))}]
           >=.05 ra k.                               (5)

Here g(t1)>.5 and the braced expression is positive. Thus the actual full-gradient first of this natural lift is unbounded as N grows, while r,a,epsilon and every original first bound remain fixed. After the actual x=cS+sU pullback, its U/U derivative is s^2 times(5), so the issue persists on the original Gaussian tape.

The x companion itself satisfies

    G_x=ra[.6*.5-.5 g(t1)]=ra[.05-O(a/k)],             (6)

whereas the native E at this same record is O(ra/k). The original pointwise bound |E|<=ra^2 epsilon|Z| is unchanged. This rules out a pointwise one-E-factor certificate for the companion as well as a uniform gradient first. It does not alone prove an Lp lower bound; a separate source-energy audit may integrate this phenomenon.

## 3. A six-VALUE ambient gradient lift avoids third jets

There is an exact finite lift that uses only original VALUES, at the price of keeping the nonlinear graph constraint explicit. Introduce EIGHT independent coordinate blocks

    Xi=(S,Z,x0,x1,p0,p1,d,lambda), each in R^d.

Define an analytical potential L (potential VALUES are not queried):

    L(Xi)/r = Vt(S+aM p0)-Vt(S+aM p1)
          + V0(x0)-V0(x1)-x0·p0+x1·p1
          + a[Vt(S)-Vt(S+epsilon Z)]-S·d
          + lambda·(x1-x0-M*d).                        (7)

Its entire gradient is an executable original-gradient VALUE program:

    G_S/r =gt(S+aM p0)-gt(S+aM p1)+aDelta-d,
    G_Z/r =-a epsilon gt(S+epsilon Z),
    G_x0/r=g0(x0)-p0-lambda,
    G_x1/r=-g0(x1)+p1+lambda,
    G_p0/r=aM*gt(S+aM p0)-x0,
    G_p1/r=-aM*gt(S+aM p1)+x1,
    G_d/r=-S-M lambda,
    G_lambda/r=x1-x0-M*d.                              (8)

All query arguments in(8) are affine in the ambient Xi. It uses exactly six original VALUES before identical-query caching: two outer gt, two g0, and gt(S),gt(S+epsilon Z). Its full first is at most C r, for a,epsilon<=1 and ||M||<=1, with a numerical C (C=16 suffices). There are no third jets. It is anchored at the full origin. Its ordinary Gaussian profile is O_p(r sqrt(d)) on a standard ambient8d tape, with caller bounded by the literal original-gradient caller sums. FIRST/adjoint verification uses the six original HVP sites; no HVP result is used to produce(8).

For clarity, among its nonzero Hessian blocks are

    G_(S,p0)=ra Ht(S+aMp0)M,
    G_(S,p1)=-ra Ht(S+aMp1)M,
    G_(p0,p0)=ra^2 M*Ht(S+aMp0)M,
    G_(p1,p1)=-ra^2 M*Ht(S+aMp1)M,
    G_(x0,x0)=rH0(x0), G_(x1,x1)=-rH0(x1),
    G_(x0,p0)=-rI, G_(x1,p1)=rI,
    G_(S,Z)=-ra epsilon Ht(S+epsilon Z),
    G_(S,d)=-rI, G_(d,lambda)=-rM,
    G_(x0,lambda)=-rI, G_(x1,lambda)=rI,

with the exact transpose blocks and the remaining diagonal S/Z blocks obtained from(8). The desired terminal/ancestor factors occur in this genuine symmetric matrix, along with real mixed constraint blocks.

## 4. Restriction to the actual nonlinear graph

At the exact actual source input define the feature map

    Gamma(S,U,Z):
      x0=cS+sU,
      d=a[gt(S)-gt(S+epsilon Z)],
      x1=x0+M*d,
      p0=g0(x0), p1=g0(x1), lambda=0.                 (9)

It uses the original aliases and finite version. On this graph (8) gives

    G_S(Gamma)=E,
    G_x0(Gamma)=G_x1(Gamma)=G_lambda(Gamma)=0,
    G_Z(Gamma)=-ra epsilon gt(S+epsilon Z),
    G_p0(Gamma)=r[aM*gt(t0)-x0],
    G_p1(Gamma)=r[-aM*gt(t1)+x1],
    G_d(Gamma)=-rS.                                   (10)

In particular the sum of the two p-companions is

    G_p0(Gamma)+G_p1(Gamma)=aM*E+rM*d.                 (11)

The second term has the actual ra epsilon scale, generally a factor1/a larger than the K3 remainder's ra^2 epsilon scale. The individual p blocks also contain the known linear x0,x1 terms and ordinary unmarked terminal profiles. G_d is a known Gaussian row, but is correlated with the original E record. These facts must be carried rather than replacing the entire lift's energy by e.

At the literal actual Gamma record, all six VALUES in(8) can be shared with evaluation of E. The three zero constraint blocks are exactly zero with the same cached arithmetic. Thus this is a genuine finite six-VALUE ambient lift and a correct algebraic readout of the original source.

## 5. The remaining descendant is the constrained input law

The actual feature map Gamma is nonlinear and its law is supported on the exact constraints in(9). It is not a standard8d Gaussian or a known affine image of one. Consequently the stationary genuine-gradient compiler for the ambient source G cannot simply be evaluated on Gamma and retain its standard-input theorem.

There are two distinct attempted uses:

1. Freeze Gamma as a coefficient root and introduce independent ambient Gaussian heat. This gives an admissible genuine-gradient source, but its selected factors are at newly perturbed, independent ambient coordinates. In particular p0,p1 no longer equal g0(x0),g0(x1), and d no longer equals the same terminal contrast. The selected product is the independently modified query law, not the K3 whole-input heat coefficient.
2. Substitute Gamma into the active source. The physical block is exactly E, but its derivative is the nonsymmetric product DG(Gamma)D Gamma. The zero constraints differentiate to the original relations Dp_i=H0(x_i)Dx_i and Dd=a[Ht(S)DS-Ht(S+epsilon Z)(DS+epsilon DZ)]; these contain exactly the correlated jets required by the open two-sided word. They are not known constant Gaussian rows. Taking the FULL gradient pullback D Gamma*G(Gamma) restores mathematical gradient symmetry only by adding those derivative companions; it requires HVP outputs and brings back unbounded higher-jet firsts in a subsequent derivative.

More explicitly, at fixed independent S,Z the x block of this FULL pullback is G_x from(2) plus r[H0(x1)x1-H0(x)x]. In the scalar fixture of §2 its x derivative is exactly Psi_xx-r(.1+.05 pi). Thus it remains unbounded with k; the additional constraint companions do not cancel the demonstrated higher-jet row.

The vanishing constraint coordinates therefore do not by themselves furnish a Gaussian-compatible active selector or a retained-root current estimate. A valid further generator would have to provide a weak joint comparison on this constrained law, or cancel/retain the extra channels with their actual energy, first and root frames. The known affine-query case is a legitimate escape when Gamma's relevant pullback is a supplied fixed linear map. It does not apply automatically to(9).

This isolates a specific new descendant, not a vague mixed-DAG difficulty: a small-first genuine-gradient ambient VALUE source together with the exact nonlinear feature constraints (9), whose physical readout is the original paired E. Treating the feature constraint as a free stationary input is the missing operation.

## 6. Diagnostics and scope of the returned component

The adjacent checker verifies1,138 cases: the entire eight-block gradient
against its potential, ambient Hessian symmetry and the16r bound, all
same-record constraint cancellations, the p-companion sum, and the exact
unbounded-first fixture. The returned constructive component is the ambient
six-VALUE gradient(8), together with its explicit graph map(9) and all
companions(10). The missing constrained-input comparison is not included in
this return. R's separate source-energy obstruction, SHA
ad5b6d24b0ab4bb2ed157b04cba99266ca75e02a4641bc52c02db590fb96c4b1,
integrates the natural companion and is separate from the pointwise/first
argument here.
