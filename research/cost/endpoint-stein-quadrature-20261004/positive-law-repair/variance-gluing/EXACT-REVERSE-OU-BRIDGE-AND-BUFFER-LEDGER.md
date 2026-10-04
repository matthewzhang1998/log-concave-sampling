# Exact reverse-OU bridge and the conditional buffer ledger

2026-10-04. A positive affine gluing component and an exact current ledger. This is not a completed high-order posterior transition.

## 1. Exact transition on the same OU endpoint path

Let `mu(dx) proportional to exp(-|x|^2/2-U(x)) dx`, with `g=grad U`, `g(0)=0`, and `0<=Dg<=A I`, A<=1/2. Its OU path rho_t is the law of `tX+sqrt(1-t^2)N`, X~mu. For `0<=r<t<=1`, couple the path by its usual forward OU Markov kernel from t down to r:

    Y_r=(r/t)Y_t+sqrt(1-r^2/t^2)N_forward.

This is a concrete joint law, not independent noise chosen separately at both times. Fix its exposed endpoint Y_r=z. Put

    s^2=1-r^2, Delta=t^2-r^2,
    a_(r,t)=r(1-t^2)/(t s^2),
    b_(r,t)=Delta/(t s^2),
    v_(r,t)=(1-t^2)Delta/(t^2 s^2).

The conditional law of X given z is exactly

    nu_(r,z)(dx) proportional to exp(-|x-rz|^2/(2s^2)-U(x)) dx.

Gaussian conditioning yields the exact positive same-endpoint transition

    Y_t = a_(r,t) z+b_(r,t) X+sqrt(v_(r,t))N,
    X~nu_(r,z), N~gamma independent conditional on z.        (1)

At t=1, the additive N term vanishes and Y_t=X. Formula (1) retains every conditional covariance and higher cumulant of X. It is not a Gaussian mean replacement.

Define the exact conditional velocity and curvature matrix

    m_r(z)=E_nu g(X),
    H_r(z)=E_nu Dg(X)-Cov_nu(g(X)).                          (2)

These are analytical objects. No HVP, covariance or conditional-mean oracle is added. At r=0 the formulas still make sense, while the derivative identity below is written only for r>0.

Conditional Gaussian integration by parts gives

    E_nu X=rz-s^2 m_r(z),
    Cov_nu X=s^2 I-s^4 H_r(z),
    D_z m_r(z)=r H_r(z), r>0.                              (3)

The inverse-Hessian covariance inequality for the strongly convex conditional potential also gives `0<=H_r<=A I`. One direct justification is

    Cov_nu(g)<=E_nu[Dg(s^(-2)I+Dg)^(-1)Dg]<=E_nu Dg;

all matrices in the pointwise middle expression are functions of the SAME local symmetric Dg. This does not commute Hessians at different locations.

Set

    c=r/t, d=Delta/t, v_0=Delta/t^2=1-c^2.

Then the exact conditional mean and covariance of (1) are

    E[Y_t|z]=c z-d m_r(z),
    Cov(Y_t|z)=v_0 I-d^2 H_r(z).                           (4)

Since d^2=v_0 Delta and H_r<=A I, the exact covariance has relative gap at least 1-A Delta>=1/2. Every conditional cumulant of order k>=3 is `b_(r,t)^k` times the corresponding cumulant of X. They are actual generated descendants of (1) and are not declared zero.

## 2. A unit-buffer mean can be glued positively

Suppose an ACTUAL fresh conditional law service M_r has

    ||W2(Law(M_r|z),N(m_r(z),I))|| <= epsilon_r(z),

with z retained and all its fine/source/clock banks integrated, as in the preceding conditional-velocity service. Execute, using one fresh independent Gaussian N,

    K^mean_(r,t)(z)=c z-d M_r(z)+sqrt(v_0-d^2)N.              (5)

This is positive for EVERY 0<=r<t<=1 because

    v_0-d^2=v_0(1-Delta)>=0.

Conditionally on z, it satisfies

    W2(Law(K^mean_(r,t)(z)),
           N(cz-dm_r(z),v_0 I)) <= d epsilon_r(z).           (6)

The proof is exactly conditional convolution and affine scaling. It does not identify the unobserved comparison Gaussian in M_r with an executed internal root. In particular, (5) does not subtract the mean program's own Gaussian carrier. It reads only its complete output, the already exposed caller z, and a new independent N.

Thus a unit Gaussian buffer is compatible with a same-endpoint affine stage without creating an unpriced variance increase. But the target on the right of (6) is a Gaussianized transition. It differs from the exact transition (1) by the REAL covariance current

    -d^2 H_r(z),                                          (7)

and all the k>=3 descendants specified after (4). Mean accuracy alone cannot remove them.

## 3. Positive room for a separately supplied covariance reserve

For a bounded increment Delta<=1/4, a stronger covariance reserve could be added independently of M_r conditional on z. Its necessary reference covariance is

    C_res(z)=v_0 I-d^2[I+H_r(z)].                           (8)

It obeys

    C_res(z)>=v_0[1-Delta(1+A)]I >=(5/8)v_0 I.

Therefore positivity and a fixed relative Gaussian gap are NOT the obstruction to a full second-moment join on such increments. If an admitted root-retained positive covariance service really targeted (8), its independent sum with `cz-dM_r` would be W2-close to a Gaussian REFERENCE with the exact mean and covariance (4). The error is at most d times the mean-service conditional error plus the covariance-service conditional error. The ACTUAL sum is not asserted to have those exact moments when either law service has nonzero error; W2 closeness does not erase that distinction.

The required source is the whole H_r in (2), including both E Dg and the true same-conditional-law covariance of g. It is not a forward square, an unconditional trace, or just the derivative of an unqualified finite conditional mean. A VALUE source for a derivative must be supplied by its legal native coefficient filter/pair theorem, and every same-root product and finite-restoration error must be retained. A covariance match still does not remove the higher cumulants of (1).

## 4. Exact quadratic test

For `g(x)=Bx`, with any fixed symmetric `0<=B<=A I`,

    m_r(z)=r B(I+s^2 B)^(-1)z,
    H_r=B(I+s^2B)^(-1).

Both (1) and its exact conditional law are Gaussian. The mean-only stage (5), even with a PERFECT M_r, therefore has the covariance defect exactly `d^2 H_r`. This is a full matrix covariance discrepancy, not a scalar trace artefact. At r=0,t=1 its mean is zero and its covariance is I, whereas the correct endpoint covariance is `(I+B)^(-1)`.

The new velocity service can consequently have arbitrarily accurate mean on quadratics while (5) still fails already at first order for a large step. The separately supplied reserve (8) would exactly repair this quadratic transition; it is not assumed available by this note.

## 5. Alternative diffusive path, with a paid stage buffer

There is another exact way to give the unit mean-law buffer a role. The OU density path rho_r solves the SDE

    dY_r=[-Y_r-(1+r)m_r(Y_r)] dr+sqrt(2)dW_r,
    Y_0~gamma.                                            (9)

Indeed `grad log rho_r=-x-rm_r`, so its Fokker-Planck right side is exactly `div(rho_r m_r)=partial_r rho_r`. The positive-law audit proves `Lip(m_r)<=r[A+(1-r^2)A^2]`, with continuous endpoint limits and finite Gaussian moments. Thus the time-dependent drift is globally Lipschitz in x, measurable in r, and has linear growth. The SDE has a unique finite-second-moment solution. The displayed density solves its Fokker-Planck equation, so uniqueness identifies its marginals as rho_r, including the terminal mu. No third derivative of U is used.

For a frozen-force exponential-Euler step with 0<h, r+h<=1, let `c_h=exp(-h)` and `d_h=(1+r)(1-c_h)`. The positive program

    c_h x-d_h M_r(x)+sqrt(1-c_h^2-d_h^2)N                  (10)

targets `N(c_hx-d_hm_r(x),(1-c_h^2)I)` with error at most `d_h epsilon_r(x)`. For `h<=log(5/3)`, its radicand is nonnegative uniformly over r in [0,1]; for h<=1/4 it has a fixed positive gap relative to h. Again only the complete M_r output is read.

Here is the complete finite-stage argument. Set q_h^2=1-c_h^2-d_h^2. Since 1+r<=2,

    q_h^2>=(1-c_h)(5c_h-3).

Therefore h<=log(5/3) makes q_h^2>=0. If h<=1/4, use c_h>=3/4 and 1-c_h>=3h/4 to obtain q_h^2>=9h/16>h/2; also q_h^2<=2h. There is then a genuine fresh numerical Gaussian reserve of order sqrt(h).

At the SAME exposed caller x, couple the complete output M_r to `m_r(x)+G_ref` at its stated conditional error. The reference G_ref is standard Gaussian conditional on every captured caller/parameter, so it may be represented as fresh independent reference noise. It is NOT declared equal to any executed fine/source/clock Gaussian inside M_r. Draw the additional N independently of this entire coupling. The reference output is

    c_h x-d_h m_r(x)-d_h G_ref+q_h N.

Its two noises have exact covariance `(d_h^2+q_h^2)I=(1-c_h^2)I`; hence their normalized sum is a fresh standard Gaussian relative to the exposed caller. The actual/reference output gap is exactly `-d_h(M_r-m_r(x)-G_ref)`. This proves the conditional W2 cost `d_h epsilon_r(x)` with no missing covariance cross term. For random callers, square and integrate at the actual entering caller law. The caller and all declared external labels remain fixed at this cut; the private banks of M_r do not.

The conditional mean service and its imported smallness/log guards are unchanged: d_h is applied AFTER its complete output, not by shrinking its buffer inside the compiler. A stage costs one complete M_r call plus one fresh D-dimensional Gaussian block and known arithmetic. Its original VALUE and first/adjoint counts are exactly those of that complete call; every stage is another full call unless it has the identical captured caller and complete semantic key. The new private tape is the entire M_r tape plus this N block, not just two Gaussians. Known square roots and small h receive their actual precision allowances.

If the completed service has actual incoming first at most L_M=Lambda A r, then the actual stage incoming first is at most `c_h+d_h L_M`. Under `Lambda A<=1/4`, this is at most `(1+c_h)/2<1`. Its fresh first is O(sqrt(h)) using the complete mean program's unit carrier and O(Lambda A) residual first. These are optional actual-graph consequences of the separately supplied finite first port, not conclusions drawn merely from a conditional Gaussian law.

This solves the BUFFER ALLOCATION at that finite stage. It does not prove a high-order time integrator for (9). Ordinary fixed-force Euler accuracy must be paid. Even granting the sharper velocity service, an unproved high-order interpolation cannot be inserted at no cost. All time steps, independent complete M_r calls, Gaussian roots, original replay counts and numerical allocations are additional work.

## 6. What remains for actual order amplification

Either route now has an explicit same-endpoint gate:

- Reverse-OU transitions require a positive compiler for the exact conditional covariance (7) and all higher conditional descendants of (1), or a different complete non-Gaussian conditional transition that handles them together.
- The diffusive formulation requires a positive high-order finite transition for (9), with its actual unit-buffer mean services and no hidden smooth-force or strong-mean oracle.

An arbitrary finer time grid supplies neither a sublinear cost recurrence nor an all-order native compiler. Conversely, the buffer by itself is no longer an undefined objection: (5), (8), and (10) state precisely where the variance goes, when the reserve is positive, what caller is retained, and which remaining currents have not been canceled.
