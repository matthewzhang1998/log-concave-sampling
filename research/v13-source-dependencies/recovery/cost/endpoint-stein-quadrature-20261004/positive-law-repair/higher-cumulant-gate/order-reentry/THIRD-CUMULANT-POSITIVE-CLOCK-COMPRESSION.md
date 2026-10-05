# Positive five-clock compression of the connected third cumulant

2026-10-04. Analytical coefficient and finite-clock result, pending independent audit. Original gradient derivatives below identify the target; only a separately admitted native VALUE tree may execute the coefficient action.

## 1. Exact compressed identity

Use R_k=(k-L_OU)^(-1)=int_0^1 q^(k-1)P_q dq. For original g=grad U, 0<=Dg<=A I, let

    F=Dg, B=R_2 F=Dv, T=DB,
    J=T dot B, U=R_3 J.

Index conventions before final symmetrization are

    J_jka=sum_b (partial_a B_jb)B_kb,
    (B dot U)_ijk=sum_a B_ia U_jka.

The connected-cumulant recurrence and full symmetry of DB give

    kappa_3(H|Z)=24 R_3 Sym_3[B dot R_3(T dot B)].       (1)

Indeed DC=4 R_3 Sym_(j,k)(T dot B), where only its TWO covariance slots j,k are symmetrized. The two product-rule terms become equal after the FINAL physical Sym_3, so the displayed raw J gives (1). Do not symmetrize its spatial derivative slot a with j,k before contracting B_ia at the outer point; that would change the matrix ordering at different points.

Expanding B=R_2F and T=D R_2F gives five positive scalar clocks: the two R_3 clocks q,s and the three primitive clocks r_1,r_2,r_3. The scalar factor is

    24 q^2 s^2 r_1 r_2^2 r_3,

and the original three-force tensor has leaf Dg, center D^2g and leaf Dg. The derivative at the center is taken after its positive Gaussian heat. No finite weighted H is cubed and no artificial same-clock diagonal is introduced.

## 2. Uniform proper-cut bounds for the derivative slot

For a unit spatial direction u and c=sqrt(1-r^2), gradient symmetry gives in the distribution/heat sense

    D_u[P_r Dg](x)
       =(r/c) E[(Dg(rx+cG)u)G^T].                    (2)

The first Gaussian-chaos Bessel inequality yields

    ||D_u[P_r Dg](x)||HS <=(r/c) A.

This is dimension free: one uses the vector Dg u, whose norm is at most A, rather than the full matrix-HS norm of Dg. Therefore

    ||T(x)||_(R^D->HS)<=A int_0^1 r^2/sqrt(1-r^2)dr
                            =(pi/4)A,
    ||T(x)||HS <=(pi/4) A sqrt(D).                    (3)

DB is fully symmetric, since v is a gradient. Every one-versus-two proper cut has the same bound. B has operator norm at most A/2. Consequently all proper cuts of T dot B, after physical symmetrization when needed, are bounded by C A^2, its HS norm by C A^2 sqrt(D), and analogous bounds for B dot R_3(T dot B) are C A^3 and C A^3 sqrt(D).

These statements hold pointwise in the captured carrier. They require only the bounded original Hessian and Gaussian heat, not a Hessian continuity modulus.

## 3. Finite positive quadrature with one marked energy

Let Q_k use the already admitted positive dyadic rule to approximate R_k, all interior nodes. Choose its ordinary uniform moment error delta0. Define

    B_Q=Q_2 F,
    T_Q=D Q_2 F,
    J_Q=T_Q dot B_Q,
    U_Q=Q_3 J_Q,
    K_Q=24 Q_3 Sym_3(B_Q dot U_Q).                    (4)

The two B slots may use the same deterministic quadrature version but their conditional expectation banks are independent when the product is realized. All appropriate source-specific versions are retained.

The exact first-moment calibration gives ||B_Q||op<=A/2. The dyadic bound sum w r^2/sqrt(1-r^2)<=C gives the same proper-cut estimate C A for T_Q as in (3), uniformly in node count. Its full HS bound is C A sqrt(D).

The derivative-aware resolvent lemma in `FIRST-RESOLVENT-CLOCK-QUADRATURE.md` gives

    ||B_Q-B||_(L2;HS)<=delta0 A sqrt(D),
    ||T_Q-T||_(L2;HS)<=sqrt(7delta0) A sqrt(D).        (5)

Set delta0<=delta^2/7, with delta<=1. Product telescoping retains one marked HS error. For example,

    (T_Q dot B_Q)-(T dot B)
        =(T_Q-T) dot B_Q+T dot(B_Q-B).

The first term costs ||B_Q||op times the first HS error. For the second, use the proper cut of T from its contracted B-index to the two other slots, costing C A times ||B_Q-B||HS. Thus

    ||J_Q-J||_(L2;HS)<=C delta A^2 sqrt(D).

A positive Q_3 has bounded L2 norm and a uniform proper-cut bound; applying (1) and telescoping the two remaining R_3 and B products similarly yields

    ||K_Q-kappa_3||_(L2(Z);HS)<=C delta A^3 sqrt(D).   (6)

There is no D^(3/2) coefficient allowance in (6). All products use one full HS factor and uniform proper-cut/operator bounds for the others. Fixed-order positive quadrature uses O(log^2(1/delta)) nodes per clock, so expanding (4) uses at most the fifth power of that public-log count. This is a coefficient-action target size; explicit D^3 tensor materialization is not authorized or charged as a gradient-query algorithm.

## 4. Literal Gaussian query genealogy and separate shields

At each five-clock node, the analytical tensor is the expectation of the three original forces on

    Y=qZ+sqrt(1-q^2)G0,
    X=sY+sqrt(1-s^2)G1,
    Q1=r1 Y+sqrt(1-r1^2)G2,
    Q2=r2 X+sqrt(1-r2^2)G3,
    Q3=r3 X+sqrt(1-r3^2)G4,

where G0,...,G4 are independent standard D-vectors. Thus the two inner sites share the SAME X, and all three share the inherited Y. This is not an independent replacement of the mixed tensor.

Each original query already has its OWN independent Gaussian floor. Split it into equal halves with

    sigma_i=sqrt(1-r_i^2)/sqrt(2),
    Q_i=Q_i,coarse+sigma_i Z_i,

where the three Z_i are independent of the complete coarse bank. This preserves the exact query law. There is no need to lower all three shields to the minimum of q,s,r1,r2,r3.

For a native vertex at the fixed coarse query, the literal gradient source

    f_i(u)=[g(Q_i,coarse+sigma_i u)-g(Q_i,coarse)]/(A sigma_i)

has private first at most one, zero by saved-anchor reuse, and coarse-caller first at most 2/sigma_i. The leaf target is Dg/A; the center's first filtered public response targets (sigma_2/A)D^2g. Therefore only the center normalization contributes 1/sigma_2 to the desired product amplitude.

An all-A three-vertex native normalization places the known coefficient proportional to w A/sigma_2 at the root and radii A at the other vertices; w is the complete positive five-clock coefficient. Source radii are small because each node weight near r2=1 is O(1-r2), so w/sigma_2 is uniformly bounded up to fixed known factors.

The root/coarse first ledger has the form

    C A sum_nodes w/[sigma_2 sigma_i], i=1,2,3.        (7)

For i=2 its one-clock inverse-square sum is logarithmic. For i=1 or i=3 the two separate inverse-first sums are bounded by geometric sums. Thus (7) is a public-log A bound. It keeps the complete coarse roots and their actual paths; a fixed-coarse private-radius statement alone would not supply it.

This section is a source/genealogy specification. The native tree's complete output-law, original-query replay count, source-index response identities, positive variance allocation and all feedback still require their independent native-packet theorem. Formula (6) does not execute a tensor oracle or prove that theorem.

## Scope

Returned: exact connected third-cumulant identity, dimension-safe proper-cut/one-HS bounds, positive polylogarithmic five-clock approximation, and a literal original-gradient Gaussian query/shield genealogy. Not returned: a full fourth-order posterior bridge, a standalone executable tensor, or a repeatable all-order finite law compiler.
