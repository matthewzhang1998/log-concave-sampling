# Known nonsymmetric affine query: exact one-pair curl extraction

New bounded addendum to the one-pair adapter (source pin0ad6959e). All matrix
adjoints are Euclidean. The map B below is known and fixed under caller
and Gaussian differentiation. Its known linear actions and fills are paid.

Let ||B-I||<=a<1/2 and let g=grad V be an original unit-first gradient. The
whole same-record source is f(W)=r[g(BW)-g(W)]. Define

    S=(B^(-1)-I)/a,     ||S||<=1/(1-a)<=2,
    h_B(x)=B* g(Bx),    Lip h_B<=||B||^2.

When the source is supplied as B=I+aK with known ||K||<=1, compute
S=-B^(-1)K. This equivalent formula avoids numerical subtraction followed by
division by a. If only B itself is encoded, its actual a^(-1) matrix-precision
amplification must be charged; no exact known arithmetic is granted for free.

h_B is an actual genuine gradient regardless of symmetry of B. At a fixed
whole-input clock cX+vZ let

    J_B(X)=E D h_B(cX+vZ),  J_I(X)=E Dg(cX+vZ).

Both are symmetric and the exact finite-source heat derivative is

    J_f=r[B^(-*) J_B-J_I].

Consequently its curl equals

    C_f=r[B^(-*) J_B-J_B B^(-1)]
       =ra[S* J_B-J_B S].                              (1)

Every factor uses the SAME X and original covariance cX+vZ. The r g(W)
terminal cancels from curl because its heat Jacobian J_I is symmetric. There
is no independent primitive heat inserted after B.

## Executed genuine-gradient block

At fixed X form the original-VALUE source

    h_(X,v)(z)=(rho/v)[h_B(cX+v z)-h_B(cX)].

Its complete active radius is at most L0 rho for numerical L0>=||B||^2,
its selected matrix is s0 rho J_B, and its old-X first is at most
2 L0 rho c/v. Use one full gradient pair on two copies of this source.
For C^2>=1+||S||^2 set

    U=[S*,-I]/C,       V=[I;S]/C.

The norms are at most one. Known input and output fills complete the
stationary channel exactly as in the one-pair adapter. Its selected matrix is

    K=(s0 rho/C^2)[S*J_B-J_BS].                          (2)

This is skew even when S is nonsymmetric. The conditional-passive covariance
is I-KK*. Only the two source gradients, known B/B* maps and known fills are
executed; no matrix J_B or its derivative is queried.

The actual marker source is [f(cX+vz)-f(cX)]/v, with all original aliases.
Its selected coefficient is M=s0 J_f and its integrated centered energy is
at most e/v, with e the energy of the WHOLE f. Hence (1)--(2) fit the supplied
positive fork with tau=ra. For a target covariance -q O_X, where
O_X=-Sym(J_f C_f), choose

    b c_side=-q ra C^2/(2 s0^2 rho).

The retained body first is O(ra/b), root first O(ra/(b v)), marked deletion
O(b e/v), and finite side prior O((ra/(b rho)) delta_pair), all with the
original positive-clock weights attached to the complete outputs. Full
root/current and variance ownership are exactly those of the held affine
fork, under its explicit old covariance-stage and independent-keep premises.

## Root coefficient frame and scope

D_X J_B[h] has HS norm at most L0 |h|/v by source-derivative commutation and
Gaussian first-chaos Bessel. Equation(1) therefore gives
||D_X C_f[h]||HS<=C ra |h|/v. These are the frames needed by the whole-root
Price current; no Hessian continuity or derivative oracle is being assumed.

This extends exact native heat matching from symmetric affine B to all
well-conditioned known near-identity affine maps. It does not authorize a
random B depending on X, a nonlinear ancestor, or substituting Dq(X) for B.
Such substitutions add derivative paths and change the gradient pullback.
Dense known inverse/fill arithmetic can be substantial; explicit rank-one or
scalar-block geometry is needed to claim the small native arithmetic cost.
The number of gradient-pair banks is one, but its duplicated source blocks
both cost original gradient calls.
