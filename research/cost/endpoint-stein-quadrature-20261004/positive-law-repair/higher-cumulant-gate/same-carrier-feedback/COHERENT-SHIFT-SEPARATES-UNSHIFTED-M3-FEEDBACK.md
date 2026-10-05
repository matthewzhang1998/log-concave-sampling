# Coherent mean shift blocks an unshifted m3 feedback shortcut

2026-10-05. Counterexample candidate, independent review requested. This rejects one exploratory shortcut, not the literal F3_Q source or the proposed covariance-only bias correction in the separate m3 gate note.

## Claim being tested

Let H=F1 and E2=F2−H on the same conditional Gaussian OU history. The tempting simplification

    R1 E[g(x−F2)−g(x−H)|x]
        ?=−R1[Dg(x) E(E2|x)]+O(A4 sqrt(D))           (1)

does not retain the coherent conditional mean shift already present in H. The smooth Hessian-bounded family below gives a remainder of order A2 at A=D^(−1/2), while A4 sqrt(D)=A3. Thus (1) is not a safe source reduction even with fixed public-log losses.

## 1. Reuse the proved coherent-shift family

Use c=1/4, epsilon=1/2, v=(1,...,1)/sqrt(D), A=D^(−1/2),

    b(t)=t+epsilon log cosh(t),
    g_D(x)=c A (b(x1),...,b(xD))+c A v tanh(v dot x).

It is C-infinity, g(0)=0, and 0<=Dg<=c(2+epsilon)A I<A I. The exact same-Gaussian-history asymptotics in `SAME-CARRIER-P3-MEAN-AND-FULL-COVARIANCE-GATE.md`, SHA256 a63c38f14206ec278f7c12a56542968a7d1675af6f550610e4539e500a653eca, give

    v dot H -> d,
    (v dot E2)/A -> e(U),
    d=c epsilon E log cosh(N)>0,
    e(U)=c[C0+integral_0^1 h(U_r)dr],
    C0=−d−epsilon c k,
    h(u)=tanh(u−d)−tanh(u),
    k=E[(integral_0^1 b(X_r)dr)tanh(X1)].            (2)

Here U_r=v dot X_r is a standard one-dimensional OU path and U1=v dot x. The first convergence is in L4, and the second in L2 with uniformly bounded L4 norms. The same Taylor/coordinate-average proof gives these moment strengthenings because b has linear growth and all Gaussian moments are finite.

Parity and Gaussian regression give

    k=(1/2)E[N tanh(N)]>0.                            (3)

In particular C0<0 and h(u)<0 for every real u.

## 2. Evaluate the exact feedback remainder

Define the actual vector remainder

    R_D=g_D(x−F2)−g_D(x−H)+Dg_D(x) E2.              (4)

Each coordinate H_i has L4 norm O(A), and each E2_i has L4 norm O(A2). The latter includes its rank-one term: A v_i=A/sqrt(D)=A2. The coordinatewise b part of (4), projected onto v, is O_L2(A3). Indeed b'' is bounded, so its ith remainder is at most C A (|H_i|+|E2_i|)|E2_i|; after the v projection, Minkowski costs A/sqrt(D) times D terms of order A3, equal to O(A3).

For the rank-one term, bounded tanh'' yields

    tanh(U1−v dot F2)−tanh(U1−v dot H)
       =−sech2(U1−v dot H)(v dot E2)+O((v dot E2)2).

Its quadratic remainder, multiplied by c A, is O_L2(A3). Therefore (2) gives

    (v dot R_D)/A2
       -> c[sech2(U1)−sech2(U1−d)] e(U)             (5)

in L2. This is the correct derivative at the shifted scalar argument U1−d, contrasted with the unshifted derivative in (1).

## 3. The conditional mean and outer resolvent cannot erase the gap

Conditioning the limit in (5) on the complete endpoint x leaves only U1, because the scalar projected Gaussian history is independent of all orthogonal endpoint coordinates. Conditional expectation is an L2 contraction, so no interchange of uncontrolled conditional limits is needed. The limiting conditional function is

    F(u)=c2 [sech2(u)−sech2(u−d)]
                     [C0+integral_0^1 P_r h(u)dr].   (6)

The bracket on the second line is strictly negative by (3), while the first bracket is nonzero except at u=d/2. Thus F is a nonzero bounded L2 function. The one-dimensional Gaussian resolvent R1 is injective on L2 because its Hermite multipliers are 1/(n+1)>0. Consequently ||R1 F||2>0.

The D-dimensional resolvent preserves the scalar projected subspace: R1[F(v dot x)](Z)=(R1^(1D)F)(v dot Z). Applying its L2 contraction to the error in (5) therefore proves

    ||R1 E[R_D|x]||_(L2(gamma_D)) >= c_* A2           (7)

for all sufficiently large D, with c_*>0 independent of D. Since A4 sqrt(D)=A3, (7) contradicts the claimed uniform remainder in (1), including any fixed polynomial-log multiplier.

## Scope

This separator concerns the unshifted derivative Dg(x) in the SMALL-FEEDBACK term m3−m2. It does not invalidate the independently proved K reduction, whose coefficient is already a conditional covariance of order A2 and whose extra weak-shift cost has a different one-energy budget. It does not refute a correctly resummed expansion around the coherent conditional mean, nor the literal finite F3_Q VALUE program that keeps the shifted ancestors.

No new finite sampler or general impossibility theorem is claimed. The safe m3 target remains R1 E[g(x−F2)|Gaussian endpoint x], with its full original history and current restoration.
