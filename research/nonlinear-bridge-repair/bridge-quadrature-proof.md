# Gaussian bridge and positive operator quadrature

## Verified conclusion

For every finite dimension D and every 0 < delta < 1, there exist positive nodes s_i > 0 and positive weights p_i with

- sum p_i = 1;
- sum p_i exp(-s_i) = 1/2;
- the number of nodes is O(log^2(1/delta)), with absolute dimension-free constants;
- || integral_0^infinity exp(-s) T_s ds - sum_i p_i T_{s_i} ||_{L2(gamma_D) -> L2(gamma_{3D})} <= delta.

Here T_s h(x,N,M) = E_L h(exp(-s)[x+2sN+2s(s-1)M]+c(s)L), with the stated c(s). The same conclusion holds for finite-dimensional Hilbert-valued outputs.

The bridge formula is correct for each fixed s. A single common residual L does NOT reproduce the simultaneous multi-time OU law. Nothing in this operator quadrature proof needs that stronger assertion.

## 1. Exact bridge

Use a stationary OU realization, for t >= 0,

X_t = exp(-t)x + sqrt(2) integral_0^t exp(-(t-u)) dB_u,

where x is standard Gaussian and independent of B. Stochastic Fubini gives

H1 = x/2 + (1/sqrt(2)) integral_0^infinity exp(-u) dB_u,
H2 = x/4 + (1/sqrt(2)) integral_0^infinity (u+1/2)exp(-u) dB_u.

Set

N = sqrt(2) integral_0^infinity exp(-u) dB_u,
M = 2sqrt(2) integral_0^infinity (u-1/2)exp(-u) dB_u.

The identities integral exp(-2u) du=1/2, integral u exp(-2u) du=1/4, and integral u^2 exp(-2u) du=1/4 show that x,N,M are independent standard D-Gaussians. Therefore H1=x/2+N/2 and H2=x/4+N/2+M/4.

At fixed s, the covariance coefficients of X_s against x,N,M are respectively

exp(-s), 2s exp(-s), 2s(s-1)exp(-s).

Gaussian orthogonal projection therefore gives an independent residual c(s)L, where

c(s)^2=1-exp(-2s)[1+4s^2+4s^2(s-1)^2].

The residual is nonnegative by projection, and the independent complex-disk certificate below also proves it is strictly positive when s>0.

## 2. A dimension-free complex disk certificate

The requested relative disk radius 1/1024 can be strengthened to eta=1/64. Put

a=65/64, b=63/64,
q(z)=exp(-z)(1,2z,2z(z-1)),
E(s)=sum_{j=0}^8 (2bs)^j/j!.

For |z-s|<=s/64,

||q(z)||^2 <= exp(-2bs) P_-(s)   if 0<=s<=1,
||q(z)||^2 <= exp(-2bs) P_+(s)   if s>=1,

where

P_-(s)=1+4a^2 s^2[1+(1-bs)^2],
P_+(s)=1+4a^2 s^2[1+(as-1)^2].

This follows directly from Re z>=bs, |z|<=as, and |z-1|<=|s-1|+s/64.

For a fully rational positivity certificate, expand R(s)=(E(s)-P_-(s))/s in the degree-7 Bernstein basis B_{k,7}(s)=binom(7,k)s^k(1-s)^(7-k). Its coefficients have common denominator

24629060462182400

and numerators, in order k=0,...,7,

48488462784921600,
26273173943091200,
15076247889510400,
12524793794396160,
16419361684717568,
24828792912732160,
36265100707349760,
50076245130694685.

Every numerator is greater than half the denominator. Hence E(s)-P_-(s)>s/2 for 0<s<=1.

Next expand E(1+t)-P_+(1+t)=sum_{k=0}^8 C_k t^k. Its coefficients have common denominator

703687441776640

and numerators, in order k=0,...,8,

1430749860876991,
4011516377289976,
3642743757137380,
228474934592712,
2444527082810,
1071611807151816,
279243854976996,
47517861678840,
3938980639167.

Every numerator is positive. Therefore E(s)>P_+(s) for s>=1. Since E(s)<=exp(2bs), these two certificates give

||q(z)||<1 for s>0 and |z-s|<=s/64.

At s=0 the disk is the singleton 0 and ||q(0)||=1.

The companion exact-arithmetic verification script verifies every coefficient above and the strict inequalities. No floating-point or Sturm assertion is needed.

## 3. Holomorphic operator extension

Complexify the Gaussian L2 spaces and decompose into homogeneous Hermite chaoses. On chaos of degree k, T_s is the symmetric k-fold tensor power of the linear embedding

j_s : C^D -> C^{3D},  v -> (q_0(s)v,q_1(s)v,q_2(s)v).

Define T_z by precisely the same chaos formula. Its degree-k block has norm ||q(z)||^k. Consequently ||T_z||=1 wherever ||q(z)||<=1 (the constants give equality).

At every z in the open union of the relative disks, ||q(z)||<1. On a sufficiently small compact neighborhood it is bounded by r<1. Truncating at chaos degree n then has operator-norm error at most r^(n+1). Each truncated block operator is holomorphic, so locally uniform convergence proves operator-norm holomorphy.

Do NOT try to obtain this result by choosing an analytic square root c(z) inside a complex Gaussian expectation. The Hermite definition avoids that issue and preserves the correct complex bilinear covariance automatically.

The operator-valued integral K=integral_0^infinity exp(-s)T_s ds exists in operator norm: T_s is norm-continuous on (0,infinity), bounded by 1, and exp(-s) is integrable. Norm continuity at s=0 is unnecessary and generally false.

## 4. Positive quadrature with exact mass

Construct first Q with ||Q-K||<=epsilon, then set epsilon=delta/9. It is enough to consider 0<epsilon<=1.

Let h=epsilon/8, H=log(8/epsilon). Replace the head [0,h] by an atom at h/2 with exact weight 1-exp(-h), and the tail [H,infinity) by an atom at H+1 with exact weight exp(-H). Their respective operator errors are at most 2(1-exp(-h))<=epsilon/4 and 2exp(-H)=epsilon/4.

Partition [h,H] into dyadic panels [A,B] with B<=2A; the number of panels is

K_panels=ceil(log_2(H/h)).

Use the n-point Gaussian quadrature for the positive measure exp(-s)ds on each panel. These rules have positive interior nodes and weights, are exact on polynomials of degree <=2n-1, and their weights sum exactly to the panel's exponential mass.

For completeness, existence and positivity are elementary: construct the monic degree-n orthogonal polynomial for exp(-s)ds. Its n roots are simple and lie inside the panel (otherwise multiplying the factors corresponding to its interior sign changes contradicts orthogonality). Interpolation at these roots and division by the orthogonal polynomial proves degree-(2n-1) exactness. For the Lagrange polynomial ell_i of degree n-1, the corresponding weight equals integral exp(-s)ell_i(s)^2 ds>0, by exactness.

Take rho=65/64. The Bernstein ellipse for [A,B] with parameter rho is contained strictly within the union of the complex disks from Section 2. Indeed, write its boundary as

z=C+d[(rho+rho^-1)cos(theta)/2 + i(rho-rho^-1)sin(theta)/2],

where C=(A+B)/2, d=(B-A)/2<=A/2. The point t=C+d cos(theta) lies in [A,B], and

|z-t|<=d(rho-rho^-1)/2<d(rho-1)<=A/128<t/64.

The ellipse interior is contained in the same constant-radius tube around [A,B]. Thus T_z is holomorphic in a neighborhood of the closed ellipse and bounded there by 1.

Banach-valued Chebyshev truncation at degree m gives an approximating operator polynomial P_m with

sup_{s in [A,B]} ||T_s-P_m(s)|| <= 2rho^(-m)/(rho-1).

This follows from the usual Cauchy integral coefficient estimate ||a_k||<=2rho^(-k), which works unchanged for bounded-operator-valued functions.

With m=2n-1, positivity, exact polynomial integration, and exact mass show that each panel's quadrature error is at most

4 [panel exponential mass] rho/(rho-1) rho^(-2n).

Summing the panel masses gives a total middle error <=4rho/(rho-1)rho^(-2n). Choose

n=ceil(log(8rho/((rho-1)epsilon))/(2log rho)).

The middle error is then <=epsilon/2. Head, middle, and tail together give ||Q-K||<=epsilon. All weights are positive and sum to exactly 1. The node count is n K_panels+2=O(log^2(1/epsilon)); constants do not depend on D.

## 5. Exact exponential moment by one positive extra node

Let m=sum_i p_i exp(-s_i). Apply Q-K to the unit-norm first Hermite h(y)=y_1. Its x_1 coefficient is m-1/2, so |m-1/2|<=||Q-K||<=epsilon.

If m=1/2, retain Q. Otherwise choose

r=1/4 and s_*=log 4 if m>1/2;
r=3/4 and s_*=log(4/3) if m<1/2.

Define lambda=(1/2-m)/(r-m) and Q'=(1-lambda)Q+lambda T_{s_*}. Then 0<lambda<1, lambda<=4epsilon, all weights remain positive, total mass remains 1, and the exponential moment is exactly 1/2. Also

||Q'-K|| <= (1-lambda)epsilon + lambda||T_{s_*}-K||
             <= epsilon+2lambda <=9epsilon.

At most one node is added. Taking epsilon=delta/9 proves the stated conclusion.

## Scope and pitfalls

1. Fixed-s Gaussian law is exact; simultaneous OU multi-time law with a common L is not asserted.
2. Analytic continuation must use Hermite blocks, not a complex square root/noise integral.
3. The interval Gauss rule must be weighted by exp(-s) if exact total exponential mass is claimed. Unweighted Gauss-Legendre followed by multiplying weights by exp(-node) is only approximately mass-exact unless normalized.
4. Operator holomorphy starts on s>0; no norm-continuity assertion at s=0 is needed.
5. Exact moment repair here concerns only sum p_i exp(-s_i)=1/2. It does not make the whole rule exact on degree-2 Hermite chaos, and that stronger property is unnecessary for the specified quadratic source.
6. All existence and accuracy statements use real-arithmetic weights. A finite-precision implementation needs a separate positive normalization/moment-correction tolerance analysis.
