# Test of intermediate Gaussian reserves

2026-10-04. Bounded test of one explicit escape from the singular retained-path Gram. This is not a posterior-law or arbitrary-algorithm lower bound.

## Answer

Adding independent small Gaussian reserves does make a transverse path covariance subtraction positive. Under the currently available **conditional C2** source contract, however, its necessary reserve size and the subsequent force-mean price reproduce exactly the existing strong earlier-layer query exponent. A graded allocation can spend more noise at remote layers, but the layer immediately before the last nonlinear force retains that exponent.

Correcting the reserve's effect by a literal same-record value telescope is legal, but the corrected observer reconstructs the original unbuffered variable. One cannot then continue to count the reserve as an unread positive gap for that observer.

An improved marginal law estimate, or a new positive joint counterpacket for the reserve-bias/current, could escape this conclusion. Neither is supplied by adding the reserves alone.

## 1. Size of the gap and its later price

Use the notation of `SCALAR-RANDOMIZED-QUARTER-FLOW-AND-ALL-LAYER-DEBTS.md`. In normalized state coordinates, layer j has a stratified clock residual with characteristic energy

    e_j=C_j a N_j^(-3/2).

Take a fixed finite set of retained path observations, whose size may depend on the requested order. Let w be a transverse zero-Gram direction such as the explicit one in `ALL-LAYER-MARTINGALE-CURRENTS-AND-PATH-GAP.md`.

If a new independent known Gaussian path reserve has standard deviation tau_j in that direction, a positive covariance budget capable of absorbing the clock error requires

    tau_j^2 >= (1+gap) Var(w dot D_j).

Here e_j is an upper/characteristic scale; a necessity claim additionally needs a nondegenerate transverse fixture at that scale. Such a family is obtained by replacing the single impulse in the companion path-gap test by N independent stratified times in a fixed interval between t_1 and t_2, each with weight h/N. The transverse scalar coefficient is a nonconstant smooth trigonometric function of time. Its stratum variances sum to Theta(a^2 N^(-3)), with constants uniform at small a. Thus this family has tau_j>=c_j a N_j^(-3/2). The strict gap, fixed-rank constants and every reserve coordinate are actual parameters, not free covariance symbols.

The present conditional C2 estimate for passing that reserve through the remaining M-j force levels is

    physical error <= C_M sqrt(a) a^(M-j) tau_j.       (1)

At a fixed smooth near-kink caller, the first subsequent force has a genuine mean displacement of order a tau_j. Additional affine small-force stages propagate it by their actual a factors. Thus a uniform conditional source bound cannot replace (1) by a tau_j^2 estimate without another theorem.

Combining the target a^R allowance with the gap condition gives

    tau_j <= C a^[R-(M-j+1/2)],
    N_j >= c a^[-(2/3)(R-(M-j+3/2))]_+.              (2)

The bracket in the exponent means the positive part of R-(M-j+3/2). This is exactly the strong earlier-layer allocation already present in the randomized scalar ledger. The nearest previous layer j=M-1 has

    N_(M-1) >= c a^[-(2R-5)/3]

when that exponent is positive. Known logarithmic allowances and fixed-order constants do not change it.

For more than one reserve, the aligned conditional mean errors may have the same sign. A sufficient sum budget therefore does not remove the largest exponent by averaging independent reserve draws. Independence centers the Gaussian input; it does not center its nonlinear force readout.

## 2. Minimal C2 mean-price fixture

Let G be a standard Gaussian, q=0, and

    b(x)=a[x/2+c(sqrt(x^2+eta^2)-eta)],  0<c<1/4.

This has b(0)=0 and a strict Hessian sandwich a(1/2-c)<=b'(x)<=a(1/2+c). It uses a smooth original potential. For eta/tau tending to zero,

    E b(tau G) =a c [E sqrt(tau^2 G^2+eta^2)-eta]
               ~a c tau sqrt(2/pi).                 (3)

Any independent final Gaussian keep leaves this mean error unchanged. Thus a final W2 comparison still pays at least its mean displacement. All records and the zero trajectory are legitimate; no discontinuous source or unavailable derivative is needed.

This is a fixed-caller conditional test. Integrating a nondegenerate, unobserved Gaussian q can improve its marginal law price. That possible gain is precisely a new observer-qualified source claim, not a consequence of (3) or a contradiction to it.

## 3. Why exact reserve cancellation does not retain a free gap

One may evaluate

    b(q+D+tau G)+[b(q+D)-b(q+D+tau G)]=b(q+D).         (4)

This is an exact same-record VALUE identity, with two original queries and their actual firsts. But the observer in (4) reads q+D, either directly or by reconstructing (q+D+tau G)-tau G. Its joint retained record is therefore the noisy node together with the reserve noise (or the original unbuffered node).

For a whole path, stacking (Q+tau G,G) leaves the reconstructible coordinate Q. If Q's reference lies in a low-rank carrier plane, the linear observer that subtracts tau G still sees its transverse zero-Gram direction. The added Gaussian is not an unread reserve for that corrected observer.

Subtracting only a linear HVP response does not remove the conditional mean (3): its Gaussian mean is zero. Differentiating that saved HVP to manufacture a higher correction would introduce the original-third-derivative problem. A genuine finite VALUE counterpacket or a higher-order marginal heat comparison would be a new contribution and must retain its own descendants and callers.

## 4. What stronger reserve theorem would change the arithmetic

If a new source theorem replaced (1) by an observer-qualified price

    C_M sqrt(a) a^(M-j) tau_j^p,  p>1,

with all mixed/root/covariance terms covered, then the formal count would improve to

    N_j >= a^[-(2/(3p))(R-(M-j+1/2)-p)]_+.

This calculation is conditional, not a proved source return. It shows why testing the reserve's **actual law effect with its observers** is the productive gate. A fixed finite p still gives a linear-in-R exponent. Order-dependent p with controlled source closure could do more, but that is the missing all-order correction problem.

The scalar final-only keep lemma in the companion note supplies p=2 for its very specific invariant-posterior boundary. It does not apply to arbitrary intermediate conditional path nodes or to their later nonlinear observers.
