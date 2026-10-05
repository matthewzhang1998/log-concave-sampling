# Independent audit: discounted OU cumulants and one-derivative clocks

Date: 2026-10-04.

## Verdict and pinned scope

**PASS for the corrected analytical connected-cumulant hierarchy and the one-derivative positive-clock lemma.** The generating PDE, averaging-symmetrization convention, binomial factors, covariance and third-cumulant specializations, and dimension-safe DC estimate are correct. The clock derivative result is a Gaussian Hilbert-operator statement and does not turn any coefficient or derivative into a VALUE oracle.

Audited final files:

1. `../DISCOUNTED-OU-CONNECTED-CUMULANT-RECURRENCE.md`, SHA256 `cd37a782bbcef3cdf5ec39a33be4874df56ff06af75170485b182c04b361d05e`.
2. `../FIRST-RESOLVENT-CLOCK-QUADRATURE.md`, SHA256 `56968a4ff1880fee1dff80cb6b8f0044c27ee8108dca4c0d8a3ec8c45ee08926`.

One local error was found in the first inspected recurrence pin `ba98b38420dce8f83a3d75861b570e0e65bd81849344791302fcb8d996e388bf`: Section 3 put a full Sym on only the right-hand side of an unsymmetrized martingale covariance identity. The author removed that inner Sym, leaving the required outer cubic Sym intact. The final pinned display was checked. The final third-cumulant formula is unaffected.

The independent checker passes **5,024 assertions**. It compares raw-moment and connected recurrences through scalar rank six, checks coupled two-output tensor factors, derives the quadratic third cumulant directly from the true conditioned Markov covariance, and tests actual positive dyadic rules. No author diagnostics were imported.

## 1. Discounted Markov generating equation

Write `H_z=int_0^infinity e^(-s)g(X_s^z)ds` for the standard OU process with generator `L=Delta-z·grad`. Splitting at t and applying the Markov property gives

    H_z=H_[0,t]+e^(-t) H_future(X_t),
    M(z,lambda)=E_z[exp(lambda·H_[0,t])M(X_t,e^(-t)lambda)].

Expanding at t=0 yields

    L_z M-lambda·partial_lambda M+(lambda·g)M=0.

Since the generator has diffusion coefficient one, `L exp(K)/exp(K)=L K+|grad K|^2`, with no factor one half. Therefore

    (L_z-lambda·partial_lambda)K+|grad_z K|^2+lambda·g=0.

The signs and diffusion normalization in the source are correct. In correlation notation r=e^(-t), the conditional OU path has mean rz and the same covariance as the original path ending at X_1=z. The integral is the true conditional Markov functional; a single common-noise substitute has another law.

For a globally Lipschitz g, the functional has finite Gaussian moments of every fixed order and finite exponential moments in each scalar direction. Its mean grows at most linearly in z. Finite-dimensional/path smoothing justifies the local generator calculation. Differentiating in lambda at zero uses only finite conditional moments; spatial identities can then be interpreted weakly. No third derivative of the underlying potential is needed for this analytical construction.

## 2. Tensor coefficients and the positive resolvent

Use `K=sum_k kappa_k:lambda^k/k!`. The product of the degree-i and degree-(k-i) spatial derivatives in `|grad K|^2` has coefficient `1/[i!(k-i)!]`. Multiplying by k! gives the binomial coefficient. Because lambda's output slots are symmetric, the tensor coefficient is standard averaging Sym, not an unnormalized sum. Thus

    (1-L)kappa_1=g,
    (k-L)kappa_k=sum_(i=1)^(k-1) binom(k,i)
                         Sym(Dkappa_i·Dkappa_(k-i)).

The contracted index is the shared spatial derivative index. Neither a physical output contraction nor a product of uncentered raw moments would give this formula. The use of log M automatically performs the connected subtractions.

The positive Laplace resolvent is

    (k-L)^(-1)F=int_0^infinity e^(-kt)P_(e^-t)F dt
                =int_0^1 q^(k-1)P_qF dq.

On Gaussian chaos n its multiplier is 1/(k+n). Polynomial/moment growth, or the original Markov representation, selects the stated solution; arbitrary homogeneous solutions are not part of the claim.

Let v=kappa_1, B=Dv and C=kappa_2. At k=2 there is one term of coefficient 2, so

    C=2R_2(BB*).

For gradient g, derivative/OU commutation gives `B=int_0^1 q P_q Dg dq`, which is symmetric. This specializes to exactly the previously audited `2 int qP_q(B^2)dq`. For a general nongradient source it must remain BB*, as the source specifies.

At k=3, the i=1 and i=2 terms both have coefficient 3. Their full averaging symmetrizations agree after reordering output slots. Consequently

    kappa_3=6R_3 Sym(B·DC).

No matrix commutation has been used. B and DC are evaluated at the same coarse state inside P_q; independent evaluations at unrelated centers would change this product target.

## 3. Independent martingale factor check and the corrected display

The conditional martingale

    M_t=int_0^t e^(-s)g(X_s)ds+e^(-t)v(X_t)

has differential `sqrt(2)e^(-t)B(X_t)dW_t`, since `(1-L)v=g`. Its bracket recovers C. Applying Itô to the centered third tensor gives

    kappa_3=6 Sym int_0^infinity e^(-2t)
                     E[(M_t-v) tensor BB*(X_t)]dt.

The factor 6 is three product-rule placements times the bracket's factor two. For an arbitrary matrix field A, the exact unsymmetrized covariance identity is

    E[(M_t-v) tensor A(X_t)]
      =2 int_0^t e^(-s)P_time-s[B·D P_time-(t-s) A] ds.

Full symmetrization may be applied to both sides when this is inserted in the cubic formula, but may not appear on only one side. A simple witness is g(x)=x in dimension two and `A(x)=x_1 e_2 tensor e_2`. The left tensor has component (1,2,2) equal to `t e^(-t)` and the other permutations zero, while full averaging Sym changes that component to one third as large. This confirms that the initial display required the correction now present in the audited file.

Set t=s+v_time, q=e^(-s), u=e^(-v_time). The old measure is `e^(-2t)e^(-s) ds dt`; it becomes `q^2 u dq du`. Combining the two factors 6 and 2 gives

    12 Sym int_0^1 q^2 dq int_0^1 u du
                       P_q[B·D P_u(BB*)].

Since `DC=2 int_0^1 u D P_u(BB*)du`, this equals `6R_3 Sym(B·DC)`. The time/correlation notation switch and every factor are therefore consistent. Rough bounded matrix fields can first be smoothed; the identity is an analytical covariance statement, not permission to execute a derivative-valued producer.

As a separate test independent of that derivation, take the scalar algebraic source g(x)=x^2-1. The exact recurrence gives

    kappa_1(z)=(z^2-1)/3,
    kappa_2(z)=2(z^2+1)/9,
    kappa_3(z)=16(3z^2+2)/135.

At z=0, Wick's formula on the true conditional Markov path gives covariance `2/9` and third cumulant `32/135` directly. On the ordered domain r<=s<=t, the cubic covariance product is `r^2(1-s^2)(1-t^2)^2/t^2`; integrating it with coefficient 48 gives 32/135. This quadratic fixture checks algebra and path identity, not the theorem's global Lipschitz hypotheses.

## 4. The dimension-safe covariance first

With the whole standard Gaussian driving tape, H has private first at most A, hence Gaussian Poincare gives `Cov(H|z)<=A^2I`. Its literal z derivative obeys

    ||D_z H||op<=int_0^infinity e^(-t) A e^(-t)dt=A/2.

For a unit u, differentiating covariance uses only this first:

    D_u C=Cov(D_u H,H)+Cov(H,D_u H).

The matrix covariance inequality

    ||Cov(U,V)||HS<=||U-EU||2 sqrt(||Cov(V)||op)

follows by applying the scalar-L2-to-vector map `a -> E[a(V-EV)]` row by row. Centering is L2 contractive, so `||D_uH-E D_uH||2<=A/2`, with no extra factor two. The two covariance terms yield

    ||D_u C||HS<=2 A(A/2)=A^2.

Thus `DC:R^D->HS` has operator norm at most A^2, and its full tensor HS norm is at most A^2sqrt(D). Also `||B||op<=A/2`. Contraction of B with DC therefore costs at most `(A^3/2)sqrt(D)`, and Sym is contractive. R_3 has total positive mass 1/3, so multiplication by 6 gives

    ||kappa_3(z)||HS<=A^3sqrt(D).

The estimate is uniform in z and uses one physical sqrt(D) energy. These constants are conservative; no sharper bound is required. It does not expose DC or kappa_3 as an executable field oracle.

## 5. Positive-clock one-derivative addendum

The ordinary clock estimate controls the multiplier errors

    e_n=Q(r^(n+k-1))-1/(n+k), |e_n|<=delta0.

Positivity and the actual dyadic geometry supply an additional high-chaos tail. On a panel at distance [a,2a] from one, all nodes satisfy r<=1-a and the total weight is a, so the full panel contributes at most `a e^(-ma)` to Q(r^m). This bound is independent of the number and arrangement of positive interior Gauss nodes.

For m>=1, set x=ma. The inequality

    x e^(-x)<=2 int_(x/2)^x e^(-t)dt

and disjoint dyadic intervals bound the panel sum by 2/m. The terminal midpoint contributes `h e^(-mh/2)<=2/(em)<1/m`. Therefore `Q(r^m)<=3/m<=6/(m+1)`. At m=0 the mass is one, so C=6 works for all m>=0.

Since k>=1, positivity and the exact coefficient give

    |e_n|<=min(delta0,7/(n+1)).

Gaussian Hilbert-valued chaos isometry gives the exact derivative norm

    ||D(R_k-Q)F||2^2=sum_(n>=1) n|e_n|^2 ||F_n||2^2.

The two scalar error bounds imply `|e_n|^2<=7delta0/(n+1)`, hence `n|e_n|^2<=7delta0`. This proves the displayed dimension-free derivative error `sqrt(7delta0)||F||2`. Asking the existing positive quadrature for delta0<=delta^2/7 therefore yields a delta first-derivative bound with the same polylogarithmic node class. No HVP-valued source instruction is needed to state this analytical result.

The stated extension to Sobolev multiplier `(n+1)^a`, 0<=a<1, follows by interpolating the same two bounds. The restriction is real: with finitely many nodes strictly below one, the quadrature multiplier decays exponentially as n grows, whereas the resolvent multiplier is 1/(n+k). Thus the twice-differentiated error has asymptotic norm one on arbitrarily high chaos. One cannot infer a uniform second-derivative approximation from this rule and plain Gaussian L2 data.

Finally `B=Dv=R_2(Dg)` follows from derivative/OU commutation. DB is one derivative of that resolvent acting on the bounded original Hessian field Dg, so the analytical addendum applies with its ordinary Gaussian Hilbert L2 norm. It does not authorize reading Dg or DB as a producer value, and does not by itself control products of independently approximated tensor fields.

## 6. Diagnostic and construction limits

`check_discounted_ou_cumulants.py` writes `discounted_ou_cumulant_checks.json` and passes 5,024 assertions. The source pins are enforced. Raw moments are independently generated from the linear MGF equation `(k-L)mu_k=k g mu_(k-1)` and converted to cumulants; these agree with the nonlinear hierarchy through rank six. A coupled two-output polynomial gradient verifies the actual tensor factors in k2 and k3. Independent ordered Markov Wick integrals agree with the quadratic covariance and third cumulant. The clock diagnostics use real positive interior Gauss nodes, verify the mass/tail/error bounds over a range of chaos degrees, and exhibit the high-chaos two-derivative obstruction. These computations supplement the proof and are not all-order executable producer tests.

The hierarchy concerns the exact discounted additive H only. Higher posterior reference orders also contain nonlinear resolvent/Picard feedback. An analytical connected recurrence and a scalar positive quadrature theorem do not price those feedbacks, construct tensor-free VALUE actions, or prove positive full-law replacement. The source correctly leaves the first/caller/root/product/numerical and same-endpoint law gates open.

**Conclusion:** the corrected recurrence and its one-derivative clock addendum are valid analytical targets at the inspected pins. They supply no finite all-order posterior program and no eventual-sublinear query exponent theorem.
