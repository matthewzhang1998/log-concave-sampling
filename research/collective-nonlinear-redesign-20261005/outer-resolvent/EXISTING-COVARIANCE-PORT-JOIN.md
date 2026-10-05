# Exact compatibility check against the existing full covariance reserve

2026-10-05. This checks the actual port in `same-carrier-feedback/FULL-COVARIANCE-CLOSURE-AND-RESUMMED-MEAN-FRONTIER.md` against the new conditional quartic lemma. It does not construct a replacement reserve or silently change its fixed variance budget.

## 1. What already matches

Rename the retained caller Z in the existing service to Y. Its target is

    R_v(Y) approximately N(0,vI+Cov(F3|Y))

with conditional integrated W2 error `Lambda A^4sqrt(D)` plus restored floors, for a fixed positive buffer v. The analytical comparison

    ||Cov(F3|Y)-Cov(F2|Y)||_(L2;HS)<=C A^4sqrt(D)

is already admitted. Thus the named covariance and retained caller agree with the delayed true F2 tail, whose covariance is q^2 Cov(F2|Y). The old mean source also supplies the needed q·m_cont(Y) conditional mean up to order A^3.

No new covariance FUNCTION is required mathematically. The remaining question is the type of actual finite return and its use, not whether the analytical covariance has been identified.

## 2. What does not match the quartic RAW port

The new lemma requests an actual b_fin(Y,W), conditionally independent of X given Y, with full private-bank first O(A), its own conditional mean, and its own conditional covariance within the stated one-energy Hilbert-Schmidt tolerance.

The reserve supplies a completed LAW with a fixed additional vI. It integrates every owned coarse, first-filter, pair, fill, mark, mean and action root. Its carrier-subtracted first is O(Lambda A), but that is not a statement that the carrier-subtracted covariance equals Cov(F2|Y). The law's Gaussian reference noise is explicitly not identified with an actual source-zero root.

This distinction is necessary even in dimension one. Let a desired small covariance be s^2 and let

    R_v=sqrt(v+s^2) G, C_v=sqrt(v) G.

Then R_v has EXACTLY the advertised N(0,v+s^2) law. Its known carrier-subtracted first is `sqrt(v+s^2)-sqrt(v)=O(s^2)` and can be made smaller than an O(A) port. Nevertheless

    Var(R_v-C_v)=(sqrt(v+s^2)-sqrt(v))^2=O(s^4/v),

rather than s^2. Perfect marginal law and a small residual first therefore do not justify subtracting the carrier to manufacture the desired small random force.

Using the whole R_v as b_fin also fails the desired small-bank premise: its Gaussian carrier gives an O(1) private first at fixed v, and its covariance includes vI. The error term in the quartic lemma would then be O(K), not O(A^4).

Nor can one reinterpret a W2 law allowance directly as the required one-energy conditional covariance error. The generic second-moment estimate pays the square-root moment of the buffered law, potentially another sqrt(D). The actual producer/cumulant port must be checked rather than inferred from the law target.

## 3. The internal Gram construction does not already export the missing raw field

The actual rectangular first-coefficient source targets B_Q(X)p through a conditional MEAN calibration. Its completed return M_X(p) targets N(B_Q(X)p,I). Integrating fresh p then supplies the covariance I+B_Q B_Q*. The program never exports B_Q as a strong matrix or B_Qp as a small deterministic action.

The raw partial-variable first filter has a conditional mean close to B_Qp but retains its own private variance. Small mean calibration does not make its actual covariance B_Q B_Q*. The native LAW comparison is exactly what removes that variance from the analytical reference, at the price of the fixed Gaussian carrier. Reading one of its integrated z,H/filter/mean roots afterward would exceed the proved comparison.

The mixed-K reserve uses the same kind of fixed positive-buffer completed-law composition. Its real covariance target is useful, but its hidden complete banks cannot be repurposed as a fresh small RAW conditional force without a new joint/strong source theorem.

## 4. A precise possible LAW-only join, with the new obligations exposed

There is a legitimate direction that uses the complete reserve without subtracting its carrier: pay the reserve's full positive buffer from the actual conditional X Gaussian variance.

At an outer correlation t, let

    X=t z+cG, c^2=1-t^2,
    Y=qX+sigma Z, sigma^2=1-q^2.

For fixed exposed z, the true joint Gaussian disintegration is

    Y|z ~ N(qt z,(1-q^2t^2)I),
    E[X|Y,z]=[t sigma^2 z+q c^2Y]/(1-q^2t^2),
    Var(X|Y,z)=v_X I,
    v_X=c^2 sigma^2/(1-q^2t^2).                           (1)

The true future F2 tail is independent of X conditional on Y and z. Therefore a completed mean/covariance LAW could, in principle, replace the conditional law of the buffered tail while retaining precisely (z,Y), with total positive mean/covariance variance shares summing to v_X. Scaling an unscaled covariance reserve by q requires its corresponding buffer to be v_X/q^2, after accounting for the independent mean branch's share.

Adding the exposed z as an unread exterior observer is compatible with an existing Y-conditional LAW when all complete private banks are freshly drawn conditional on Y and independent of z. That bookkeeping is not itself an obstruction, and the integrated Y marginal is still standard Gaussian. The forbidden move is retaining or reusing a private root that the existing LAW comparison already integrated.

This is a genuine conditional-law composition, not subtraction of a source-zero carrier. However it is NOT the currently admitted fixed-gap join:

- v_X tends to zero as t approaches one. The fixed-v native radius, covariance-gap, precision and current bounds cannot be used unchanged.
- The reserve's covariance-mixture theorem has its explicit factor v^(-3/2). At shrinking v this is an intrinsic current allowance, not an assignable numerical floor.
- The mean and covariance services must use complete independent banks conditional on the SAME exposed (z,Y), with exactly the variance shares in (1). All caller and row chains remain live.
- The true source/current comparison must preserve the short-prefix contribution or its actual nonlinear correction, rather than replace it by an unshifted covariance at x.
- Any endpoint cutoff, missing Gaussian room and retained-variable law must be proved at the actual finite quadrature nodes and charged in the full bill.

The new outer-resolvent quartic lemma instead keeps its dual test tied to X and asks for a small RAW covariance-matched source. It cannot simply consume the reserve's marginal LAW as if those were the same port. The LAW-only disintegration route above belongs to a different, explicitly joint construction and needs its own proof.

## 5. Exact remaining join, rather than a vague missing producer

The active covariance service already supplies the correct conditional covariance as a fixed-buffer positive terminal LAW. To use it here, one must establish ONE of the following concrete conversions:

1. A small RAW force return with the true F2-tail conditional mean/covariance, O(A) complete-bank first, and the full same-caller covariance error. The existing known-carrier subtraction does not establish this.
2. A joint LAW join through (1) that budgets the reserve's carrier from the actual conditional Gaussian variance, proves the shrinking-buffer guards and intrinsic errors, and retains (z,Y) until the final original-g consumption.

Neither conversion follows from the existing reserve's published terminal LAW contract. Thus the covariance producer is not absent in an undifferentiated sense: the exact terminal covariance LAW exists; the RAW or shrinking-buffer JOINT conversion is the still-unproved part. This is the first relevant input-type debt for the quartic outer-resolvent bootstrap.
