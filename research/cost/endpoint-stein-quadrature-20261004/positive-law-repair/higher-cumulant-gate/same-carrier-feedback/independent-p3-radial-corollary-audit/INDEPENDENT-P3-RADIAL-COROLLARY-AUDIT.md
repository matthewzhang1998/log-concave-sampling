# Independent audit: radial obstruction to the unshifted P3 mean candidate

Date: 2026-10-05 UTC.

## Verdict

**PASS, with the public-log scope qualification in Section 6.** The frozen corollary extends the independently audited single-history separator to the actual two-level force inside the three-root raw P3 source. It refutes equation (9) of the frozen P3 raw-source note: the unshifted first-displacement covariance correction has an eventual L2 mean discrepancy of order A² at A=D^(-1/2), against an advertised allowance of Lambda A³.

No mathematical blocker remains for this extension. It does not refute the raw F3_Q source, its guarded OWN-mean compiler, a centered or exactly shifted/resummed current, or every possible fourth-order sampler. No full endpoint-law conclusion follows solely from this conditional-mean separator.

Primary reviewed text, left unchanged:

- `../RADIAL-COROLLARY-FOR-THE-UNSHIFTED-P3-MEAN-CANDIDATE.md`
- SHA256 `bae3fff77a2cfe75e0888637b4d39a1346bb64c1eec1ef8217f3e151af1eecc0`

The canonical radial proof, its independent audit, the P3 raw note, the earlier nested-mean theorem and independent review, the P3 component review, and the positive-quadrature source were read directly. All eight exact input pins appear in `MANIFEST.json` and the checker. The new algebra and estimates below were checked independently; the canonical smooth radial estimates and covariance witness are admitted at their specified reviewed pins.

## 1. Correct imported mean and clock tolerances

The old nested theorem's F_cont, conditional on its exposed Gaussian endpoint, is precisely F2_H here. Its finite F_Q becomes the displayed F2_Q after renaming its outer root G to H and its inner root H to J. It retains the same H across all middle nodes and the same J across all innermost occurrences. No Markov/shared-root law equality is used.

The actual intrinsic theorem is

    ||E F2_H - E F2_Q||2
       <= C[A³ + delta_mid A + delta_in A²] sqrt(D).

The source writes one delta dominating both tolerances; separating them as above follows from its inner and outer estimates. With delta_in,delta_mid<=A² and sqrt(D)=1/A, this is at most C(2A²+A³). Its constant is absolute, so the corollary's weaker Lambda A² allowance is valid. This is an L2 estimate under x~gamma_D, not a pointwise-in-x bound or a strong retained-root coupling. First moments alone would not suffice. The corollary expressly requires the missing actual Hermite multiplier tolerances.

## 2. Exact decomposition and conditional moment profile

Put c=d=1/4 and g=c A id+d A f, with |f|<=1. The exact nonlinear coefficients are

    K_H = -c integral e^(-t) I_t dt
          +d integral e^(-t) f(X_t-I_t)dt,
    K_Q = -c sum_j v_j L_j +d sum_j v_j f(y_j-L_j).

Thus both decompositions have A K_a, not d A K_a. The leading terms are c A(x/2+sigma_a G_a), with sigma_H=1/2 and sigma_Q=beta=sum_j v_j sqrt(1-tau_j²). The true linear OU integral has mean x/2 and conditional covariance I/4. The finite leading term follows from the middle rule's exact first moment. Each G_a is standard and independent of x, while K_a may depend on that same Gaussian root and all its ancestors.

For fixed p, positivity and Minkowski give the common conditional bound

    ||K_a||_(Lp|x) <= c A(|x|/4+C_p sqrt(D))+d
                   <= C_p(1+A|x|).

For H, use the marginal of X_(t+u) conditional on x in the double positive integral. For Q, each z_jk has mean sigma_k tau_j x and conditional Gaussian covariance at most I; both exact first moments give the factor 1/4 after averaging. This uses only marginals and positivity, not independence between nodes. Integrating x~gamma_D yields ||K_a||p<=C_p when A sqrt(D)=1.

The leading Gaussian conditional means agree exactly. Therefore

    E F2_H-E F2_Q=A(E K_H-E K_Q),
    ||E K_H-E K_Q||2<=C A.

The division by A is essential and is performed correctly in the corollary. No pointwise boundedness of K_H is claimed or needed.

## 3. Taylor legality and correlations

The specific radial family is globally smooth: its perturbation vanishes on the ball r<=sqrt(D)-2. Its second derivative as a bilinear vector map is uniformly bounded. This follows directly from the radial derivative formula in the canonical proof, with all nonzero denominators satisfying r>=sqrt(D)-2>=sqrt(D)/2. Taylor is being used for this fixed smooth admissible family, not for every member of the general C2 class.

Let y0=(1-cA/2)x and y_sigma=y0-cA sigma G. The imported radial estimates hold for every fixed p, in particular p=4:

    ||Df(y_sigma)-Df(y0)||_(L4;op)<=C A,
    E f(y_sigma)=f(y0)+(c² sigma²/2)A psi'(S_D-c/2)n_x
                 +O_L2(A²).

Expand only in the additional vector A K_a. The integral Taylor remainder is bounded by C A²|K_a|², whose joint L2 norm is at most C A²||K_a||4². Conditional expectation contracts L2. For the correlated linear term, joint Hölder gives

    ||E[(Df(y_sigma)-Df(y0))K_a | x]||2
       <= ||Df(y_sigma)-Df(y0)||4 ||K_a||4 <=C A.

This works on the enlarged true-history probability space as well as the two-root finite space. It makes no independence assertion about K_a and G_a. Consequently

    E f(y_sigma-A K_a)
       =E f(y_sigma)-A Df(y0)E K_a+O_L2(A²).

Since Df(y0) is uniformly bounded and the K-mean mismatch is O(A), its difference contributes O(A³) after the outer d A. The outer linear component contributes -c A²(E K_H-E K_Q), also O(A³). Every new nested-ancestor contribution is therefore below the A² leading radial bias.

## 4. Coefficient, covariance, and resolvent witness

Combining these estimates gives, with an absolute remainder constant,

    psi2_H-psi2_Q=k_Q A² psi'(S_D-c/2)n_x+O_L2(A³),
    k_Q=(d c²/2)(1/4-beta²).

The correction in raw equation (9) deliberately uses Cov(F1_H|x)-Cov(F1_Q|x), not the full second-force covariance. It is exactly the covariance pair in the canonical radial proof. That proof controls its actual nonlinear covariance remainder by Gaussian first-chaos Bessel and the HS-to-vector Hessian contraction bound. It gives

    (C_H-C_Q):D²g(x)/2=k_Q A² psi'(S_D)n_x+O_L2(A³).

Hence no new covariance approximation is required in this extension. Degree-two multiplier accuracy implies beta>=2/3-delta_mid>=7/12, so |k_Q|>=13/18432, uniformly in the admitted finite middle rules.

Set h(s)=psi'(s-c/2)-psi'(s) and T_D(x)=A x. Then ||T_D||2=1 and R1 T_D=T_D/2. As D tends to infinity,

    <R1[h(S_D)n_x],T_D>
       =(1/2)E[h(S_D)|x|/sqrt(D)] -> (1/2)E h(S),
    S~N(0,1/2), E h(S)<0.

The strict sign follows by integrating the symmetric bump's interval superlevel sets: every nonzero translation decreases a centered Gaussian interval's mass. Boundedness of h and |x|/sqrt(D)->1 in L2 justify the limit. Contractivity of R1 controls the O(A³) remainder. The true P3 candidate discrepancy is therefore at least c_* A² for all sufficiently large D. This is the additional full-m3 comparison that the earlier single-history audit correctly did not claim.

## 5. Actual top quadrature and OWN-mean compiler

Under independent standard x,H,J, every unshifted affine site y_j,z_jk has standard Gaussian marginal. Positivity gives

    ||F2_Q||2<=A(1+A)sqrt(D),
    ||psi2_Q||2<=A(1+A+A²)sqrt(D).

The raw F3_Q uses one frozen middle/inner source at every outer node. Thus its conditional mean is exactly Q_out psi2_Q, even though its private roots are shared across nodes. The actual top multiplier certificate delta_out<=A³ gives

    ||(Q_out-R1)psi2_Q||2<=C delta_out A sqrt(D)<=C A³.

The imported, guarded raw OWN-mean compiler has integrated conditional error at most Lambda A⁴ sqrt(D)=Lambda A³, plus its actual separately paid floors. An ideal exact unshifted counterterm still leaves the larger A² conditional-mean discrepancy. A compiler within the declared allowance of that wrong mean-law target cannot repair it: conditional W2 dominates conditional mean distance, and triangle inequality preserves the separator. This conclusion concerns the stated retained-Gaussian-carrier target and does not assert an unconditional endpoint-law impossibility.

## 6. Necessary public-log reading

For the advertised fixed-log refutation, choose an explicit admitted sequence

    A=D^(-1/2), delta_in=delta_mid=A², delta_out=A³,

with ordinary polynomial-accuracy numerical floors and versions satisfying all inherited guards. The positive dyadic rules have O(log²(1/A)) nodes; their actual errors may be smaller than their declared tolerances. On this sequence every public logarithm is O(polylog(1/A)), so Lambda A tends to zero and Lambda A³ cannot absorb c_* A². The checker additionally uses a harmless stricter inner tolerance A²/2 to test unequal rules.

This is the correct reading of the corollary's sentence about fixed public-log polynomials. It is not a claim that deliberately exponentially over-resolved auxiliary versions have logarithmic logs. A uniform theorem is refuted by this one admitted bounded-version sequence. The analytical A² lower bound itself is stronger: using the intrinsic old estimate with absolute C makes it uniform over all rules with the stated inner/middle tolerances. No frozen input was edited to record this clarification.

## 7. Reproducible checks and limits

Run `python check_p3_radial_corollary.py`. Its 200 passing assertions verify eight exact pins; positive finite clock examples and their analytical all-degree certificates; selected Hermite multipliers; the actual two-level same-root nonlinear decomposition with distinct inner and middle rules; quadratic mean normalization; true OU linear normalization; exponent bookkeeping; and the scalar limiting witness.

The largest literal decomposition residual is below 9.8e-17. The diagnostic values are E h(S)=-0.00182950693274536 and its linear-resolvent half=-0.000914753466372681. These floating-point values supplement the analytical strict-sign proof. The checker does not simulate a true OU path, certify a finite-D asymptotic threshold, prove the imported analytical theorems by computation, or construct a centered/resummed consumer.

Only this independent audit directory was written. `MANIFEST.json` pins sources and artifacts; `SHA256SUMS` seals the completed packet.
