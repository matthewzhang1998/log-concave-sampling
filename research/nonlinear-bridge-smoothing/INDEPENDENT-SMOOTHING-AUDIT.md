# Independent audit of the C2 stage-two smoothing transfer

2026-10-05. Status: **PASS, with the stated guarded-import and dimension qualifications.**

The main theorem has a valid analysis-only route to grade 8/3 for the unchanged finite original-VALUE graph. The entropy orientation, nonsymmetric Jacobians, conditional Gaussian inputs, genuine history, finite outer replacement, original caller captures, and execution bill have been checked independently. No blocking defect was found. This is not an independent implementation or reproof of the imported native own-mean compilers, and it does not establish a full endpoint join or arbitrary-order closure.

## 1. Audited inputs and exact conclusion

The sealed stage-two manifest matches the supplied SHA256

    cd1111e28dc399643f27bcdd0f97dc735ada196f8a43f220211500bcb45be3ae.

Its main theorem SHA256 is

    c7275de9261bce2eeaee1ba054351055b1e236590340ee7b1554ff064ffa2f13.

The new files reviewed are:

- C2-STAGE-TWO-SMOOTHING-TRANSFER.md, the main entry;
- BASELINE-SMOOTHING-TRANSFER.md, including the finite-stencil extension in Sections 7–8;
- FINITE-SMOOTHING-IMPLEMENTATION.md, an optional implementation that is not required by the main theorem.

The valid main conclusion is a guarded finite canonical-m3 mean-law component for anchored g=grad U, U in C2, 0<=Dg<=A I, and 0<A<=1/36. Under a declared orthogonal block bound b, the new target-bias certificate is at most

    110 A^(8/3) b^(2/3) sqrt(D),

plus exactly the displayed inner-rule, finite-outer, completion, and restored numerical terms. The source shape may vary with A; no uniform Hessian modulus is required. The block decomposition remains an assumption, not something inferred from finitely many VALUES. Its basis need not be disclosed because the executed graph is unchanged. For a completely general D-dimensional source, b=D is valid and gives 110 A^(8/3) D^(7/6), not a dimension-uniform sqrt(D) result.

## 2. Entropy: correct orientation and no symmetry gap

For X~N(m,c^2 I), c>0, write T(x)=x-F(x), ||DF||<=l<1. The contraction equation x=y+F(x) gives a global inverse. The matrices I-s DF stay invertible for 0<=s<=1; their determinants therefore remain positive. No symmetry assumption on DF is needed.

The relevant relative entropy is KL(T#Law(X) || Law(X)), in that order. Change of variables gives

    KL = E[-(X-m).F/c^2 + |F|^2/(2c^2) - log det(I-DF)].

Gaussian integration by parts converts the first term into -E tr DF. The real log-determinant branch continued from I satisfies

    -tr K-log det(I-K) = sum_(n>=2) tr(K^n)/n
                         <= d l^2/[2(1-l)].

This is valid for nonsymmetric real K: |tr(K^n)|<=d||K||^n suffices. Positivity of each individual trace is neither assumed nor needed. Pinsker, with total variation converted to the L1-density convention, then gives exactly

    |E f(T(X))-E f(X)|
      <= ||f||_infinity sqrt(E|F|^2/c^2+d l^2/(1-l)).

There is no missing factor two. The vector-valued version follows by duality or by integrating the Euclidean norm against the signed density. The proof differentiates F once only. It does not differentiate DF or an inverse-Jacobian coefficient, and therefore does not insert a hidden Hessian-Lipschitz assumption.

## 3. Genuine-history target transfer

The conditional OU path is represented as X_s=e^-s x+xi_s, with the whole process xi independent of x. Its actual future dependence is retained in both F1,a and F2. The independent derivative calculation gives

    ||D_x F1,a|| <= A e^-a/2,
    ||D_x F2|| <= A/2+A^2/4 = L.

The stationary marginal and anchoring give ||F2||_2<=A(1+A)sqrt(d). For two coherent sources with sup|g-h|<=delta, their first-history terms differ by at most delta, and their second-history terms by at most (1+A)delta. These estimates preserve all nested arguments.

The remainder decomposition

    [g(x-F2,g)-g(x)]-[h(x-F2,h)-h(x)]
      = (g-h)(x-F2,g)-(g-h)(x)
        +h(x-F2,g)-h(x-F2,h)

is exact. The second line's last term is already at most A(1+A)delta. To treat the first term, fix the caller z, outer time t, and entire private future xi, then apply the entropy estimate to x=tz+sqrt(1-t^2)G. Only after that conditional estimate are xi averaged and the L2 norm in z taken. Conditional Jensen and Minkowski justify this order.

The only endpoint singularity is 1/sqrt(1-t^2), whose integral is pi/2. No uniform bound as t approaches one is claimed; the point t=1 has zero Lebesgue measure. This yields the main entry's exact C_d(A), and C_d(A)<5sqrt(d) on A<=1/2. The target transfer therefore has the required extra A, while replacing the full m3 target would leave an unsuppressed baseline term.

## 4. Finite-stencil weak stability and baseline cancellation

This is the decisive additional step. The original finite residual is F_Q(q)-q(x). After expanding its terminal source VALUES, the terminal signed mass is one and absolute mass is 7/2. The four families have absolute masses 1,1,1,1/2: the base, single centered, linear-pair centered, and polarized quadratic terms.

For fixed private roots, each original U or V is affine in x with x coefficient of magnitude at most one. The exact noncommuting derivative of e is

    De=[Dq(U)-Dq(U-q(V))]U_x
       +Dq(U-q(V))Dq(V)V_x
       -Dq(q(w))Dq(w)/4.

The two symmetric matrices in its bracket are in [0,A I], so their difference has norm at most A, not 2A. The resulting terminal derivative-defect bounds are

    S:                  A/2+A^2/4,
    S +/- d:            2A+A^2/4,
    S +/- e:            3A/2+3A^2/2,
    S +/- d1 +/- d2:    7A/2+A^2/4.

All are bounded by L_*=(7/2)A+A^2/4<1 when A<=1/36. The generally nonsymmetric nested products are retained in their original order.

For stationary x and private roots, put k=sqrt(3/8), s0=1/sqrt(2)+Ak, d0=1+1/sqrt(2), and e0=1+k. Then

    ||S-x||_2<=A s0 sqrt(d),
    ||d||_2<=A d0 sqrt(d),
    ||e||_2<=A^2 e0 sqrt(d).

Weighting each terminal family, rather than multiplying one worst displacement by 7/2, gives exactly

    H_sum=(7/2)s0+2d0+A e0
         =2+11/(2sqrt(2))+A[1+(9/2)k].

The source-difference identity is

    E_g-E_h = sum c_i[(g-h)(T_i,g)-(g-h)(x)]
                 +sum c_i[h(T_i,g)-h(T_i,h)].

The baseline term cancels because sum c_i=1. Coherent terminal-argument differences satisfy

    |S_g-S_h|<=(2+A)delta,
    |d_g-d_h|<=2delta,
    |e_g-e_h|<=(3+2A)delta.

Summing their A-Lipschitz terminal costs gives exactly A delta[14+(11/2)A]. The first sum is controlled by conditional entropy terminal by terminal, using the weighted H_sum and total Jacobian mass 7/2. This proves the stated H_d(A). At A=1/36, H_1=36.486920877..., and the exact rational upper certificate

    4647443/127008 < 37

follows from pi<22/7, sqrt(2)>7/5, k<5/8, L_*<1/10, and sqrt(1-L_*)>18/19. Thus H_d(A)<37sqrt(d) is safely justified.

These estimates require identical source-independent finite rules on both sides. They are uniform in rule length and node positions, so rules may depend on A and the declared b. Shared roots across distinct nodes are never asserted to form a single genuine OU history. Only each node's true marginal/conditional row law is used.

## 5. Analytical smoothing and the exponent

The anchored Gaussian smoothing h_e is a proof object. It has the same anchored PSD Jacobian interval and preserves the original orthogonal block decomposition. On a d-block,

    sup|g-h_e|<=2 A e sqrt(d).

Centering Dg by A I/2 in the Gaussian derivative kernels gives the dimension-free induced-multilinear bounds

    Lip(Dh_e)<=A/(sqrt(2pi)e),
    Lip(D2h_e)<=A/(sqrt(2)e^2).

The second kernel has variance 1+(u.v)^2<=2. These are bounds for the ideal smooth source only; no regularity is imputed to a finite sum of shifts of a C1 source.

Writing J_q=R1(psi_Q,q-q) and R_q=m3(q)-R1q, the exact telescoping identity is

    J_g-R_g=(J_g-J_he)+(J_he-R_he)+(R_he-R_g).

Both outside terms have the weak extra-A estimate. The middle term is the sealed smooth stage-two bias applied to h_e on the same rules. Squared block norms sum correctly: replacing d_j by b in constants leaves sqrt(sum d_j)=sqrt(D), not an extra ambient dimension factor per block.

The exact resulting terms are those in main equation (4.3). Choosing e=A^(2/3)b^(1/6) gives A^(8/3)b^(2/3)sqrt(D). No restriction e<=1 is needed. The p_A and q_A coefficients remain fully priced; at A=1/36 they are approximately 1.311961379 and 2.978850902. Even the conservative upper calculation

    84+1.32sqrt(3)+3sqrt(15) = 97.905257105... <98

is below the declared 110. Rational upper checks are also available: use sqrt(2)>707/500, k<613/1000, and sqrt(2pi)>1253/500 to get p_A<=167467602683/127565424000<1.32 and q_A<=208417693252713926614471/69941736046820736000000<3. The last radical bound follows from the classical pi>223/71.

The coarse cap is valid independently: ||J_g||_2<=A^2[3/2+5/(2sqrt(2))+A(1+2k)]sqrt(D) and ||R_g||_2<=A^2(1+A)sqrt(D). Their sum is less than 5A^2sqrt(D). One may take the minimum of this cap and the refined exact-outer bound, then add the finite-outer, completion, and numerical allowances.

## 6. Outer norm, native reentry, and execution accounting

The new comparison is an exact-outer mean statement. It does not imply the same rate for the unaveraged inner graph, a raw sample, or a single outer quadrature node. The actual finite outer rule is introduced only afterward through its sealed L2 Gaussian operator bound, applied to the actual unsmoothed psi_Q,g. Its norm bound A C_F2 sqrt(D) gives precisely delta_out A C_F2 sqrt(D).

The actual original source and every original stage-two port are unchanged. The sealed first/curl/energy/caller proofs require the anchored PSD Hessian interval, not B or C. Thus the guarded own-mean completion remains applicable with private dimension 5D, normalized first <=9A, curl <25A^2, the original O(A^2) energy, actual caller-origin capture, independent complete baseline/residual banks, and all native guards. The output-law comparison is a comparison of completed own means with the canonical mean, followed by the conditional Gaussian W2 identity; raw roots or source records are not appended to the observed output.

No smoothed VALUE, smoothed HVP, new smoothing root, or hidden convolution call is executed. The original occurrence count, source parameter dependencies, live original caller anchors, one-HVP-per-recorded-site sweeps, full replay rules, and original numerical floors therefore remain unchanged. At a nonzero caller, its actual captured origin is still restored. At the original all-zero caller/root input the actual original graph has its original literal zero descendants; at zero source the imported Gaussian zero-source row remains. The proof parameter e does not create an executed division by a vanishing energy or source parameter.

The weak estimate concerns a coherent change from one exact source to another everywhere in the graph. It does not improve arbitrary independent VALUE/HVP/row/weight errors. Those keep the sealed absolute floors, including finite-precision moment and mass defects. This distinction is explicitly preserved in the main text.

The conditional normalization corollary has also been checked. For f(y)=s[g(a+s y)-g(a)], Df(y)=s^2 Dg(a+s y), so its anchored gradient interval is [0,alpha I], alpha=As^2, and the same blocks are preserved. For 0<alpha<=1/36 the theorem applies with alpha and all actual new guards. The alpha=0 case uses the literal zero-source branch instead of the positive-radius proof. Both original source sites, live anchor/scale firsts, full-key reuse, replays, and physical output rescaling remain explicitly charged. This is structural closure, not an endpoint induction.

## 7. Pair-count and optional finite smoothing checks

The stated unrestricted-dimension pair count is valid for 0<delta<=1. Since log C_D=O(D), the logarithm of the internal inverse tolerance is O(D+log(1/delta)). The sealed quadrature degree has that size, and its dyadic-panel count is O(log(1/epsilon)+log log(1/epsilon)), hence also O(D+log(1/delta)). Marginal count is therefore O((D+log(1/delta))^2), and tensor pair count O((D+log(1/delta))^4). The very small fixed holomorphy radius produces very large absolute constants. This is a source-node count, not a practical-runtime or dimension-free-arithmetic assertion. Native occurrences and coefficient setup/precision are still charged separately.

The optional finite smoothing note also passes its relevant checks. The clipped product grid's centered independent errors have scalar covariance bounded by tau^2 I, which is sufficient for every unknown block projection; rotational invariance of a finite law is not claimed. The rational positive-weight repair preserves mass and symmetry. Every finite shifted source remains anchored, convex-gradient, and only C1 in general. Its smooth bias is transferred from the ideal source, rather than assigning ideal derivative bounds to the actual finite graph. Unknown block structure does not lower the generic tensor exponent D; disclosed block coordinates can use the common b-dimensional rule with correct block marginals. Original shifted anchor leaves may be nonzero even when a macro source value is zero, and the note correctly keeps their dependencies and charges.

The main result does not need this optional grid. Its one-node interpretation is exactly g(x+e*0)-g(e*0)=g(x), or directly the telescoping comparison above. Consequently a vanishing-grid-tolerance bill cannot be imposed on the analysis-only theorem, and an exact smoothing oracle has not been concealed.

## 8. Independent diagnostics and limitations

Run:

    python independent_smoothing_checks.py

The script presently passes 23,511 assertions and records its result in independent_smoothing_checks.json. It checks the supplied sealed-manifest pin; exact affine-Gaussian KL formulas against the transport formula for nonsymmetric matrices; determinant orientation and the log-series upper bound; terminal-map derivatives against finite differences on genuinely noncommuting convex-gradient examples; every declared terminal derivative/displacement bound; signed/absolute terminal masses; weighted coherent-source displacement; the rational H_d constant; p_A/q_A and interpolation exponents. The maximum affine identity discrepancy is about 2.3e-13, and maximum terminal finite-difference discrepancy is about 1.9e-10.

The separate author check_smoothing_constants.py was also inspected. Its 14,522 assertions include exact rational bounds H_d/sqrt(d)<37 and 19587/200<98 for the conservative main coefficient. Its parameter checks are consistent with the main entry and do not claim a numerical target-law proof.

These are diagnostics, not a replacement for the analytic arguments above. No Monte Carlo rate experiment, executable native compiler, high-dimensional dimension-free result, general-C2 order-four theorem, or all-order recurrence is claimed. All conclusions remain conditional on the same imported finite positive own-mean services and their actual native guards.
