# Reconstructed matched-source calibration obstruction

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

## Provenance and status

This is a NEW reconstruction. It is not claimed to reproduce the original bytes or hashes of either:

- the independent audit `850ca174f7d44ab181ffa6718f6ac4ecbe4dc4733c6288962faac04cd300d99d`;
- C's source `238acc863fed07c0d4e08723cff51995609573f65a9bca056887fc7639b5dd94`.

The former audit had 1,471 high-precision regressions. Those old execution artifacts are not restored by this note. A new independent checker adjacent to this file passes **542 checks** and reproduces the same certified sign constant and numerical limit. 

The mathematical result is preserved: a literal two-node nonlinear curl VALUE surrogate has a correct zero jet and kappa-small VALUE envelope, yet its weak calibration can have an order-kappa-times-actual-energy defect with an intact finite first-response input. Every deterministic whole-input scale sigma(A)<=1 can be defeated uniformly over the bounded-Hessian source class. This is a constructor/interface result, not an identified output of a named weak-proxy remainder factory.

## 1. Literal source and one original strongly convex potential

Fix c=d=1/10, epsilon=A^(39/10), nu=epsilon/A, any deterministic sigma in(0,1], and k=1/(sigma epsilon). Set

    Q(x)=x/2+d sin(kx)/k,
    g0(x)=Q(x)-x/2=c sin(kx)/k,
    gt(x)=Q(x)-(1/2-nu)x=nu x+d sin(kx)/k,
    q(x)=x+A g0(x),
    f(x)=A gt(q(x)),  kappa=A^2.

Q is the SAME original gradient in all calls, with Hessian in[.4,.6]. The adapters use only known linear subtraction and original VALUES. Their own potentials need not be convex. f is scalar, odd, and source-zero at zero; its true curl is identically zero.

The exact bound

    sup_x |f(x)-epsilon x| <= sigma epsilon A(d+epsilon c)

gives actual centered Gaussian energy e_f/epsilon ->1, uniformly over sigma<=1. Also Lip f<=A(nu+d)(1+Ac). The normalization therefore uses the actual tiny source mark, not an unrelated first radius or dimension estimate.

This scalar graph has P=C=1, so its nonterminal reads the physical terminal coordinate. It does not satisfy the separate protected-terminal nonterminal-independence requirement. Original-oracle realizability alone is not named weak-proxy admission.

## 2. Exact finite response, with the normalization retained

Use the literal common-node first-response definition

    I_f(p,z)=sum_j d_j [f(sqrt(1-b_j^2)z+b_j p)
                       -f(sqrt(1-b_j^2)z-b_j p)]/(2b_j).

The finite Chebyshev--Lobatto filter has sum d_j=1, sum |d_j|<=4, and b_j^(-1)<=2K. All nodes use the SAME z, and all f calls use the same full finite source. With S_K=sum |d_j|/b_j<=8K,

    |I_f(p,z)-epsilon p|
       <=S_K sigma epsilon A(d+epsilon c).                (1)

This is unnormalized physical I_f, equivalently A I_(f/A). Inserting I_(f/A) without the factor A would be a different source. The formula is p-odd. One f query uses two original Q VALUES, so its K-node response costs4K VALUES.

## 3. The tested finite curl surrogate

At the saved root W, for input u define

    D0=g0(W+u)-g0(W),
    Dt=gt(q(W)+u)-gt(q(W)),
    C_W(u)=A[gt(q(W)+u+A D0)-gt(q(W)+u)]
                -A^2[g0(W+Dt)-g0(W)].

All saved centers and repeated VALUES are cached identically. The program uses six distinct original Q queries, including its base values. Its zero input and input derivative at zero vanish exactly, as the scalar curl does. The general small VALUE envelope is |C_W(u)|<=2kappa|u|.

For u=sigma v, theta=kW, and w=v/epsilon, the exact normalized expression is

    R(theta,w)=C_W(sigma v)/(sigma kappa epsilon)
      =nu c[sin(theta+w)-sin(theta)]
       +(d/A)[sin(theta+w+Ac sin(theta+w))
                   -sin(theta+w+Ac sin(theta))]
       -c[sin(theta+nu w+d{sin(theta+Ac sin(theta)+w)
                           -sin(theta+Ac sin(theta))})-sin(theta)].

A special stronger input-first estimate holds for this explicit source:

    |partial_w R|<=2nu c+4cd.                             (2)

To verify it, the first term costs nu c. Differentiating the divided sine difference gives a cosine difference bounded by2Ac, plus an inner-sine term bounded byAc; after multiplication by d/A the cost is3cd. The last term costs c(nu+d). No generic small-first theorem is inferred from this trigonometric identity.

Thus the actual rescaled map v -> C_W(sigma v)/sigma is kappa(2nu c+4cd)-Lipschitz. Combining(1)--(2), uniformly in EVERY W,z,p,

    |C_W(sigma I_f)/sigma-C_W(sigma epsilon p)/sigma|
       /(kappa epsilon)
       <=(2nu c+4cd)S_K sigma A(d+epsilon c).              (3)

This is o(1) whenever KA->0, uniformly over all sigma<=1, including fixed sigma. W and z may coincide; the new checker deliberately uses W=z. Only p is fresh and independent of W in the phase limit below.

## 4. Strict nonzero first Gaussian coefficient

Since k>=1/epsilon tends to infinity, theta=kW tends to uniform phase. Its nonzero wrapped-Gaussian Fourier coefficients are exp(-m^2 k^2/2). The exact R above is uniformly bounded, and dominated convergence gives its limiting phase mean at w=p:

    F(p)=c[(d/2)sin p-J1(2d sin(p/2))cos(p/2)],
    J1(t)=E_phi[cos(phi) sin(t cos(phi))],

where phi is uniform on[0,2pi]. This definition suffices; no external special-function identity is needed.

Taylor's remainder under the phase integral gives

    J1(t)=t/2-t^3/16+R5(t),  |R5(t)|<=|t|^5/384.

The linear term in F cancels. Using sin^3(p/2)cos(p/2)=sin p/4-sin(2p)/8,

    E[p F(p)]
       >=(c d^3/8)[exp(-1/2)-exp(-2)]
             -(c d^5/12)sqrt(2/pi)
       =0.000005823451825883353700049... >0.               (4)

The numerical integral is0.000005884780025222952937554..., inside this analytic error interval. The strict sign does not rely on quadrature.

Equation(3), multiplied by |p| and integrated, proves the SAME nonzero limit with the actual finite I_f:

    E[p C_W(sigma I_f(p,z))/sigma]/(kappa e_f)
          ->Lambda_small>5.82e-6.

It follows by Cauchy--Schwarz that the conditional mean's L2_p norm is at least a fixed multiple of kappa e_f. The correct scalar curl target is zero. Therefore this held source cannot have a uniform O(A kappa e_f) calibration error. Oddification in the p-odd input leaves this first Gaussian coefficient unchanged.

The quantifiers are essential: for each prescribed sigma(A), the original potential is retuned through k=1/(sigma epsilon). This refutes a uniform bounded-Hessian-class rate from whole-input shrinking; it does not assert failure of sigma->0 for each fixed, bounded-frequency potential.

## 5. Boundaries retained from the original reports

- The finite VALUE envelope and the zero jet are valid; the failed claim is stronger weak calibration.
- The incoming mark is the literal actual finite first response of the SAME f, at unit width. No independent generic Gaussian is substituted in the final result.
- The theorem does not address an unspecified fully integrated Riesz kernel with different clock weights.
- A changed internally attenuated or variance-preserving circuit requires its own analysis.
- The source uses one strongly convex original gradient via affine adapters, but is not identified with a named protected weak-proxy E batch.
- This result does not refute the newer nonlinear-shear conditional-heat transfer: that construction averages a transverse Hessian against an exposed scalar Gaussian with an intact conditional marked Jacobian, a different genealogy and calibration mechanism.

## New reproduction

    python research/audit-matched-curl-reconstructed/check_reconstructed_matched_curl.py

The new checker runs literal original-Q callbacks, common-node filters, same-root surrogate VALUES, strong-input-first tests, the transfer bound, and the analytic sign constant. Its542 checks pass. New hashes in the adjacent manifest establish the current reconstructed artifacts only.
