# Independent adversarial audit of the nonlinear-shear weak transfer

Publication copy: nonmathematical context and/or local paths were sanitized. Original and public SHA-256 values are recorded in `INVENTORY.json`; source/audit pins below identify their historical versions, not these edited bytes.

## Verdict and audited scope

**Scoped PASS.** The claimed coefficient replacement

    ||O_X-O_X,eta||_(L2_X;HS) <= C kappa eta e/s^2

is valid for the stated shear geometry and the intact conditional marked Jacobian. It requires neither a modulus of the scalar derivative nor a terminal third derivative. With eta_j=b s_j^2, its positive weighted allowance is C b kappa e. The full-square genealogy supplement preserves the needed independent conditional banks and the orientation order/transposes.

I found no counterexample satisfying those precise hypotheses. There are explicit same-source counterexamples if the marked conditional Jacobian is replaced by an unaveraged same-innovation derivative, or if a fast observer of the exposed scalar is silently retained. Those are useful regressions, not objections to the supplied callback graph.

This report does not promote a coefficient comparison to an executed repair. C/R separately own the finite pair, conditional-reference, retained-body, and old-endpoint joins. The inverse-eta source paths and the actual positive-clock width must still be paid there.

New mathematical work; no earlier artifact hash is asserted for this report.

### Source pins

- Candidate initially reviewed at `677340ee75a280f84a90affab7ecc6a27b75cf3de8dbf71254072f44d2612749`.
- The clock qualification identified in this audit is fixed in `p-shear-heat/NONLINEAR-SHEAR-CONDITIONAL-HEAT-CANDIDATE.md`, current SHA256 `81c01db6986b0f86fe7315b07928fbb7a3ad3d8c42487a3d371141bf24de34fd`.
- `p-shear-heat/FULL-SQUARE-SHEAR-SOURCE-GENEALOGY-AND-COST.md`, current SHA256 `a6542cddd67a9ecec004b764291569d431f9b0a0621a2f6d058382f64a3b62f6`.
- `r-true-orientation/POSITIVE-COMMON-INPUT-ORIENTATION-CLOCK.md`, SHA256 `67ca0e5c2df108bdc7a7ddf1f00918274241c3c57b665d691aa3d75a59b68bc4`.

Paths above are relative to `research/`.

## 1. Exact genealogy and orientation order

Let u,v be orthogonal unit vectors, Pi=I-vv*,

    q(W)=W+a u g0(v*W),
    f(W)=r[g(q(W))-g(W)],  kappa=ra,
    |g0'|<=1,  g=grad V,  ||Dg||<=1.

At c^2+s^2=1, define J_f(X)=E_U Df(cX+sU). Because the covariance clock is the diagonal whole-input OU representation,

    O_node=-Sym E_X[(E_U Df(cX+sU))(E_V Curl f(cX+sV))].    (1)

U and V are independent COMPLETE standard innovations conditional on X. This is an identity for a product of conditional means. It does not detach the ancestor from the terminal within either occurrence of f. The factor order is J_f times Curl_f, with the minus sign and final symmetrization as displayed.

Only the V occurrence is split into V_perp+Z_v v. The marked callback is

    f_X(w)=[f(cX+s w)-f(cX)]/s.

Its ancestor is g0(c v*X+s v*w). Its read set contains X,w and the actual old caller; it does not contain Z_v or the exposed z=c v*X+s Z_v. Thus J_f(X) is independent of Z_v conditional on X for a concrete reason in the executed callback graph. The marker's moving origin is also independent of Z_v and remains a complete source occurrence.

Direct differentiation gives

    Curl f(W)=kappa g0'(v*W)
       [H(q(W))u v*-v u*H(q(W))],  H=Dg.

The v-parallel component of H(q)u cancels from the bracket. Hence the exact terminal consumer is the transverse block Pi H Pi, not a separately heated full primitive Hessian.

## 2. Uniform transfer proof, with the derivative and boundary issues exposed

At fixed X,z let

    q_z(Z_perp)=c Pi X+s Z_perp+z v+a u g0(z),
    Hbar_X(z)=E Pi H(q_z(Z_perp))Pi,
    b_X(z)=Hbar_X(z)u.

For unit alpha,beta perpendicular to v, Gaussian integration by parts gives

    alpha*Hbar_X(z)beta
       =s^(-1) E[(alpha*Z_perp) beta*g(q_z(Z_perp))].       (2)

The same transverse Gaussian is used at z and z'. Since g is unit-first and g0 is unit-first,

    ||Hbar_X(z')-Hbar_X(z)||op
       <= (1+a)|z'-z|/s.

Also ||b_X||<=1, b_X is perpendicular to v, and Lip b_X<=(1+a)/s. The estimate tests arbitrary fixed unit vectors before taking their supremum, so no entrywise sum or dimension factor appears. Equation(2) uses only original first derivatives in an analytical Gaussian identity; it does not emit a differentiated HVP source.

Write p for the N(mu,s^2) density, mu=c v*X, and m_eta(z)=E_x g0'(z+eta x). For every fixed unit t, B(z)=t*b_X(z) is bounded by one and Lipschitz by (1+a)/s. The function Bp is absolutely continuous, and

    ||(Bp)'||_1 <= [(1+a)+sqrt(2/pi)]/s.                  (3)

There is no missing boundary term: B is bounded and the Gaussian density vanishes at both infinities. Translation of this L1 density gives

    |E[(m_eta-g0')B]|
       <= eta E|x| ||(Bp)'||_1
       <= C_a eta/s,
    C_a=sqrt(2/pi)[1+a+sqrt(2/pi)]<2  for a<=1/2.          (4)

Only |g0'|<=1 is used. No g0'' or derivative continuity modulus is needed. The proof remains valid by approximation for the stated C1 gradient source. Taking the supremum over t proves the vector version without a dimension loss.

For a vector d perpendicular to v, the skew matrix d v*-v d* has operator norm exactly |d|. Therefore the curl replacement satisfies, uniformly at every X,

    ||C_f(X)-C_f,eta(X)||op <= C_a kappa eta/s.             (5)

The complete marked source supplies the remaining energy factor. Conditional Gaussian first-chaos Bessel gives

    ||J_f||_(L2_X;HS)^2
       <= s^(-2) E_X Var_U[f(cX+sU)]
       <= e^2/s^2,                                      (6)

where e=||f-Ef||_2 is the actual complete-source energy. No energy of either terminal, no dimension radius, and no derivative of a small VALUE error replaces e. Constants in f and its stored anchor vanish from this analytical Bessel pairing; they are not dropped from execution.

Combining(5)--(6) yields the claimed C_a kappa eta e/s^2. The common X is retained throughout, so correlations between the marked conditional Jacobian and the scalar-side conditional coefficient at X are not discarded. The required independence concerns only the newly exposed scalar innovation within the product of conditional means.

## 3. The extra structure that defeats fast-phase objections

The transverse Gaussian average is substantial. A generic bounded matrix field could have an arbitrarily rapid dependence on z, but it need not be the transverse Hessian of a globally unit-first gradient. The integration-by-parts estimate(2) uses both that integrability and the actual common-input transverse width s.

For example, with u=e1,v=e2, take one original uniformly convex potential with gradient

    g(x,z)=(h x, h z+d sin(kz)/k),  h=1/2,d=1/10,
    g0(z)=v*[g(zv)-g(0)]-h z=d sin(kz)/k.

This is a literal same-potential chart. Its complete shear remainder is exactly

    f(W)=kappa h d sin(kW2)/k e1,
    e=kappa h d sqrt((1-exp(-2k^2))/2)/k.

Put alpha_eta=exp(-k^2 eta^2/2), and

    M(X)=kappa h d exp(-k^2s^2/2) cos(kcX2).

Then J_f(X)=M e1e2*, C_f(X)=M(e1e2*-e2e1*), and

    O_X=M(X)^2 e1e1*,  O_X,eta=alpha_eta O_X.              (7)

The L2_X error in(7) is an explicit Gaussian cosine-fourth moment. Retuning k=1/eta does not resurrect the previous raw-surrogate resonance: the intact marked and curl means carry exp(-k^2s^2), which suppresses it when eta=b s^2. The actual energy in(7) remains O(kappa/k), so the test does not hide the defect in an ambient mark.

A second family tests genuine terminal/ancestor frequency beating. Let

    V(x,z)=h(x^2+z^2)/2-d cos(lambda x+kz)/(lambda^2+k^2),
    g=grad V,
    g0(z)=v*[g(zv)-g(0)]-h z
          =d0 sin(kz)/k,
    d0=d k^2/(lambda^2+k^2).

The original Hessian lies in [(h-d)I,(h+d)I] uniformly. The actual f is odd and

    (h-d)kappa ||g0||_2 <= e <= (h+d)kappa ||g0||_2.

The transverse consumer is exactly

    b_X(z)=[h+d lambda^2/(lambda^2+k^2)
                 exp(-lambda^2s^2/2)
                 cos(lambda cX1+kz+lambda a g0(z))]e1.

Large ancestor frequency k is accompanied either by transverse heat damping or by the integrability factor lambda^2/(lambda^2+k^2). The displayed consumer cannot behave like an arbitrary unit-amplitude cos(kz) matrix when k greatly exceeds 1/s.

The adjacent checker uses the exact Jacobi--Anger Fourier representation of this entire f, retains its actual J_f, and computes conditional orientation and actual energy from Gaussian Fourier moments. It tests 108 frequency/width combinations, including k eta of order one. The omitted Bessel coefficient tail is bounded below 1e-40; the numerical diagnostics do not substitute for the uniform proof in section2. The maximum tested error divided by kappa eta e/s^2 was approximately 0.00166.

## 4. Hostile regressions for forbidden ancestry changes

The same simple source in section3 produces a concrete failure if one replaces the marked conditional mean J_f(X) by raw Df(cX+sZ) and reuses the exposed Z in the curl factor. At k eta=1, the wrong product sees a cos^2 phase rather than two independently averaged conditional means. Its averaged error is

    (kappa h d)^2 (1-exp(-1/2))
                    [1+exp(-2k^2)]/2.

For s=.25, eta=.001, k=1000, its ratio to kappa eta e/s^2 is approximately **869.45**. This is a genuine same-source, one-energy failure of the altered genealogy. It is not a counterexample to(1), whose independent conditional innovations are explicitly supplied by the new callback graph.

Likewise an additional retained observer cos(kz) defeats a claim based only on |B| and Lip b. Even with B=1,

    E[(m_eta-g0')cos(kz)]
       =-(1-exp(-k^2eta^2/2))E cos^2(kz)

for g0'=cos(kz). At k eta=1 and ks large its magnitude approaches 0.1967 rather than O(eta/s). The observer adds its own z derivative to(3). Thus a retained body reading Z_v must first be replaced through its admitted conditional stationary/reference law, or its observer current must be paid. The scalar coefficient calculation alone cannot perform that join.

## 5. Clock and inverse-width audit

For the positive orientation clock, eta_j=b s_j^2 gives exactly

    sum_j w_j C_a kappa eta_j e/s_j^2
        =C_a b kappa e sum_j w_j <= C b kappa e.

The executed complete side-root row is different. With side source active radius rho and actual captured-label row rho/eta, the fork readout kappa/(b rho^2) and full clock readout sqrt(w_j) give

    (kappa/(b rho)) sum_j sqrt(w_j)/eta_j
      =(kappa/(b^2 rho)) sum_j sqrt(w_j)/s_j^2
      <= Lambda kappa/(b^2 rho s_min).                   (8)

The dyadic positive panels prove the last sum by a geometric inverse-width sum, with the declared sqrt(M) logarithmic factor from per-panel Cauchy--Schwarz. The individual inequality w_j=O(s_j^2) alone would not establish a dimension/node-count-free version of(8).

The initial wording delta<=b was insufficient to infer s_min comparable to b. This audit reported that issue to P. Current `81c01db6` explicitly chooses delta=c0 b and retains the general kappa/(b^2 rho delta) price if a smaller unrelated delta is used. With that correction, (8) gives the claimed Lambda kappa/(b^3 rho). The exponent condition 3 zeta+gamma<g is the correct ordinary-A guard for kappa=A^(1+g), b=A^zeta, rho=A^gamma.

Five independently generated positive-clock rules, at tolerances .2 through .003, pass the stated mass/width identities and actual sqrt-weight path checks. Choosing a much smaller clock tolerance increases the path sum and cannot be hidden as scalar precision.

An eta-scaled source call itself needs finitely many VALUES, not eta^(-1) Gaussian replicas:

    (rho/eta)[g0(z+eta y)-g0(z)].

Its active first stays O(rho), but its captured-z first can be O(rho/eta). The full-square supplement charges this path, the real minimum query width, tighter VALUE/row tolerances, anchor replay, and log(1/eta). There is no inverse-eta sample multiplicity forced by this displayed source or the admitted fixed-order pair interface. This audit does not replace an expanded compiler census: the actual inherited pair count and numerical source profiles must still be supplied by that interface, at the chosen finite order. A matrix coefficient bound is not a cost certificate for an alternative empirical implementation.

## 6. Full-square algebra and remaining implementation boundary

For the supplement's source, the selected block matrix is

    D=Hbar+m_eta vv*,   J=u v*-v u*,
    K=c_sel^2 rho^2 D J D.

Because Hbar v=0 and u is perpendicular to v,

    D J D=m_eta[Hbar u v*-v u*Hbar],  K*=-K.

The known middle fill Pi0=I-uu*-vv* satisfies JJ*+Pi0 Pi0*=I. With the marked selected matrix M=c_sel J_f and O=-kappa Sym(J_f D J D),

    M K*+K M*=2 c_sel^3 rho^2 O/kappa.

Thus b c_side=-q kappa/(2 c_sel^3 rho^2) produces the intended negative covariance correction. The factor order and transposes are correct; no forward-order substitution or omitted two-curl term is being used in this algebra. Random nonsymmetric marked matrices through dimension32 were used as independent sign/order regressions. Those algebra tests do not assert positivity for arbitrary unscaled random test matrices.

The source is a genuine full-n gradient on two fixed orthogonal blocks. Neither a hidden unknown basis nor a Hessian-valued source is needed. The actual original VALUE implementation and its first sweep retain both anchors; original HVPs belong only to that first/adjoint sweep.

The remaining distinction is substantive: after the finite pair comparisons, the exposed coefficient roots and the retained unmarked body still require the exact conditional-reference/keep chronology. In particular, using the favorable orientation-only clock is legitimate only with the independent old constant-reference join stated in the candidate. Regenerating the old forward/Cov clock at the same nodes changes its accuracy requirement and the width bill. This audit certifies the transfer and its required weighting, not that separate endpoint theorem.

## Reproduction

    python research/audit-shear-transfer/check_shear_transfer.py

All **625 regressions pass**: 108 nonlinear Fourier cases, 27 exact-formula terminal cases, five actual positive clocks, and full-square/transposition tests. The invalid genealogy and observer fixtures are deliberately recorded as failures of altered hypotheses. New artifact hashes are in the local manifest.
