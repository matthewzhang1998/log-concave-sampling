# Independent analytical audit: a shrinking-buffer, LAW-only join

2026-10-05. This audit checks Gaussian geometry, direct conditional weak comparisons, retained variables, integration, finite outer quadrature, the variable-buffer native mean/covariance extension, actual final source ports, and complete finite count recurrences. Native claims remain subject to their literal imported guards; this audit does not numerically execute the entire LOW30 compiler.

## Verdict

The proposed analytical route is sound. A genuine conditional LAW join can use the complete mean/covariance services without subtracting their Gaussian carriers, provided their positive variance shares sum to the actual conditional variance v_t below. A bulk/near-endpoint split with h=w=A^(4/5) gives analytical target bias O(A^(16/5) sqrt(D)); the unmatched-covariance near branch costs only O(A^(17/5) sqrt(D)). No small RAW covariance-matched source is needed for this route.

The initial analytical obligations were actual variable-buffer service guards/errors/costs and native ports for the new composite finite source. Section 8 below and the full-draft review verify their guarded extensions and finite count recurrences; ordinary fixed-buffer contracts alone would not have established them. The finite outer rule must be applied to the fixed exact delayed target first, followed by nodewise replacements; the t-dependent surrogate is not P_t applied to a common inner source.

## Final review of the executable draft

**PASS under the declared imported native guards.** The draft `../LAW-ONLY-SHRINKING-BUFFER-JOIN.md` was checked in full, including Sections 9–11. Its proof closes this specific guarded LAW-only canonical-m3 mean join. It does not establish endpoint accuracy, an all-order recurrence, or a numerical execution of the whole native compiler.

The reviewed draft's SHA256 at the diagnostic run is `40891db38221e3b02f2e657a4ef93c3e4fd06c69b29cc03e7db9c5761fe2fc0a`. The independent diagnostic `check_law_only_join.py` passes **18,374 assertions**; `law_only_join_checks.json` records scope, the original LOW30 pin, and the tested draft hash. Tests cover exact Gaussian geometry, common-carrier/branch rotations, native scaling substitutions, noncommuting ordered K covariance, positive endpoint quadrature, square-root shield-weight estimates, and full root-count algebra. They do not replace the proof or execute imported LOW30. Float64 quadrature tests are intentionally restricted to A>=10^(-4); smaller A requires the separately declared increased scalar precision, because fixed machine precision can round strictly interior endpoint nodes to one.

No substantive correction to the reviewed draft was required. The only requested wording clarification distinguishes the independent, unused analytical eta in the true-target comparison from a service's actual reconstructed source-zero carrier after the common-row rotation. For the draft's plus-carrier convention an explicit rotation is

    U=(q c/sqrt(s))G+(sigma/sqrt(s))Z,
    eta=(sigma/sqrt(s))G-(q c/sqrt(s))Z,
    a_t(z,Y)+sqrt(v_t)eta=t z+cG,
    Y=q(t z+cG)+sigma Z.

Both coefficients are bounded by one. The complete service banks stay independent of (z,Y), and the separate branch banks stay conditionally independent of one another. They need not be independent of the reconstructed X after its Gaussian room has been consumed. The Sections 9–11 carrier pinning, final raw-split ports, count formulas, and guard powers A^(3/5), A^(7/10), A^(11/15) are correct.

The draft retains the logarithmic term Lambda A^(19/5) separately and absorbs it only under Lambda A^(3/5)<=1. This is essential and correct; a log-dependent constant must not silently become numerical.

## 1. Exact conditional geometry

Let z,G,Z be independent standard D-Gaussians and, for fixed t<1, define

    c^2=1-t^2, X=t z+cG,
    Y=qX+sigma Z, sigma^2=1-q^2,
    d_t=1-q^2t^2=c^2 q^2+sigma^2.

Conditional on z,

    Y=q t z+sqrt(d_t)U,
    X=a_t(z,Y)+sqrt(v_t) eta,
    a_t(z,Y)=[t sigma^2 z+q c^2Y]/d_t,
    v_t=c^2 sigma^2/d_t.

Here U and eta can be taken independent standard Gaussians, and eta is independent of (z,Y). In executed Gaussian-row form,

    a_t=t z+(q c^2/sqrt(d_t))U.

The true future force b(Y,W) is independent of X and z conditional on Y, with fresh complete future bank W. Thus eta is also independent of W conditional on (z,Y). The true disintegration of X-b is the conditional law of

    a_t(z,Y)+sqrt(v_t)eta-b(Y,W).

For every t the unconditional Y marginal is standard Gaussian. This fact permits exactly the same integrated-in-Y service errors to be used after z is retained.

## 2. Direct pointwise one-energy Gaussianization

Let b=b(Y,W) have full-bank Lipschitz constant L uniformly in Y, and let

    m=E[b|Y], u=b-m, Sigma=Cov(b|Y).

Freeze (z,Y). The conditional centered Stein fields in the already audited quartic lemma satisfy

    ||tau-Sigma||_(2,HS|Y) <= L e(Y),
    ||T3||_(2,HS|Y) <= L^2 e(Y),
    e(Y)=||u||_(2|Y) <= L sqrt(D).

The last inequality uses the physical-rank bound ||D_W b||_HS^2<=D L^2, not the private-bank dimension. For an independent standard H, interpolate

    C_lambda=lambda u+sqrt(1-lambda^2)Sigma^(1/2)H.

For a K-Lipschitz vector F, the covariance-preserving Stein identity gives

    d/dlambda E F(a+sqrt(v)eta-m-C_lambda)
        =-lambda^2 E[T3:D^3F(a+sqrt(v)eta-m-C_lambda)].

All coefficients in T3 are independent of eta after (z,Y) and the private sources are fixed. Integrating two derivatives against eta gives

    T3:D^3F  ->  v^(-1) DF [T3:(eta eta^T-I)].

The contracted object inside DF is a D-vector. Conditional Gaussian isometry yields

    E_eta |T3:(eta eta^T-I)|^2 <= 2 ||T3||_HS^2.

It follows, without a vector-output dimension loss, that

    |E F(a+sqrt(v)eta-b)
      -E F(a-m+sqrt(v)eta-Sigma^(1/2)H)|
       <= (sqrt(2)/3) K L^3 sqrt(D)/v.                  (2.1)

This is pointwise in (z,Y), uniformly in a. It spends the actual conditional Gaussian buffer, not any derivative of an outer resolvent test. There is no need to bound D^2b or a derivative of the Stein field. Euclidean/Gaussian Sobolev smoothing followed by strong limits justifies C1/Lipschitz sources, degenerate Sigma, and the nonexecuted true infinite future bank.

## 3. Direct mean-only near-endpoint comparison

The first centered Stein identity gives

    d/dlambda E F(a+sqrt(v)eta-m-lambda u)
        =lambda E[tau:D^2F(a+sqrt(v)eta-m-lambda u)].

One integration by parts in eta, followed by

    E_eta|tau eta|^2=||tau||_HS^2,
    ||tau||_(2,HS|Y)<=L e(Y)<=L^2 sqrt(D),

proves

    |E F(a+sqrt(v)eta-b)-E F(a+sqrt(v)eta-m)|
       <= K L^2 sqrt(D)/(2sqrt(v)).                    (3.1)

For two sources b1,b2 with full-bank firsts L1,L2 and conditional means m1,m2,

    |E F(a+sqrt(v)eta-b1)-E F(a+sqrt(v)eta-b2)|
       <= K|m1(Y)-m2(Y)|
          +K(L1^2+L2^2)sqrt(D)/(2sqrt(v)).             (3.2)

Conditional averaging over Y|z and Jensen give the same L2_z bound with the mean discrepancy replaced by its L2_gamma(Y) norm. Hence an actual RAW source with mean error O(A^3 sqrt(D)) and complete-bank first O(A) suffices in this branch, despite having unmatched conditional covariance.

## 4. How the completed LAW services may be consumed

Suppose complete fresh outputs M(Y) and C(Y) are independent conditional on Y and jointly compared, with all their private banks integrated, against

    M_ref |Y ~ N(q m_cont(Y), v_m I),
    C_ref |Y ~ N(0, v_c I+q^2 Sigma_3(Y)),
    v_m+v_c=v_t.

Then the only terminal use is

    F(a_t(z,Y)-M(Y)-C(Y)).

The reference conditional input law is

    N(a_t-q m_cont(Y), v_t I+q^2 Sigma_3(Y)).

It is the Gaussianized target of Section 2, up to the admitted Sigma_3 versus Sigma_2 covariance mismatch. Scaling a service by q scales both its covariance and its buffer by q^2, so unscaled service buffers must be budgeted accordingly.

If each service has retained-Y integrated conditional W2 error epsilon_j, the terminal L2_z mean error is at most K sum_j epsilon_j. A coupling conditional on Y can be extended by drawing z|Y independently of all its coupled service roots. This retains precisely (z,Y); it does not recover or retain any private root that a LAW comparison integrated. No pointwise-in-Y error bound is necessary if the service contract is the L2_Y integrated conditional W2 contract.

An L2_Y Hilbert-Schmidt covariance discrepancy DeltaSigma is separately bounded after the positive buffer by

    W2(N(0,vI+Sigma1),N(0,vI+Sigma2))
       <= ||Sigma1-Sigma2||_HS/(2sqrt(v)).              (4.1)

This follows from the common Gaussian square-root coupling and the Sylvester identity for the two PSD square roots. Thus an admitted O(A^4 sqrt(D)) force-covariance mismatch adds O(A^5/sqrt(v_t)) after K=O(A), smaller than the proposed covariance-mixture allowance.

## 5. Exact singularity integrals

The simplifying identity is

    1/v_t=q^2/sigma^2+1/(1-t^2).                       (5.1)

Consequently

    integral_0^(1-h) dt/v_t
      =q^2(1-h)/sigma^2+(1/2)log((2-h)/h),            (5.2)

    integral_(1-h)^1 dt/sqrt(v_t)
      <=q h/sigma+arccos(1-h)
      <=h/sqrt(w)+2sqrt(h),                          (5.3)

where w=1-q and sigma^2=w(2-w)>=w. Also

    integral_0^(1-h) v_t^(-3/2)dt
      <=sqrt(2)[q^3(1-h)/sigma^3
                 +(1-h)/sqrt(h(2-h))]
      <=sqrt(2)[w^(-3/2)+h^(-1/2)].                  (5.4)

Finally the full integral of v_t^(-1/2) is at most q/sigma+pi/2.

For K=O(A), L=O(A), these imply

- bulk non-Gaussian remainder: O(A^4[w^(-1)+log(2/h)]sqrt(D));
- near raw mean-only remainder: O(A^3[h/sqrt(w)+sqrt(h)]sqrt(D));
- near raw mean mismatch: O(h A^4 sqrt(D));
- bulk LAW covariance-mixture error epsilon=O(A^4 v_t^(-3/2)sqrt(D)): O(A^5[w^(-3/2)+h^(-1/2)]sqrt(D));
- Sigma_3 versus Sigma_2 covariance mismatch: O(A^5[w^(-1/2)+1]sqrt(D));
- bulk force-mean discrepancy O(A^3 sqrt(D)): O(A^4 sqrt(D)).

The actual F2 short-prefix replacement and coherent move into F(u)=g(u-wg(u)) retain their independent O((A^2 w^(3/2)+A^3w)sqrt(D)) allowance.

At h=w=A^(4/5), the leading prefix and bulk Gaussianization terms are A^(16/5). The near raw branch is A^(17/5), the covariance-mixture term is A^(19/5), and A^4 log(2/w) is strictly smaller than A^(16/5) uniformly on the stated small-A interval.

## 6. Finite quadrature is compatible, with an explicit comparison order

Let psi be the fixed exact delayed, prefix-reduced inner target and write m_delayed=R1 psi. Keep the independent prefix allowance ||m3-m_delayed||_2=O((A^2 w^(3/2)+A^3w)sqrt(D)). Take the already certified positive dyadic-endpoint OU rule Q with operator tolerance epsilon_out=O(A^(11/5)), exact mass one, and strictly interior nodes. First write

    ||m_delayed-Q psi||_2 <= epsilon_out ||psi||_2
                     =O(A^(16/5)sqrt(D)).

After retaining the separately charged coherent short-prefix comparison, disintegrate the exact P_(t_i) psi at each actual node t_i and replace it using the bulk or near branch according to whether 1-t_i>=h. Positivity then sums the nodewise bounds. Do not apply ||Q-R1|| to the t-dependent surrogate, which need not have the form P_t mu for a common mu.

On a dyadic panel 1-t in [d,2d], the quantity 1/(1-t^2) changes by a factor at most two. Thus v_t^(-s) changes by at most 2^s for s>0, and any positive mass-correct panel quadrature is bounded by 2^s times its panel integral. If h cuts a panel, the bulk union lies in 1-t>=h/2 and the near union lies in 1-t<=2h. These constant enlargements preserve all rates in Section 5.

Choose the final midpoint-panel width a<=min(h,epsilon_out/4), as the existing dyadic rule does for the stated parameters. Its near contribution is bounded directly by

    a/sqrt(v_(1-a/2)) <= q a/sigma+sqrt(2a),

so the endpoint panel does not create a hidden divergence. No node at t=1 is necessary. The number of nodes remains O(log^2(1/A)); scalar arithmetic and node/weight precision remain separately charged.

## 7. Native obligations and their audited resolution

The analytical join initially required the following actual-producer and final-source checks. They are now verified under the imported guards in Section 8 and in the full-draft review of main Sections 9–11:

1. Mean and covariance services available at their declared positive shares of v_t, with all guards checked at every bulk node.
2. Shrinking-buffer intrinsic errors, numerical floors, precision, and native occurrence counts fully restored; fixed-v constants cannot be frozen.
3. Fresh complete independent service banks conditional on the same retained Y, with z only an exterior observer.
4. Actual first/curl/caller fields of the composite source if it is subsequently completed by a native own-mean compiler; law closeness alone gives no derivative-port guarantee.
5. Preservation of the original nonlinear short prefix and coherent source descendants.

With those guarded specializations and full count recurrences verified, the bulk/near construction is admitted as the claimed guarded finite canonical-m3 mean theorem and optional own-mean LAW. Actual numerical admission of all guards at a chosen A, D, precision and source remains an execution-time requirement; no universal fixed-A numerical admission or full native runtime test is claimed.

## 8. Independent check of the variable-buffer native extension

The following checks use the original pinned LOW30 source at
`/workspace/shared/v9-curation-work/frozen/prerequisites/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex`, especially `b27:raw:split`, `b27:compiler:mean`, `b27:compiler:pair`, and `b27:compiler:paths`, together with the admitted finite three-C0 reserve and rectangular first-coefficient/Gram source.

### 8.1 Scaled gradient/near-gradient mean service

Normalize a physical variance-v mean service by 1/sqrt(v). For the existing actual gradient-plus-near-gradient source, take

    r=rho=Lambda A/sqrt(v), delta=a=Lambda A,
    mu=A, fixed gradient order b=4.

The LOW30 raw split error, multiplied by the physical readout sqrt(v), has the following powers (public logarithms suppressed):

    A^4 v^(-3/2), A^4 v^(-1),
    A^4 v^(-1/2), A^4 v^(-1/2),
    A^(9/2) v^(-3/2), A^5 v^(-3/2).

For 0<A,v<=1 these are bounded by Lambda A^4 v^(-3/2) times the source's one-energy/caller profile. Its independently proved actual residual first, after the same readout, is

    Lambda [A+A^(3/2)/sqrt(v)].                       (8.1)

For v>=c A^(4/5), this is Lambda A with the fixed share constant c restored. The LOW30 proof explicitly applies the physical-path argument to external-center derivatives and fixed zero-tape source constants. The small-radius guard must nevertheless be checked at rho=Lambda A/sqrt(v), and all finite clocks, filters, precision and source versions remain live.

The rectangular Gram source has the needed literal gradient/near-gradient split and actual first/curl/energy ports. Its private-H anchors must stay inside each complete source occurrence; only actual retained-(X,p) origins may be captured externally. Scaling all those recorded sources and origins coherently gives exactly the foregoing bounds. Adding independent complete mean banks at their positive shares is legitimate after (X,p) is retained, and p is integrated only after the completed conditional law comparison.

### 8.2 Rebalanced three-C0 mixed-covariance reserve

Let the existing four-clock weights satisfy w_j>0 and W=sum_j w_j=1/16. Set

    nu_j=v_path w_j/W,
    r=(4 W A^3/v_path)^(1/3),
    rho_Y=rho_V=r, rho_U=-r,
    c_j=a_j=sqrt(nu_j/2).

Then

    nu_j rho_U rho_V rho_Y=-4 w_j A^3,               (8.2)

so the exact ideal node covariance is unchanged:

    nu_j I-4w_j A^3 Sym(H_U H_V H_Y).

All three actual normalized source radii are kappa0 r times the imported fixed logarithmic factors. No inverse individual w_j enters a source radius. The native chronological-envelope bound gives

    sum_j c_j Lambda [r^b+r^(b+1)+r^(b+2)]sqrt(D)
      <=Lambda sqrt(v_path) r^b sqrt(D),             (8.3)

where sum sqrt(w_j/W) is explicitly public-logarithmic. At b=5 this is Lambda A^5 v_path^(-7/6)sqrt(D), smaller than the proposed common Lambda A^4 v^(-3/2) law allowance.

The source-path theorem gives the actual residual/private/coarse/retained-caller first bound

    Lambda sqrt(v_path) r
       sum_j sqrt(w_j/W)(1+sigma_j^(-1))(1+r+r^2).
                                                               (8.4)

The square-root weighted shield sum really is public-logarithmic. For one endpoint dyadic panel of width d containing n Gauss nodes, Cauchy-Schwarz gives sum_panel sqrt(beta_i)<=sqrt(n d). On that panel c_i>=c sqrt(d), so sum_panel sqrt(beta_i)/c_i<=C sqrt(n). The terminal midpoint has beta~d and c~sqrt(d), hence contributes O(1). Therefore the unshielded one-clock sqrt-weight sum is O(sqrt(n)), and its shielded version is O(J sqrt(n)) for J panels. Since

    sigma_j^(-1)<=sqrt(10)(c_q^(-1)+c_r^(-1)+c_t^(-1)),

factorization of sqrt(w_j)=sqrt(beta_q q)...sqrt(beta_t t) bounds the four-clock sum in (8.4) by a fixed public-log polynomial. The missing shield on the s clock causes no problem. Thus (8.4) is

    Lambda A v_path^(1/6),                            (8.5)

and is at most Lambda A for v_path<=1. Its absolute path norm back to original g VALUES similarly has factor sqrt(v_path)r/A=v_path^(1/6), with the recorded inverse-shield/filter precision still charged.

The covariance mixture and clock-restoration errors are unchanged as covariance quantities because (8.2) leaves the target and coarse matrix Delta unchanged. Their law prices at shrinking v are respectively O(A^6 v^(-3/2)sqrt(D)) and O(delta_clock A^3 v^(-1/2)sqrt(D)). The former is smaller than the Gram mixture's A^4 v^(-3/2). Positive-gap guards include A^3/4<=v/4 and all actual pair radii; v comparable to A^(4/5) leaves ample algebraic margin, subject to the imported logarithmic guards.

Every three-C0 node retains its leaf passive until the final readout, replaces its complete chronological path, then integrates its complete owned coarse bank. The known zero carrier remains c_j G_(root,j)+a_j p_j, independent of the retained caller; child carriers have zero source-zero transmission because each strict ancestor contributes its residual factor. No LAW reference Gaussian is substituted for those actual carriers.

### 8.3 Pin the actual carriers before applying a final near-gradient mean compiler

This detail is necessary if the new finite source is to have a common outer-gradient baseline. At each bulk node, let the complete mean/covariance returns have independent known actual affine source-zero carriers of variances v_m,v_c, with v_m+v_c=v_t. Their combined actual carrier can be pinned to sqrt(v_t)eta by an orthogonal change of their entire Gaussian input bank; split the two branches by an auxiliary independent Gaussian so they remain independent conditionally on (z,Y). All perpendicular directions stay on the tape.

Let beta_t=q c^2/sqrt(d_t). Since beta_t^2+v_t=c^2, choose

    U=(beta_t/c)G+(sqrt(v_t)/c)Z,
    eta=-(sqrt(v_t)/c)G+(beta_t/c)Z.

Then U and eta are independent standard roots and

    a_t(z,Y)-sqrt(v_t)eta=t z+cG.

Indeed the Y row becomes exactly

    Y=q t z+q cG+sigma Z=qX+sigma Z.

Thus the actual completed-law input is the common baseline X=t z+cG minus the actual combined service residual. This carrier pinning is a deterministic known Gaussian-bank rotation; it does not subtract a carrier to assert a new force covariance and does not identify any comparison noise with an executed root.

Each service bank is independent of (z,Y), and the mean/covariance service banks are mutually independent conditional on (z,Y). They are intentionally not independent of the reconstructed X conditional on Y: their source-zero carrier consumes its conditional variance. The whole LAW comparison retains only (z,Y), exactly as required.

Reusing the same common G across outer nodes now gives the genuine-gradient baseline sum_i omega_i g(t_i z+c_i G). If the actual combined service residual has full first Lambda A and the recorded small energy/origin profile, the terminal original-g difference has first Lambda A, square-lift curl Lambda A^2 under P_G, and energy Lambda A^2 times the one-energy caller profile. For F(u)=g(u-wg(u)), the extra composition contributes only the already charged w A^2 first/curl terms. These are actual chain-rule bounds on the pinned source, never derivatives inferred from a law error.

The common carrier pinning uses no inverse c amplification: beta_t/c=q c/sqrt(d_t) and sqrt(v_t)/c=sigma/sqrt(d_t) are at most one. The construction should be declared explicitly in any final native-completion statement.
