# Independent audit: partial-u filter and positive true-Gram covariance return

Date: 2026-10-04.

## Verdict, dependencies and scope

**PASS as a bounded, source-qualified finite original-VALUE true-Gram covariance service**, under the actual stated imported finite-mean/radius/filter/numerical guards and fixed positive buffer. Combining the source with the separately audited mixture lemma and resolvent/clock comparison supplies a positive return targeting

    N(0,v0 I+Cov(H_j|Z))

at integrated fresh-standard-Z error `Lambda A^4 sqrt(D)` plus separately restored floors.

This closes the covariance service for the deterministic conditional-mean field j along the true Markov history. It does **not** close the mixed K-current, same-carrier mean m3, full F2 covariance, or full order-four endpoint. In particular, it does not claim that H_j alone is F2.

Pins directly inspected:

- `../../order-reentry/RECTANGULAR-FIRST-COEFFICIENT-MEAN-AND-GRAM-RETURN.md`: `13b0fad9f69355c511c8884ae1e7691066f740dccf387ec453344479193154ea`.
- Coefficient source: `359cdfffa1ce4b7467d24ed24d6b8efdedb85a16c5de0b650d3b4ea37aa52cfd`.
- Corrected mixture proof: `7476219bacabfacc44030a14710408735d3e31b78103376d57f525a3518726dc`.
- Original imported LOW30 `/workspace/scratch/e90bd698ae04/research-source/High Acc Ideas/ai-bucket/30_low_acc.tex`: `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`.

The LOW30 portions read directly are `t30:eq:seed-data`, `t30:eq:first-response`, `t30:eq:first-contract`, `t30:lem:value-mean` and its actual program/first/caller bounds, plus `b27:compiler:mean` and its original-label/small-radius conditions. The high-order compiler is imported as the already established finite theorem, not re-proved or numerically instantiated here.

## 1. Literal split and stronger near-gradient VALUE energy

The source's split `V_X=V_g,X+V_e,X` is a pointwise identity on the same `(u,H)` record. V_g is the anchored original gradient packet. Its u Jacobian is a positive sum of original symmetric Hessians and is symmetric. No gradient claim is needed jointly in X and u.

For `E_Q=J_Q-g`, the elementary bound is

    |E_Q(x,H)|<=A|H_Q(x,H)|
               <=A^2(|x|/2+beta_in |H|).

Consequently the anchored V_e difference has the explicit bound

    |V_e,X(u,H)|
       <=A^2[S2|X|+|u|/4+2 beta_in S1|H|],

where `S1=sum w_i r_i/c_i` and `S2=sum w_i r_i^2/c_i` are uniformly bounded for the admitted dyadic rule. This proves its advertised small VALUE profile. Its actual u and X firsts remain O(A), and its lifted curl/H first remain O(A^2), by the same literal chain rule as in the coefficient-source audit. No derivative of this small energy has been used.

The anchor `J_Q(r_iX,H)` is private-H dependent and remains inside each complete raw occurrence. V_e(0,H)=0 by the same exact two-branch cancellation.

## 2. Partial-u Hermite filter and unchanged H

Fix H temporarily and expand `u -> V_X(u,H)` in orthonormal vector-valued Hermite chaos. For one b, averaging over the standard z in

    [V_X(sqrt(1-b^2)z+bp,H)-V_X(sqrt(1-b^2)z-bp,H)]/(2b)

kills even degrees. Its degree `(2h+1)` multiplier is `b^(2h)`. The degree-one coefficient is `[E_u D_u V_X(u,H)]p`.

The finite coefficient rule has `sum d_l=1`, annihilates h=1,...,K-1, and has `sum |d_l|<=4`, `b_l<=1/2`. Hence every remaining odd multiplier is bounded by `4*4^(-K)`. Orthogonality in the actual Gaussian p gives the exact dimension-free conditional error bound by this constant times the centered L2_u energy. That energy is uniformly at most `(L_j/2)sqrt(D)`, from the u first. Jensen in the unchanged H gives the claimed target `B_Q(X)p` and the uniform-in-X error `C A sqrt(D)4^(-K)`.

H is held unchanged for every sign and filter node. Scaling H with the u Mehler contraction would evaluate a different finite conditional response. The independent checker explicitly detects this difference on a nonlinear two-variable fixture. The source correctly uses the partial-variable rule rather than assuming an unproved restriction estimate on a full-space filter evaluated at the nongeneric public point `(p,0)`.

For clarity, the original printed Chebyshev-Lobatto rule allows its endpoint `b=1/2`; the response is perfectly well-defined there. If the new source's strict `b<1/2` convention is desired literally, choose the squared-node upper endpoint 0.24 and lower endpoint `1/(4K^2)` for K>=2. The same interpolation identities hold. The absolute extrapolation weight sum is the Chebyshev value `cosh(2(K-1)atanh(1/(2K sqrt(0.24))))`, bounded by `cosh(1/sqrt(0.24))<4` for K>=2; all the displayed inverse-width/tail bounds remain valid. At the required small-A grades K can always be chosen at least two. The checker uses this strict-interior version, so this harmless endpoint convention does not block the construction.

K grows only logarithmically for `4^(-K)<=A^3` after the chosen public-log error allocation. Each response is p-odd pointwise and is exactly zero at p=0 with identical-query reuse. This is conditional-mean calibration, not a strong approximation to B_Q or B_Q p.

## 3. Actual gradient and lifted near-gradient consumer ports

For fixed captured X,p, the full z derivative of F_g is a signed sum of differences of symmetric D V_g matrices. It is symmetric even though d_l need not be positive; convexity is not required by the imported genuine-gradient mean interface. Its actual radius is O(KA).

For F_e, each z derivative costs `|d_l|c_l/(2b_l)` and each H derivative costs `|d_l|/(2b_l)`. Therefore

    Lip_(z,H) F_e <= C K A,
    ||Curl(P_z*F_e)||op <= C K A^2,
    ||F_e||_Lp <= Lambda A^2(|X|+|p|+sqrt(D)).

The off-diagonal curl block is its actual H derivative; it is not dropped because the output has dimension D. The square lift is on the actual 2D-dimensional `(z,H)` source, with `P_z=(I,0)`. The imported square near-gradient mean may be run there and projected by P_z; the resulting output covariance is I_D. Its active dimension and complete tape must be charged accordingly.

At z=H=0 and retained nonzero p, the source's private origin need not vanish. The complete response executed at those zero private inputs is a genuine caller-only X,p origin. Its bound follows pointwise from the V_e estimate, giving `C Lambda A^2(|X|+|p|)`. Subtracting and restoring this recorded value gives the required anchored energy profile without assuming that a Gaussian Lp bound controls an arbitrary zero value. The checker exhibits nonzero retained-p origins.

The F_g origin is likewise captured from its actual program. All its raw H-free gradient anchors and all of F_e's remaining H-dependent anchors keep their appropriate source/version keys. The source-zero and p=0 cancellations are literal, not inferred from a target law.

The actual X and p firsts are O(Lambda A). In particular, the small F_e VALUE estimate does not improve the X first to O(A^2). Origins' caller derivatives are restored through their own recorded graphs.

## 4. The imported mean error really reaches order four at these ports

LOW30's near-gradient mean theorem has the literal error

    Lambda e {ell(a+mu)+ell^3(1+mu^(-1/2))}+floors,

where its curl majorant is `ell*a`. In the new source take

    ell=C K A, kappa=ell*a=C K A^2,
    e=Lambda A^2(|X|+|p|+sqrt(D)), mu=A,

including the fixed variance normalization in the constants. The bracket is

    kappa+ell*mu+ell^3(1+A^(-1/2))
      <=Lambda A^2

for `0<A<=1`, with the displayed K factors retained in Lambda. Thus the near-gradient branch error is `Lambda A^4(|X|+|p|+sqrt(D))`.

The genuine-gradient branch at any fixed order at least four has error `Lambda sqrt(D)(C K A)^4`, which fits the same allowance. Its mean/constant origin is restored exactly. Independent complete branches normalized by fixed shares 1/2+1/2 yield target covariance I_D and target mean E(F_g+F_e)=E F.

This admission requires every imported threshold at the actual normalized C K A radius, actual C K A^2 curl, lifted active dimension, padding, fixed shares, finite clocks and precision. A<=1/2 alone is not sufficient. The source explicitly imposes these guards.

The imported actual residual first is

    Lambda(ell+ell^2 mu^(-1/2))=O(Lambda A),

and its caller multiplier is `Lambda(1+ell mu^(-1/2))`. At the actual O(Lambda A) X,p source first and fixed guards this gives the stated O(Lambda A) caller port. These are finite-program first bounds; no law error was differentiated.

Combining with the partial-u calibration and integrating fresh p proves the source's conditional mean-law statement. Only captured X,p remain retained at that cut. All z,H, filter and mean-bank roots are integrated, and none is reattached as an observer.

## 5. True Gram orientation and positive variance allocation

After the preceding comparison, at fixed X the analytical target is

    B_Q(X)p+N,

with p and the unit Gaussian N independent. Its covariance is `I+B_Q B_Q*`, exactly. N is a comparison variable; the proof does not identify it with an old actual source carrier or differentiate that law replacement.

For the final positive join, let `v_mean=v0/2`, `v_free=v0/2`, `nu_j=v_mean a_j`. The mean source is scaled once by `1/sqrt(v_mean)`, independent of j. Each returned contribution has reference

    sqrt(nu_j)[B_Q(X_j)p_j/sqrt(v_mean)+N_j]
       =sqrt(a_j)B_Q(X_j)p_j+sqrt(v_mean a_j)N_j.

All complete banks and p_j are independent conditional on their captured X_j, and the final `sqrt(v_free)Z_free` is untouched and independent. Therefore, conditional on all owned X_j, the full covariance is exactly

    v0 I+sum_j a_j B_Q(X_j)B_Q(X_j)*.

There is no inverse-a_j radius. The weights sum to one and all variances are positive. The actual source-zero carrier is the sum of the scaled actual mean carriers and Z_free, with covariance exactly v0 I. The independent checker verifies these block-Gaussian row identities for nonsymmetric B_j.

The conditional W2 errors add with their literal `sqrt(nu_j)` multipliers. The bound `sum sqrt(a_j)<=sqrt(N_q)` is a public-log factor, not a claimed independent-error cancellation. Averaging the actual owned `X_j=q_jZ+sqrt(1-q_j^2)G_j` moments gives the stated `Lambda A^4(|Z|+sqrt(D))` profile.

The separately audited mixture theorem then replaces this Gaussian covariance mixture by `N(0,v0 I+C_Q(Z))` with a uniform `C_v0 A^4 sqrt(D)` error. Its actual Gaussian-to-HS Lipschitz hypothesis holds by the row-Bessel derivative calculation. This is the substantive law step; matching the covariance alone would not suffice.

Finally the audited resolvent stability and clocks restore C_Q to `Cov(H_j|Z)` in L2 of the actual fresh standard Z. The service does not promote this last bound to uniform arbitrary-caller control.

## 6. Complete VALUE recurrence, roots and precision

The safe raw costs are

    Q_g=2N_r, Q_V=2N_r(N_in+1),
    Q_Fg<=2K Q_g,
    Q_Fe<=2K(Q_V+Q_g).

Every sign/node calls the complete raw source, including all nested H_Q ancestors and terminal original g VALUES. Exact within-record cancellations may reduce costs; the safe recurrence assumes none beyond explicitly identified keys. The checker independently counts all 2K signed V occurrences.

For the actual imported complete occurrence counts N_B,N_E, one conditional mean call has the advertised recurrence

    Q_captured(X,p,origins)+N_B Q_Fg+N_E Q_Fe
                               +known/numerical/replay work.

Every q bank owns a fresh X_j and its own p_j, origins and complete gradient/near-gradient banks. The source's original counts do not treat a completed mean as one original gradient query. The final tape includes all replay source dimensions, mean/filter/covariance clocks, p_j, G_j and variance fills. At fixed grade, products of the admitted public-log counts remain public-log; no uniform-in-growing-order conclusion follows.

The filter's absolute VALUE error multiplier is at most `sum_l |d_l|/b_l<=8K`. This multiplies the already audited raw V absolute error ledger, followed by the finite mean's actual readouts and `sqrt(nu_j)` weights. There is no division by realized source energy. Numerical derivative/first accuracy, original physical mode/caller restoration and replay are their separate imported obligations, exactly as in the coefficient-source audit. HVPs occur only in requested first/adjoint sweeps at recorded original VALUE sites.

## 7. Evidence and remaining gates

The companion checker passes 277 new assertions and reruns the 1,930 coefficient-source assertions. Its filter/Gram checks cover admissible signed coefficients, killed Hermite moments and tail bounds, the unchanged-H distinction, literal signed source counts, p oddness/zero, nonzero captured origins, partial-variable curl bounds, true-Gram orientation, positive total variance and actual source-zero carrier variance. It does not numerically execute the imported high-order mean programs.

The new finite covariance service and the mixture lemma together close the first two construction obligations listed as open in the earlier audit of the frozen coefficient-only source. They do not retroactively turn that earlier source alone into a covariance algorithm.

K, the same-carrier m3 mean, original restoration floors and the complete order-four endpoint join remain separate. No all-order typed closure or full endpoint PASS is claimed.
