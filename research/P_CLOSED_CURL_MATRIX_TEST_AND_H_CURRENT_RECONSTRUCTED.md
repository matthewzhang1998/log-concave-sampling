# Closed-curl matrix-test frame and native outer-H current: reconstructed proof

NEW RECONSTRUCTION from the available mathematical context. Not asserted byte-identical. Original P matrix-test source ae934c3a0c9ddbf4d4f61b5ff9884e2b2db766c0a2dd4d31b80fff6a6fd699a3; original P current-frame source b94cefaf9dfc1ee09e8c542504933952ca6bf0f8e0955915bdbaf872043c2270. Original C swap source ba28738df7cc0724b2787a3cd59f4d75027ae10cd38db25a476cd64a60293911. C/R audits were00ceb87e5e343c6d684983f53cd2ebe000b9eb8580a0d39e6fde57a9e48b790e and b98d63cd67340b2760deb55475b8e898fac922fd2de8f58fcaa5e90d1260e08a. R's literal consumer pin678450180a0be1cf883faa97e1ab16efb93230d45836e55f51dec52c5c16c6cf is a separate source.

## 1. Gaussian swaps and closedness

Let Z be standard Gaussian, N=−L its nonnegative number operator, R=∇N^−1 on centered scalar functions, and δ_j=Z_j−∂_j. For a matrix field C, on positive chaos define

    (T_in C)_ab=Σ_j δ_j R_b C_aj,
    (T_out C)_ab=Σ_j δ_j R_a C_jb.

The operators kill constants. On each fixed chaos these are averaged coefficient-index swaps, hence contractions in the appropriate vector-row L2 space; T_in is self-adjoint. If C is a closed skew2-form,

    ∂_j C_ab=∂_a C_jb+∂_b C_aj,

then multiplying by δ_j and dividing by the chaos degree gives

    T_out C+T_in C=C−EC.

Thus A_out C:=(I−E−T_out)C=T_in C. Skewness also gives T_out C=−(T_in C)*.

For a fixed output row u, the rowwise contraction yields

    ||u* T_in C||L2≤κ|u| whenever ||C||op≤κ.

## 2. The matrix-test frame

For any fixed HS matrix U and centered scalar φ,

    T_in(Uφ)=(UZ)⊗Rφ−U D(Rφ).

This follows directly from δ_b R_jφ=Z_b R_jφ−∂_bR_jφ and the index-swap definition. Self-adjointness gives

    E〈U,T_in C〉φ=E〈C,T_in(Uφ)〉.

The first term is bounded by κ E|UZ||Rφ|. Since E|UZ|²=||U||HS² and ||Rφ||2≤||φ||2, it costs at most κ||U||HS||φ||2. For the second use

    ||U D(Rφ)||_*≤||U||HS ||D(Rφ)||HS,
    ||D(Rφ)||L2,HS≤||φ||2.

It has the same bound. Therefore

    ||〈U,A_out C〉||L2≤2κ||U||HS.                       (1)

The transpose identity gives the same bound for T_out C and consequently

    ||〈U,C−EC〉||2≤4κ||U||HS,
    Var〈U,C〉≤16κ²||U||HS².                            (2)

Known coisometric pullbacks and fixed right contractions preserve the matrix-test frame with their actual operator norms. This is substantially stronger than an entrywise or raw matrix-HS bound and does not assume a derivative bound for C.

## 3. Exact native outer square source

At an original covariance clock c²+v²=1, use the actual square source f=P*E, first A, curl κ, energy e. Let

    f_X(z)=[f(cX+vz)−f(cX)]/(Av),
    I=I_(f_X)(p,z1),     p=P*G+QZ_perp.

The SAME I is reused in every sign query. With fixed β=1/2, γ=sqrt(1−β²), L_δ=[δβI,γI], define on the complete outer square V

    H(V;X,I)=A/(4β²sv) Σ_(ε,δ) ε L_δ*
                  f(cX+vL_δV+vδεβsI).                  (3)

At s=A the actual source has first O(A), curl O(κ), and pointwise VALUE bound ΛA²|I| at EVERY outer V argument. Its p/I first is A², X first A/v and physical caller sqrt(A)/v. The integrated energy is A e/v. Both complete outer blocks, all original anchors, and the original I record remain.

Its exact covariance-versus-forward orientation is

    O_H=Cov_V H−Sym ∫E J_H,t² dt
       =−Sym ∫E[J_H,t Curl_H,t]dt.

The physical projected stencil is

    B O_H B* =−A²/(16β²s²) ∫Σ δεδ'ε' r_(δ,δ')
                     Sym E[P J_(δ,ε) K_(δ',ε') P*]dt,
    r_(δ,δ')=γ²+δδ'β²,

where K is the heat-averaged original Curl f at the actual signed query. Within each factor both ε signs stay grouped. The shared V,X,I remain; only the expressly independent heat shields are independent. Refreshing I separately in the two factors changes the coefficient.

The covariance discrepancy contributed by the old whole-H mean has a PLUS coefficient +q2 α²Σw² B O_H B*, q2=1/3. O_H is not assumed positive; positivity belongs to the complete gapped Gaussian kernel.

## 4. Homotopy and marked defect

For a smooth vector f, a fixed h and q,

    f(q+h)−f(q−h)
      =∇_q ∫_(−1)^1 h·f(q+th)dt
         +∫_(−1)^1 Curl_f(q+th)h dt,

with the corresponding consistent curl convention. Thus the actual H can be written analytically as a genuine gradient plus a remainder r_H satisfying |r_H|≤ΛAκ|I|. This is a proof decomposition, not an executable new gradient oracle. A_V=I−E_V−T_V kills that gradient, and its L2 norm is at most2. Hence W=A_VH has energy at most Λκe/v.

Differentiating H in p, use original Hessian symmetry for the gradient portion and apply the closed-curl matrix-test result to the remaining curl factors. At fixed p and original I labels, all right D_pI factors are fixed in V and have their actual bounded first. This yields the UNIFORM fixed-p frame

    ||〈U,D_pW〉||L2(V)≤Λ Aκ||U||HS.                    (4)

No arbitrary nongradient matrix-field version is claimed.

H and W are odd in FULL p, not merely the physical G with its perpendicular complement fixed. The p-first of H is A². Consequently the fully p-averaged orientation has HS at most ΛA²κe/v and ordinary op at most ΛA³κ. This is an averaged coefficient result, not a license to discard its retained-public field.

## 5. Retained-public current without a dimension trace

Let C=R_p W. The uniform frame(4) passes through the p-Mehler integral, so its matrix-test scale at fixed p remains τ=ΛAκ. Define

    B_ij=E_V Σ_a (D_aH_i) C_ja,
    Γ_ija=E_V H_i C_ja.

Since C has mean zero in V, H may be centered in V in Γ. The key HS factorization is an operator composition: the map from L2(V) to H coordinates has HS norm≤||H||L2(V), while the map U→〈U,C〉 into L2(V) has operator norm≤τ. Therefore

    ||Γ||HS≤τ||H||L2(V),
    ||Γ||L2(p),HS≤ΛA²κ e/v.                            (5)

This is not a raw bound on |H||C|, which could lose dimension. The rank2 B term is bounded by the p-first A² times ||W||2, hence the same A²κe/v. The corresponding ordinary matrix cuts follow from the row covariance of centered H and the matrix-test frame.

Now suppose the already calibrated endpoint, conditional on p, is Gaussian with mean m(p) and covariance Σ(p), and the auxiliary V defining the coefficient is independent of that endpoint. Apply one Gaussian p-Riesz transfer to the W factor. The exact terms are

    rank2: B contracted with D²φ,
    rank3: Γ D_pm contracted with D³φ,
    rank4: (1/2)Γ D_pΣ contracted with D⁴φ.

Including the initial covariance interpolation factor1/2 gives coefficients c/2,c/2,c/4. If D_pm is bounded and D_pΣ=O(A³), their one-energy allowances are A²κe/v, A²κe/v and A5κe/v. At κ=A² these are A4e/v and A7e/v.

All original G,X,I observers are retained through the complete endpoint. Only auxiliary V is an independent coefficient representation. A marginal old-G law cannot be appended after G has been integrated. This exact chronology is the separate R consumer's responsibility.
