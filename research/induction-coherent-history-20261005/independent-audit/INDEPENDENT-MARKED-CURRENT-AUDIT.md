# Independent marked-current rank-two audit

2026-10-05. The original construction and every prior file were preserved. This audit concerns `REVIEWED-CURRENT-SNAPSHOT.md`, SHA-256 `440dff2bb37da96f65b58765422cf4e107b6b1538780be8c7589ce64a02ec114`. Subsequent author revisions need their own comparison against this snapshot.

## Verdict

**PASS for the stated local rank-two weighted-current construction and the genuine retained-mark coupling, after the assumption/count clarifications below. No native admission follows from the likelihood identity.** The exact Gaussian formulas and the stated Hermite truncation constants are correct. The mark's correlations with the old source bank are preserved inside the weighted identity. The unweighted new mark does not retain the old mark law. The construction does not supply the missing base/singleton history source, an unrestricted-dimensional history comparison, or the next mean induction.

Four clarifications should be incorporated:

1. The fresh Gaussian must satisfy `Law(G | Y,Omega)=N(0,I)` and be conditionally independent of the **complete** source/mark bank. The statement `M independent of G given (Y,Omega)` alone is vacuous when M is Omega-measurable. If M uses additional randomness, that randomness belongs in Omega.
2. Ordinary `W2` on the product Euclidean space requires finite conditional second moments of M and of any retained labels. With heavier-tailed labels the same zero-label-displacement coupling cost can still be written, but it should not be called a comparison in the usual `P2` Wasserstein space without qualification.
3. The prior fourteen-VALUE source charges at most `17d` Gaussian coordinate roots **plus scalar clock randomness**, and separately charges any inactive-slot buffer. The new mark adds its genuinely new roots. These qualifications are absent from the snapshot's abbreviated total `17d+d_M`.
4. Rank-two, dimension-independent likelihood energy is a **single rank-two update** result. For independent physical blocks, a global likelihood is the product of the block likelihoods; its rank and energy amplification accumulate. The genuine retained-mark blockwise W2 estimate does not have this defect, but the one-block likelihood constant cannot silently be reused for that global product.

## 1. Exact Gaussian computations

Fix the complete source bank. With `B=H/(2v)`, `K=I+B`, `E=K^2-I`, the conditional law of `Z=sqrt(v)KG` is `N(0,vK^2)`. The ratio against `X=sqrt(v)G` is exactly

    L(X)=det(K)^(-1) exp[X*(I-K^(-2))X/(2v)].

The positive determinant has at most two nontrivial factors. The sign of an eigenvalue of H is immaterial to positivity of K under `||B||<=h<=1/8`. The covariance is exactly

    E[ZZ* | Y] = vI + C(Y) + E[H^2 | Y]/(4v).

Conditional change of variables proves the identity with any measurable test of `(Omega,X)`, subject to actual weighted integrability. That is the correct way to preserve the correlations between M and H. Merely saying a test is integrable against one unweighted carrier would not guarantee the two marked expectations exist.

For `p>=1`, the exact Gaussian integral is finite if and only if `I-(p-1)E` is positive definite, and then

    E_G L^p = det(I+E)^(-(p-1)/2) det(I-(p-1)E)^(-1/2).

The author's sufficient condition `(p-1)e<1`, with `e=2h+h^2`, is valid. At the stated guard, `e<=17/64`, and p=3 is safely allowed. At p=2 this becomes `det(I-E^2)^(-1/2)` as claimed. Since rank(E)<=2,

    E_G L^2 <= 1/(1-e^2),
    E_G (L-1)^2 <= e^2/(1-e^2).

Multiplying these conditional inequalities by `|M|^2` before integrating the source bank gives (5.2)-(5.3). No independence between M and H, no fourth moment of M, and no bound on its first are needed.

A uniform higher-moment constant usable below is

    C_p = [(1-e)^(-(p-1)/2) (1-(p-1)e)^(-1/2)]^2,

for `(p-1)e<1`; it is conservative but rank-two and dimension independent.

## 2. Genuine joint Wasserstein comparison

Under full-bank freshness, `Zstar=(vI+C)^(1/2)G` is independent of Omega conditionally on Y. Use the same `(M,Omega,G)` on both sides. Inserting `sqrt(v)(I+C/(2v))G` gives the exact first difference norm

    || (H-C)G/(2sqrt(v)) ||_2 = ||H-C||_(L2;HS)/(2sqrt(v)).

For every c>=0, the scalar concavity remainder lies between zero and `c^2/(8v^(3/2))`. Spectral calculus and Minkowski give (3.1), with exactly its factor 1/8. The analytical square root is not a producer leaf. The mark displacement is exactly zero.

For physically independent blocks, `E||H_j-C_j||HS^2<=R_j^4`. With each integer `d_j>=1`,

    (sqrt(d_j)+sqrt(2L))^4 <= C(d_j+L)^2,
    sum_j (d_j+L)^2 <= C(b+L^2)D.

This establishes the loose blockwise fluctuation term in (3.2). If `C_j<=ell^2 I`, its square-remainder term is at most `ell^4 sqrt(D)/(8v^(3/2))`. The upper bound `C_j<=ell^2 I` requires the stated Gaussian conditional-centering/Poincare structure, not merely a cap of size R. For the imported source it follows by retaining the appropriate non-inner labels and using its conditionally centered inner-Gaussian first bound.

Independence of blocks is a sufficient way to realize the asserted product Gaussian reference. The first HS sum itself does not require independent H_j when the constructed map is explicitly block diagonal.

## 3. Hermite formulas and a sharper optional tail bound

In an eigenbasis, the likelihood is a product of one-dimensional Gaussian variance-ratio expansions. For eigenvalues `lambda_1,lambda_2`, the degree-2k component has squared norm

    a_k = sum_(a+b=k) c_a c_b lambda_1^(2a) lambda_2^(2b),
    c_a = binom(2a,a)/4^a.

Consequently `sum_k a_k t^k=det(I-t E^2)^(-1/2)`. This proves (6.1)-(6.2), including all factorials and Wick normalizations. The generating identity is understood in its convergence range, in particular at the author's selected `t=1/(2e^2)>1`.

At that t each determinant factor is at least 1/2, so the generating sum is at most 2. The author's tail bound `sqrt(2)(sqrt(2)e)^(p+1)` is correct. An optional stronger uniform rank-two estimate is

    ||L-L_[p]||_2 <= e^(p+1)/sqrt(1-e^2).

Indeed `sum_(a+b=k)c_a c_b=1`, so `a_k<=e^(2k)`. This improvement is not needed for validity.

Conditioning on Omega before multiplying by M proves the relative-energy estimate without taking a derivative of M. Bounded tests, and bounded-gradient current tests, then follow by Cauchy-Schwarz. The truncated polynomial is generally signed; calling its underlying pushforward a positive probability law does not make this polynomial a positive density. The snapshot distinguishes these correctly.

## 4. Readsets, positivity, dimensions, and costs

The executing exact packet owns the old complete source bank and one new Gaussian carrier. Its source-to-mark map is explicit, and no completed LAW is converted to RAW noise. Retained labels must be unread by the new Gaussian draw, and changing a source argument must replay every dependent original call. Correlation of M with H is allowed and essential.

Rank-two linear algebra uses O(d) projections/inner products and fixed-size matrix operations. Once the original H and M graphs are supplied, likelihood evaluation adds no original gradient VALUES; output multiplication costs O(m). The degree-p expansion adds O(p^2) scalar work. The first-five-call mark from the imported shifted-current note therefore adds at most five VALUES, for an upper bound of nineteen, with fewer only by exact-key reuse.

This cost statement is not permission to execute a true integrated history at unit cost. Nor is it a certified bit-complexity analysis against an unspecified original-g oracle. The snapshot is appropriately explicit that finite-precision source/scalar/clock floors must be charged separately.

For block-diagonal H with B physical blocks, the total likelihood second moment can be as large as `(1-e^2)^(-B)` under independent blocks. A local mark cannot simply be multiplied by its own block likelihood for arbitrary global tests. Such a simplification needs a block-separable test/readout theorem or a separate conditional construction.

## 5. Why exact likelihood does not imply a bounded native first

Already for deterministic `F=F'=R e_1` and M=1,

    L(G)=(1+h)^(-1) exp[a G_1^2/2],
    a=1-(1+h)^(-2)>0.

Both L and its G-first are unbounded. Small h and dimension-independent L2 energy do not change this fact. A finite Hermite truncation is also generally a polynomial with unbounded first.

Moreover, the rank-two projected-root event in section 8 bounds values, not all source/caller firsts. Let `K(theta)=I+b u(theta)u(theta)*`, with `u(theta)=(cos(theta),sin(theta))`. At theta=0,

    D_theta log L = [1-(1+b)^(-2)] G_1 G_2.

The currently active projection `|G_1|` may be bounded while this derivative is arbitrarily large through G_2. Differentiating only Ksmall misses motion of its basis. Thus the section-8 small-matrix derivative observation is correct only for that restricted derivative; it is not a complete derivative certificate. An invariant full-matrix derivative avoids this issue:

    D_K log L[V] = -tr(K^(-1)V)
                  +(K^(-1)G)* V (K^(-2)G).

No derivative of eigenvectors is needed in this formula.

## 6. One valid finite clipping certificate

This section supplies a conservative local certificate, not a native admission result. Assume fixed v and frozen scalar modes, global source firsts `Lip(F)<=L_F`, `Lip(F')<=L_F'`, and cap `|F|,|F'|<=R`. Let `tau_B` be the radial projection of the **complete carrier** G onto the radius-B ball, and clip the likelihood value at T:

    w_(T,B)(Omega,G) = min{L_H(tau_B(G)), T}.

This finite graph is globally Lipschitz and has an almost-everywhere first. Smooth radial/value caps, identity below their lower thresholds and constant above their upper thresholds, provide a C1 version with fixed-factor changes in all displayed constants.

Write

    alpha=(1-h)^(-1),
    a_h=(1-h)^(-2)-1,
    J_K=R(L_F+L_F')/(2v).

For every unit source direction, both the nuclear and operator norms of DK are at most J_K. Consequently a valid complete first majorant is

    Lip(w_(T,B)) <= T[J_K(alpha+alpha^3 B^2)+a_h B] =: L_w.

For a caller direction use its actual source firsts instead; varying v, caps, or clocks adds their actual paths. This bound includes rotation of the source span. It is independent of a Gram determinant or basis condition number.

Let `q_B=P(|G|>B)`. For p>2 with `(p-1)e<1`, condition on Omega to obtain

    ||M(L-w_(T,B))||_2
      <= ||M||_2 [sqrt(C_p) T^(-(p-2)/2) + T sqrt(q_B)].

The first term is the likelihood-value tail, using `L^2 1_(L>T)<=T^(-(p-2)) L^p`; the second compares two weights in [0,T] that differ only outside the carrier ball. Only an L2 mark moment is needed. A full-d carrier cap has a chi-d tail and can introduce dimension-dependent derivative constants. It is not the original rank-two value cap.

To bound the product's first, also require a bounded mark `|M|<=S` and `Lip(M)<=L_M` (or explicitly cap a mark with such a first). Then

    Lip(M w_(T,B)) <= T L_M + S L_w.

For any supplied square lift m of the mark on the complete input, with `|m|<=S`, first L_M and skew-Jacobian bound kappa, the product has

    ||D(wm)-D(wm)*|| <= T kappa + 2 S L_w.

This follows from the product rule. Mark clipping preserves the first upper bound but need not preserve a sharper preexisting curl bound; the generic bound `kappa<=2L_M` is always available. L2 mark energy alone supplies neither a global first nor a global supremum. These constants may be too large for an intended native theorem, whose own source normalization, numerical guard, complete square lift, caller ports, and positive allocations still must be checked.

## 7. Audit of the proposed finite first-correction extension

In correspondence after the snapshot was pinned, the author proposed using only the first current correction

    N = M T1(E,G),  T1=(G* E G-tr E)/2.

For norm-one clips of M at R_M and the complete G at B, and a uniform E-first L_E, the claimed bound is valid:

    L_N <= L_M e(B^2+2)/2
           +R_M L_E(B^2+4)/2 + R_M e B.

The derivative DE in any one direction has image contained in `span(F,F',DF,DF')`, hence rank at most four. Also

    L_E <= (1+h)L_H/v,   L_H<=R(L_F+L_F').

The trace bound in the middle term is therefore legitimate even as the rank-two span rotates. The un-clipped conditional Gaussian norm satisfies `||T1||_(2|Omega)<=e`, since `E_G T1^2=tr(E^2)/2`. Clipping can be paid by Gaussian tails; alternatively its crude energy remains at most a fixed multiple of e times the mark energy. For retained clipped G, use the explicit conditional profile bounded by `e(B^2+2)/2` instead.

Retuning the imported secant from `delta=e_smooth A^2 sqrt(tau)` to `delta=e_smooth sqrt(tau)` changes its scaled mean-field error from `O(A^5 n)` to `O(A^3 n)`, hence its covariance error to `O(A^5 n^(3/2)+A^6 n^2)`. At the old dyadic cutoff `d_K~A^(10/9)`, the complete source first improves from `A^(-26/9)` to `A^(-8/9)`. For fixed physical block size and charged logarithms, this gives

    L_N=O(A^(37/9)/v),   energy(N)=O(A^7/v).

These grades concern the retuned finite local current source, not a construction of the missing true-history mark.

If this actual bounded square source is normalized by `1/sqrt(u)`, set `l=L_N/sqrt(u)<=1`. Padding its declared radius to `r=sqrt(l)` is legitimate: first <=r and generic relative curl <=2r. With `mu=r` and relative energy `delta_native=energy(N)/(sqrt(u) r sqrt(n))`, the five non-baseline terms of the imported raw-split bound reduce to

    C energy(N)[r^2+r^2+2r^2+r^(5/2)+r^3]
      <= C energy(N) l
      = O(A^(100/9)/(v^2 sqrt(u))).

At `u=v` of order `A^2` times charged logs, r has the positive power `A^(5/9)`. This validates the arithmetic and the proposed padding majorants. It does not waive the imported numerical native/capture/source/caller guards or the complete occurrence/root count. A genuine-gradient baseline identically zero can be emitted exactly; that simplification must be part of the actual compiler rather than silently deleting a nonzero theorem term.

**Target warning:** `E_G[N | Omega,Y]=0` exactly before clipping. Applying an own-mean compiler that integrates G yields a near-zero Gaussian mean and does not preserve the current `E[N grad(phi)(sqrt(v)G)]`. For a useful conditional current contract, G must remain a retained caller, Omega must be the private integrated bank, and the later test coefficient must read only retained callers. A Gaussian auxiliary mark noise contributes zero to such a linear test only under its proved conditional centering/independence. That does not preserve the old mark law, create a small-energy raw output, or authorize reading the consumed private bank. The author's planned revised contract was not yet available when this audit section was written.

## 8. Executed diagnostics

`check_current_resummation.py` uses deterministic 96-node Gauss-Hermite quadrature in each of two dimensions. It checks p=1,2,3 likelihood moments, degree-0 through degree-14 Hermite norms and truncation tails, an exactly finite correlated mark/source bank, the retained-mark coupling bound, the rotating-range derivative example, and rank(DE)<=4 on deterministic seeded test matrices.

All checks pass. In the finite-bank example the same-G endpoint coupling cost is `0.1300132926885364`, below the stated bound `0.1306322047603437`; its correlated-mark polynomial identity residual is approximately `2.22e-14`. Details are in `current_resummation_checks.json`. The analytic arguments above, rather than these diagnostics, establish the general claims.
