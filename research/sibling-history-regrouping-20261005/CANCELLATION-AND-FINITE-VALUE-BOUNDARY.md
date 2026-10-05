# Shared-history cancellation: an exact rank-eight packet, a finite scalar VALUE producer, and the remaining boundary

2026-10-05. Local mathematical work. No external calls, uploads, Picard construction, or cold restart. This is a new check of the pinned rank-eight obstruction, not an arbitrary-order native theorem.

## 1. Useful conclusions

1. There is an exact, finitely generated 64-history packet containing the obstructed double star. Its shared-clock/common-bank sum cancels the problematic **difference-frequency component** before absolute values.
2. This cancellation does **not** remove the old random coefficient's obstruction. Its remaining common-frequency components still have variance asymptotic to a positive constant times K^4. The exact constant is computed below.
3. Integrating by parts in the deepest shared structural innovation gives an original-g VALUE conditional-mean producer. Primitive-clock importance weights remove all primitive inverse shields from its actual VALUE expression and its retained-caller first. In the scalar interior-parent sector it uses 16 g evaluations per complete sample, has explicit radius/energy/first bounds, and has a direct positive polynomial consumer with an explicit finite Monte Carlo bill.
4. This is not an equality of the old random coefficient's law. It can preserve the joint law of retained callers/observers only through a stated coupling error and only when the integrated roots are internal and unobserved. The primitive shield improvement does not resolve the deepest structural endpoint, multidimensional one-energy requirement, or admission to the imported native pair/filter architecture. No eventual-sublinear complexity statement follows.

## 2. Exactly which generated histories are grouped

Start with the literal rank-two sibling AST

    R2[(R2 g'_1)(R2 g'_2)].

For n=3,...,8 use the i=n-1 term of the connected generator, append one fresh rank-one leaf, and allow its derivative on the old subtree to hit either of the TWO ORIGINAL siblings. Do not pick one favored split. This makes 2^6=64 actual labelled generated histories. Each has integer coefficient

    2 * 3 * ... * 8 = 8! = 40320.

If p of the six hits reach sibling 1, its primitive resolvent is R_(p+2) and its spatial force order is p+1. Sibling 2 has R_(8-p) and force order 7-p. All seven structural resolvents are R8. The six appended leaves are R2 and have force order one. There are binomial(6,p) histories of this type, with their distinct actual edge labels intact. The p=3 family contains the (4,4,1,1,1,1,1,1) obstruction.

Use one common clock/base-measure realization. Let rho_2,...,rho_8 be the structural clocks and r_1,...,r_8 the primitive clocks. The common positive base measure is

    dmu = product_(j=2)^8 rho_j^7 d rho_j * product_(i=1)^8 r_i d r_i.

The extra primitive powers are r_1^p r_2^(6-p). They belong to the common-parent derivative, below. All structural derivative transport is ALREADY in the rho_j^7 densities; multiplying an additional transport factor would be a mistake.

The same Gaussian column dictionary is used for every history. In particular the two siblings share the deepest parent

    X = rho_2 Y + s G,      s^2=1-rho_2^2,

where Y is the next outer parent and G is the deepest structural innovation. The six added leaves are outside this subtree, so none depends on G. Define

    sigma_i^2=(1-r_i^2)/2,
    q_i=r_i X+sigma_i W_i,  i=1,2,
    B_i(q_i)=E_Z g'(q_i+sigma_i Z).

Write L for the product of the other six privately averaged g' factors. Conditional on all clocks and all roots except G and the private shields, L is independent of G. The actual 64-history sum is

    8! * L * D_X^6 [B_1(r_1 X+sigma_1 W_1) B_2(r_2 X+sigma_2 W_2)],

integrated against dmu. This is literal Leibniz regrouping, not an expectation/cumulant oracle. Every physical mark is retained. The scalar formula is sufficient for the obstruction test.

A pre-existing finite quadrature can be regrouped this way only if its different primitive R_k orders share a common base grid with weights w_j r_j^(k-1), or if it is rebuilt into such a rule with its numerical debt paid. Different order-dependent clock grids cannot silently be paired.

## 3. What cancels, and what does not

For g_K(z)=a z+(b/K)sin(Kz), set

    beta_i=b exp(-K^2 sigma_i^2/2),
    u_i=K r_i X+K sigma_i W_i.

The exact sum is

    D_X^6(B_1 B_2)
      = -K^6 {a beta_1 r_1^6 cos u_1 + a beta_2 r_2^6 cos u_2
          +(beta_1 beta_2/2)[(r_1-r_2)^6 cos(u_1-u_2)
                            +(r_1+r_2)^6 cos(u_1+u_2)]}.

At equal primitive clocks, the cos(u_1-u_2) term is EXACTLY zero. This is the component whose positive interior OU covariance was isolated in the old single-history obstruction. At unequal half-heats sigma_i^2=c_i/K^2, with fixed c_i in a compact positive interval, r_1-r_2=O(K^-2); its unweighted coefficient is O(K^-6), and the two endpoint masses make it O(K^-10). Thus the cancellation is stable under comparable unequal endpoint clocks, although exact vanishing requires equality.

However the single-frequency and common-double-frequency modes remain. Take sigma_1=sigma_2=1/K, center clock masses K^-2 each, and put all other primitive widths and structural clocks in fixed interior panels. The six actual outer g'_K factors equal a^6 plus a uniform exponentially small error. Suppress only the common positive 8!, fixed clock masses and the base factor r_1 r_2, the last of which tends to one. With T=K r X and Var(T) tending to infinity,

    L_group/K^2
      = -a^6 r^6 [a b e^(-1/2)(cos(T+W_1)+cos(T+W_2))
                  +32 b^2 e^(-1) cos(2T+W_1+W_2)] + o_L2(1).

Consequently

    lim Var(L_group)/K^4
      = a^12 [a^2 b^2 e^-1 (1+e^-1) + 512 b^4 e^-2] > 0.

For a=1/2,b=1/4 this constant is 0.00006800129310109894 before the displayed suppressed fixed factors. The original source remains smooth and has a-b <= g'_K <= a+b uniformly.

Therefore this natural complete Leibniz packet is a real cancellation result, but it is NOT yet cancellation of the full coefficient-covariance obstruction. A claim that summing its actual derivative-hit histories makes the old-bank L2 budget cutoff-uniform is false.

## 4. A legitimate conditional-mean change of representation

For every fixed interior deepest clock,

    E_G D_X^6 F(rho_2 Y+sG) = s^-6 E_G[H_6(G) F(rho_2 Y+sG)],

where H_6(z)=z^6-15z^4+45z^2-15. No derivative of the six outer leaf factors appears because G does not enter their rows.

At every force occurrence use an independent private standard normal Z_i and the original-g VALUE score

    S_i=Z_i [g(q_i+sigma_i Z_i)-g(q_i)]/sigma_i.

Then E_Z S_i=E_Z g'(q_i+sigma_i Z_i), exactly, and |S_i|<=A Z_i^2 under |g'|<=A. Thus the conditional expectation of the full packet has the original-VALUE representation

    8! integral dmu s^-6 E[H_6(G) product_(i=1)^8 S_i].

There are 16 original g queries in one complete sample: one value and one same-center anchor per force occurrence. There are no g'', high derivative, expectation, covariance, or cumulant tensor inputs. All q_i are the original ancestry rows. A single shared sample realizes the grouped expression; its occurrences are not independently re-centered.

At a fixed clock list the normalized sample before its 8! and clock weights has bounds

    L2 <= A^8 s^-6 sqrt(720) 3^4,
    caller derivative L2
       <= 2 A^8 s^-6 sqrt(720) 3^(7/2)
                    sum_i |D_y q_i|/sigma_i.

The private-score proof uses only original first derivatives as analysis; implementation uses only VALUE. In particular the caller bound does not secretly assume a bounded g''.

This is equality of CONDITIONAL EXPECTATIONS over the integrated roots, not equality of random coefficients and not a joint-law identity with G retained. If G, any replaced private root, or another integrated bank is an actual retained observer, this transformation is not licensed without a new conditional target or a separate observer-current analysis. Merely being stored for replay does not make an internal root a law observer.

## 5. Primitive-clock importance sampling cancels inverse widths in the actual source

Fix a deepest structural interior sector 0<rho_2<r_*<1. Let

    I_* = integral_0^r* rho^7/(1-rho^2)^3 d rho.

This is an explicit finite scalar normalizer. For u_*=1-r_*^2,

    I_* = 3/4 + 1/(4u_*^2) - 3/(2u_*) -(3/2)log(u_*) + u_*/2.

Sample its deepest structural clock with density rho^7/[(1-rho^2)^3 I_*]. Sample the other six structural clocks independently with density 8rho^7.

For each primitive common base measure r dr, use the probability density

    p(r)=r/sqrt(1-r^2),  0<r<1.

Equivalently take U uniform on (0,1), r=sqrt(1-U^2), sigma=U/sqrt(2). Its exact likelihood ratio for the original r dr is U. Multiplying the normalized source BEFORE bounding gives

    J_i = U_i S_i/A
        = sqrt(2) Z_i [g(q_i+sigma_i Z_i)-g(q_i)]/A.

The division by sigma has disappeared from the executable expression. Hence

    |J_i| <= U_i Z_i^2 <= Z_i^2,
    |D_q J_i| <= 2 sqrt(2)|Z_i|,
    |D_Z J_i| <= 2 U_i |Z_i| <= 2|Z_i|.

The original coarse innovation has |D_W q|=sigma<=1, and every retained-caller/original-root row has norm at most one. The primitive endpoint now has no inverse-shield source-radius or caller-first cost.

Put

    C_* = 8! A^8 I_* / 8^6,
    T = C_* H_6(G) product_i J_i.

The expectation of T is EXACTLY the scalar 64-history packet restricted only in its deepest structural clock. This includes all eight primitive clocks, all other structural clocks, the common ancestry roots, and every private shield. Source-zero clocks of probability zero are never queried.

Uniformly in the retained caller, with clocks fixed under caller differentiation,

    ||T||_2 <= E_* := C_* sqrt(720),
    ||D_y T||_2 <= L_* := 16 sqrt(2) C_* sqrt(720).

For the energy, E[U_i^2 Z_i^4]=(1/3)*3=1, so the primitive-clock moments remove the earlier factor 3^4. For a differentiated factor use |D_y J_i|<=2sqrt(2)|Z_i| |D_y q_i|; each other factor still has squared envelope moment one. The latter displayed bound uses |D_y q_i|<=1. It is a genuine complete retained-caller first for this scalar source, not a derivative of a law estimate. It is also uniform in all primitive endpoint widths. Mixed correlated q_i cause no difficulty because the envelopes use only independent private Z_i and the Hermite score G.

In Sections 5–7, y is the original retained caller; it is not the intermediate outer parent Y used in Section 2. All seven structural roots are now internal sample roots. The sampled clocks are randomized frozen program parameters in these bounds. A native Gaussian-bank argument that also differentiates a Gaussian encoding of the clock choices requires additional clock-map rows and their cutoffs. That is not silently included in L_*. A direct coupling consumer below does not require a Gaussian Poincare estimate in the clock choices.

## 6. Explicit finite radius, original-root first, precision, and sampling bill

Choose an even Sobolev cutoff chi_R equal to one on [-R,R], zero outside [-R-1,R+1], in [0,1], with derivative magnitude <=2. Multiply H_6(G) by chi_R(G), and each J_i by chi_R(Z_i). Do not change the original affine query rows, callers, clocks, or source origins. Write M=R+1 and

    B_6=M^6+15M^4+45M^2+15,
    B_5=M^5+10M^3+15M.

The truncated actual coefficient has uniform radius

    B_* = C_* B_6 M^16.

A conservative complete first with respect to all original Gaussian query-bank roots and the eight private Z_i, with clocks frozen, is

    F_* = C_*[(2B_6+6B_5) M^16
              +16 sqrt(2) B_6 M^15
              +8 B_6 (2M^2+2M) M^14].

The three terms are respectively the score/cutoff derivative, all affine q_i caller/root derivatives (their individual row norms <=1), and all eight private-score/cutoff derivatives. Triangle bounds deliberately avoid an independence assumption about q_i. The retained caller alone has bound

    F_caller = 16 sqrt(2) C_* B_6 M^15.

Multiplying the target coefficient by an external amplitude lambda multiplies these bounds by |lambda|. For a stated radius gap b and first gap f it is sufficient to require

    |lambda| <= min(b/B_*, f/F_*).

If the eight-force amplitude is alpha^8, this gives the explicit scalar guard alpha <= min((b/B_*)^(1/8),(f/F_*)^(1/8)). This is a source/pushforward guard. It is NOT a proof of the imported native pair program's full side-spine, calibration, or gradient/curl contract.

Gaussian union bounds give a uniform conditional-mean cutoff bias

    |E T_R - E T| <= E_* sqrt(18) exp(-R^2/4).

Thus choosing R >= sqrt(4 log(max(1,E_*sqrt(18)/epsilon_bias))) makes it <=epsilon_bias (take R>=1 if desired). The source retains L2 <=E_* after cutting.

For N independent COMPLETE samples, average T_R with its same clock/ancestry construction inside each sample. Let K(y)=E T and Khat_N(y)=N^-1 sum T_R^(j)(y). Use the same sample tapes when evaluating Khat at different callers, so its actual caller derivatives exist. Then

    ||Khat_N-K||_(L2(Y,samples)) <= epsilon_bias + E_*/sqrt(N),

for any law of retained Y, since the bounds were uniform in Y. A sufficient finite choice for total epsilon is epsilon_bias<=epsilon/2 and N>=4E_*^2/epsilon^2. No oracle is hidden in that average.

Complete source bill: 16N original g VALUE evaluations, 23N original/private standard Gaussian draws (seven structural, eight coarse, eight private), 15N scalar clock draws, scalar inverse-CDF geometry for the deepest clock, and the displayed finite cutoff/polynomial/arithmetic work. No old call is assumed free merely because its site looks similar. Exact repeated anchors may be cached only with the same source/caller/precision version. All histories in a sample share their literal ancestry. Replication is across complete grouped samples, not across individual derivative histories.

If every g VALUE is returned with absolute error <=delta_g, the clipped importance-weighted factor has error <=2sqrt(2)M delta_g/A. For this to be <=1, the product telescope gives coefficient error at most

    C_* B_6 * 8 * [2sqrt(2)M delta_g/A] * (M^2+1)^7.

For A>0 choose delta_g to make this <= the requested arithmetic allowance. The A=0 source is identically zero. Scalar clock/root arithmetic needs its own propagated absolute error; near endpoint shields are not divided in J_i, so small sigma does not create a hidden VALUE-precision inverse width. The clock inverse CDF and affine roots must be evaluated to that separately chosen precision. This is an explicit finite VALUE model, not an assertion of an already available hardware precision oracle.

## 7. A direct positive consumer and its cost

In one physical dimension let X be Gaussian of variance v>0, independent of the coefficient sample bank, and let H_7^(v)(X)=v^(-7/2)H_7(X/sqrt(v)). Keep every existing retained caller/observer and an independent positive Gaussian buffer unchanged. Form the actual random pushforward

    W_N = base + [Khat_N(Y)/8!] H_7^(v)(X) + buffer.

This is a positive probability law for every finite N and R. Couple it with the identical expression using K(Y), with the same retained callers/observers, base, X and buffer. Provided those retained observers do not include integrated coefficient roots,

    W2(joint law W_N, joint law W_exact)
      <= sqrt(7!)/(8! v^(7/2))
                       [epsilon_bias+E_*/sqrt(N)].

If the coefficient carries an additional scalar amplitude lambda, multiply the right side by |lambda|. This is ordinary coupling of two actual random maps; it is not an expectation tensor sampler or formal cancellation. A requested law allowance eta can be met with the preceding cutoff choice and

    N >= [2 |lambda| sqrt(7!) E_*/(8! v^(7/2) eta)]^2.

A cutoff of the physical polynomial, if required by a bounded native prior, has its own explicit Gaussian-tail/caller-first bill and is not included for free here. The polynomial pushforward itself needs no boundedness assumption to be a positive finite-second-moment law.

This consumer certifies only approximation to the stated positive rank-eight reference map. The old all-rank raw-cumulant consumer has lower-rank moving-covariance and map-feedback currents; adding this coefficient does not cancel them. Nor does the Monte Carlo bill establish any desired asymptotic native query exponent. For fixed lambda=alpha^8 and target eta=alpha^P, direct variance reduction alone can cost order alpha^(-2(P-8)).

## 8. The remaining structural and multidimensional debts

### Structural endpoint

As r_* tends to one, I_* grows like [4(1-r_*^2)^2]^-1. The importance sampler merely normalizes this mass. It cannot turn the absolute rho^7(1-rho^2)^-3 kernel into a finite full-endpoint measure. The omitted parent-endpoint sector is not known small from this representation. It must be grouped with a larger ancestor cluster or handled by a different identity, with that sector's actual first and law budgets proved.

### Multidimensional one-energy requirement

The naive vector Hermite replacement is not dimension-free in Hilbert energy. Even with g(x)=a x, whose exact sixth common derivative is zero, the labeled score coefficient contains

    a^8 s^-6 I_2 tensor H_6(G),

before outer scalar weights. Its squared L2 Hilbert norm is

    a^16 s^-12 D * 6! * binomial(D+5,6),

of order D^7. Symmetrizing the eight physical slots does not reduce its leading dimension power: components with exactly one repeated coordinate and six other distinct coordinates already contribute at least

    a^16 s^-12 D(D-1)...(D-6)/28.

Thus the naive source has energy of order D^(7/2), while the actual cumulant has one-energy order sqrt(D). A constant-in-G subtraction kills this particular linear fixture, so this is NOT a no-go theorem for a more careful centered/composite source. It is a concrete reason the displayed scalar construction cannot simply be asserted to satisfy the multidimensional native one-energy ledger.

### Actual native admission and joint law

The finite VALUE estimator/direct positive consumer is not the imported C_k selected-pair/filter compiler. The radius and first estimates above are for the actual new scalar source; they do not certify native selected matrices, physical readout shares, full chronological feedback, or the terminal gradient/curl condition. If a native compiler is used to execute it, its complete returns and all integrated/clock-bank observers must be regenerated, rather than borrowed from the old per-history packet. In particular a same-mean IBP replacement must never be called an exact joint-law replacement of its old random coefficient.

## 9. Extension that is proved, and extension that is not

At rank n, append n-2 fresh leaves to the original siblings and retain both sibling choices at every stage. The exact packet has 2^(n-2) labelled histories, coefficient n! per history, binomial multiplicities, all n-1 structural resolvents R_n, and the common derivative order n-2. On a fixed interior deepest-parent sector the same scalar construction uses H_(n-2), 2n original g VALUE calls per complete sample, and joint-clock normalized energy bound sqrt((n-2)!). The scalar prefactor is n! A^n I_n/n^(n-2), with the structural normalizer I_n displayed below. Primitive importance weights remove every primitive inverse shield just as above. The deepest structural normalizer becomes integral rho^(n-1)(1-rho^2)^(-(n-2)/2)d rho and ceases to be absolutely integrable at the endpoint once n>=4.

This proves a computable finite family of cancellation-preserving original-VALUE conditional-mean producers. It is not yet a producer for every rank-n history, and not an all-rank native positive algorithm. The next missing result is specifically an ancestor-cluster regrouping that removes the transferred structural endpoint debt while retaining one Hilbert energy and the complete native/observer first census.

## 10. Checks and provenance

Input main pin:

    /workspace/shared/arbitrary-ancestry-heat-20261005/FINITE-ANCESTRY-COMPILER-AND-BOUNDARY.md
    SHA256 e758d05d9901093701bd6962477a6ffb967860ca534e7742624e19de6f2c50ef

Generator pin:

    /workspace/shared/rank-indexed-positive-returns-20261005/generate_rank_terms.py
    SHA256 0fbff5a4d57a0d89ff84373422d72913e759061156faa09f4f0afcc6c564c74f

`check_sibling_regrouping.py` checks literal generated memberships and resolvent indices through rank ten; exact Leibniz identities through order eight; nontrivial polynomial Gaussian-VALUE/IBP identities through order eight; and the exact Gaussian cosine covariance asymptotics. It writes the actual 64 rank-eight ASTs to `rank8_packet.json`. Polynomial fixtures test algebra only, not bounded-Hessian source admissibility. The oscillatory admissible source supplies that separate analytic test.

Independent audit: `/workspace/shared/sibling-regrouping-audit-20261005/`. It independently reconstructed the 64 ASTs, verified clock powers, oscillatory modes/variance, and the original VALUE/first identities. Its Section 9 and separate `check_importance.py` / `importance-checks.json` also independently verify the importance-sampling extension and sharper joint-clock constants.

`finite_scalar_value_producer.py` executes the finite estimator using only a g VALUE callback. `producer-runs.json` records 4096-sample runs for a linear null packet and bounded-Hessian oscillatory sources K=4,16,64,256. These runs validate executable bookkeeping, not a deterministic realization error claim. The rigorous stochastic certificate is the stated RMS bound. Floating clock geometry has a separately disclosed precision debt. The implementation does not execute imported native source/pair programs.
