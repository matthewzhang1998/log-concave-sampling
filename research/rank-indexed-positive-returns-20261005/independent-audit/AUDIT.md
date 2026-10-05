# Independent audit: all-rank Appell bounds and positive polynomial consumer

2026-10-05. Scope: analytical finite-rank identities, dimension-free tensor bounds, and completeness of the positive-consumer currents in `ALL-RANK-APPELL-AND-POSITIVE-CONSUMER.md`. No native original-VALUE compiler or source-admission conclusion is asserted.

## Verdict

**PASS within the stated analytical scope.** No mathematical correction to A1–A8 was required. In particular:

- A1's covariance-operator induction is dimension safe for mixed, dependent Gaussian-bank vector sources.
- The one-mark inequality and telescoped all-rank stability bound have the stated factorial constants.
- The proper-cut recursion uses genuine cross-partition identities, and never obtains a self-trace from an abstract operator bound.
- A6/A7 include every Faà di Bruno partition and both covariance-path contributions.
- The surviving-physical-output graph argument proves a one-energy bound even when Wick pairings create cycles.
- A8 includes the moving-root/time term, covariance sampler term, all polynomial-map feedback, and the terminal centered Stein current.
- The very large displayed `C_m` is a valid conservative choice. An explicit counting domination is supplied below.

The important limitation is correctly stated: rank extension improves the Stein remainder, but the raw-cumulant polynomial reference retains `L² h₃`, so this particular uncorrected consumer does not establish arbitrary-order accuracy.

## 1. Covariance induction and one-mark factorization

Let `a_m` denote the covariance-operator bound for a labelled rank-m Appell tensor. For a deterministic tensor `T`, write the derivative of its contraction with `W_I` as the sum of m terms. In the j-th term, first form the vector `v_j` by contracting all slots other than j against `W_(I\j)`. Then the term is `(DB_j)^* v_j`; pointwise operator control gives

    ||(DB_j)^* v_j||² <= L_j² ||v_j||².

No independence between `DB_j` and the other Appell factors is used. Slice `T` along its j-th physical slot. The covariance bound for `W_(I\j)` gives an expectation bounded by `a_(m-1)` times the sum of squared slice norms. That sum is exactly `||T||HS²`, with no dimension multiplier.

Gaussian Poincaré and the elementary inequality `||sum_j u_j||² <= m sum_j ||u_j||²` yield

    a_m <= m² a_(m-1) L²,    a_1 <= L².

Thus `a_m <= (m!)² L^(2m)`. With unequal bounds, the same calculation gives `(m!)² product_j L_j²`.

For a centered vector `R`, the matrix `E[R tensor W_I]` is a cross-covariance operator. Cross-covariance factorization, or coordinatewise duality against its rows, gives

    ||E[R tensor W_I]||HS² <= E||R||² ||Cov(W_I)||op.

The first factor is the sole physical energy. Taking square roots proves A2. A finite moment-cumulant expansion proves the identity with the one-mark cumulant even when R has no exponential moments: Gaussian Lipschitz factors have moments of all finite orders, and Cauchy–Schwarz supplies every required mixed moment with `R in L²`.

For A3, telescope n slots and apply A2 with `m=n-1` to each term. This gives precisely `n (n-1)! = n!`. The remaining factors may be arbitrary mixtures of U and V on their original common bank. No independent recopy is used or needed.

A5 is also an immediate special case of A2 with `R=B` and `m=n-1`; its proof does not actually require a separate appeal to the centered Stein recurrence.

## 2. All-proper-cut recursion

For a fixed p|q cut, normalized exponential generating functions give

    E[W_p tensor W_q]
      = sum over partitions whose every block crosses the p|q cut
          tensor products of the corresponding cumulants.

The one-block partition contributes `κ_(p+q)`. Every other participating block has at least one slot on each side. Because a second block also needs a slot on each side, every such block has size at most `p+q-2`. The induction is therefore genuinely lower rank.

Cross-covariance factorization yields

    ||E[W_p tensor W_q]||_(p|q op) <= p! q! L^(p+q).

Each proper cross-partition product is, after separate input and output permutations, a tensor product of proper lower-cumulant flattenings. Its operator norm is at most the product of their bounds. The triangle inequality proves exactly the stated recursion.

This reasoning is dimension free. It does not use a cut estimate to evaluate a self-trace, nor does it apply the recursion to arbitrary abstract tensors.

An independent rooted bipartite-block dynamic program, rather than the author's partition generator, gives:

    n:    2  3  4   5    6    7     8      9       10
    c_n:  1  2  6  24  132  864  7248  67920  760320

The author's available values through rank six agree.

## 3. Positive interpolation and full current census

The covariance path is `C_t = I/2 + (1-t²)Σ` at unit buffer, so it commutes with Σ and

    x'_t = -t Σ C_t^(-1) x_t.

The direct derivative of the map splits into:

1. the explicit coefficient derivative `-sum_r r t^(r-1) q_r`;
2. `D_t`, which includes covariance-score differentiation and the changing visible root;
3. the standalone derivative of `x_t`.

Gaussian integration by parts on item 3 yields `-t Σ` against the Hessian, plus `-t Σ(Dp_t)^T`. The former cancels the source Stein covariance term; the latter is the A2 covariance sampler current.

Integrating the score defining q_r leaves the complete derivative of order r−1 of the composite first derivative of φ. The multivariate chain rule has exactly one term for every labelled partition. A block of size d supplies `D^d(x+p_t)` and therefore emits one physical derivative output, regardless of d. This proves A6 with current rank `1+|pi|`, including the largest block.

The source-bank Stein expansion applies at the same actual W_t because the map and both Gaussian roots are independent of that source bank. Its constant rank-r coefficient is `κ_r/(r−1)! = r K_r`. This cancels exactly one term: the all-singleton term with all its H factors equal to identity. The rest remain as written in A7.

### Rank-five census, explicitly

In scalar notation, put `J=p'`, `H=1+J`, `E=p''`, `F=p'''`, and `G=p''''`. The complete six currents are

    A1 = D_t,

    A2 = -t ΣJ -3t² K3 E -4t³ K4 F -5t⁴ K5 G,

    A3 = -3t² K3(H²-1) -12t³ K4 EH
         -5t⁴ K5(4FH+3E²),

    A4 = -4t³ K4(H³-1) -30t⁴ K5 EH²,

    A5 = -5t⁴ K5(H⁴-1),

    A6 = t⁵ T5.

In multiple dimensions, the products denote the corresponding labelled tensor contractions; identical scalar profiles group permutations that need not be literally identical unsymmetrized tensors. The rank-five partition multiplicities are 1, 6, 4, 3, 1 for profiles `1111`, `211`, `31`, `22`, `4`. In particular both the fourth map derivative G and the double-Hessian term E² must remain. They are present in A6/A7.

## 4. Why the graph proof survives Wick cycles

After expanding every H into identity plus map derivatives, choose the original K_r as anchor. Each added K_s factor comes from a derivative of order d≥1. Its d chain-rule inputs attach directly to the anchor. It has one final physical output and `s−1−d` Gaussian slots. Its physical output is never a Gaussian slot.

Use the multivariate Wick product formula on the Gaussian slots. A Wick pairing always joins two different initial factors. The transformed coefficient graph consequently has:

- a star of chain-rule connections to the anchor;
- additional arbitrary cross-factor Wick edges;
- at least one permanently open physical output at every coefficient vertex.

Start with the anchor in Hilbert–Schmidt norm. Add vertices in any order. When adding a vertex, contract all its edges to the already built group at once. There is at least one such edge, because the anchor is already present. There is at least one uncontracted slot, its final physical output. The added tensor is therefore used through a genuine proper cut. The standard matrix product inequality gives

    ||new partial contraction||HS
       <= ||old partial contraction||HS × ||new coefficient||_(proper-cut op).

Every edge is contracted precisely when its later endpoint is added. Thus no later self-trace of a previously bounded aggregate is required, including when the graph has Wick cycles or multiple edges.

After the deterministic contractions, the remaining Gaussian slots contract a Hermite tensor. Wick isometry supplies `sqrt(N!)`, because symmetrization on Gaussian slots is an orthogonal contraction. This proves one anchor Hilbert scale and one proper-cut scale for every added coefficient.

The argument would fail for arbitrary graphs that consume a coefficient's last output. The consumer graphs do not do that; the restriction is substantive and correctly stated.

## 5. Explicit constant domination and A8

Write `N=(m−1)(m−2) < m²` and `J_m=2^m m!`.

- There are at most `m^m` labelled partitions at each anchor rank.
- Each partition has at most m−1 nonidentity coefficient factors.
- Their derivative/covariance transforms contribute at most `J_m^(m−1)`.
- At most N Gaussian slots occur. Their partial matchings inject into the permutations of N labels, as involutions with fixed points. Thus at most `N!` diagrams occur. Restrictions against within-factor pairings only reduce this number.
- Hermite contraction costs at most `sqrt(N!)`.
- Transfer of a rank-ell current to the untouched unit-buffer half costs at most `2^(m/2) sqrt(m!)` for `ell<=m+1`.
- The factors r, grouping by number of added vertices, and any harmless rank-choice overcount fit inside a factor `m^(2m+3)`.

Consequently all feedback constants fit below the deliberately loose quantity

    D_m = m^(2m+3) 2^(m²+m) (m!)^(m+1) ((m²)!)².

The separate D_t/covariance constants, and the remainder transfer constant, are smaller than this same bound for m≥3. For example, differentiation of the standardized covariance-Hermite transform gives `(r−1)` choices, one matrix of norm at most `2 sqrt(2) L²`, the remaining matrices of norm at most `sqrt(2)`, and the finite Hermite norm. The covariance sampler uses `||Σ||op<=L²` and the first-derivative version of the same bound.

For m≥3,

    m^(2m+3) <= (2m²)!,
    (m!)^(m+1) <= (m(m+1))! <= (2m²)!,
    ((m²)!)² <= ((2m²)!)²,
    2^(m²+m) <= 2^(10m²).

The first inequality follows already from the last m² factors of `(2m²)!`, whose product is at least `m^(2m²)`. Thus `D_m <= C_m` for the displayed `C_m = 2^(10m²) ((2m²)!)^4`.

Each nonleading feedback monomial has between one and r−1 added coefficient factors. Summing their rank choices gives `h_r sum_(j=1)^(r−1) s^j`, up to the preceding fixed-rank constant. The two covariance/time terms give `L²h`; the terminal source current gives `L^m e`. These are precisely the three groups in A8.

All coefficients are independent of Z2. Gaussian integration by parts on Z2 leaves a first-derivative current, and conditional expectation given W_t can only decrease its L² norm. The resulting integrable velocity gives the usual dynamic Wasserstein bound. Polynomial moments, the uniformly positive keep covariance, and finite-rank Sobolev approximation justify passage from smooth truncations. Neither positivity nor a global Lipschitz bound on the polynomial reference map is required beyond its finite second moment and the exhibited velocity bound.

The general-v rescaling is correct: physical output variables scale by `sqrt(v)`, while rank-r cumulants scale by `v^(r/2)`. The stated small-radius grades follow from `s=O_m(L³)`, `h=O_m(L²e)`, the terminal `L^m e`, and the surviving `L²h₃=O(L⁴e)`. The added explicit `C′_m=C_m(3+2M_m)`, with `M_m=sum_(r=3)^m c_r/r!`, is safe: for `L<=1/2`, `h<=(2/3)L²e` and `s<=M_m L³<=1/2`, so the feedback is at most `2hs`. These bounds are strictly dominated by the stated enlarged constant for both m=3 and m>=4.

## 6. Independent exact checks and their limits

`check_appell_consumer.py` imports no author generator. It uses exact integer/rational/symbolic arithmetic and records 2,660 successful assertions in `independent_checks.json`.

The checks include:

- rank-1–8 Appell derivatives and centering;
- all p,q≤4 cross-partition identities;
- one-mark moment-cumulant identities and bounds;
- stability for same-bank scaled bounded sine sources;
- proper-cut constants through rank ten from a separate dynamic program;
- exact Faà di Bruno polynomial identities through rank seven;
- full rank-five A7 tests at t=1/3 and t=2/3, for φ(w)=w^n, n=1,…,7;
- the scalar reference variance `1+2a²+6b²+24c²`;
- 722 expanded rank-five feedback monomials, their permanent physical outputs, and maximum Gaussian degree 12.

The admissible scalar inequality fixture is `X=sin(sqrt(log(2)) G)`, with exact rational even moments. `log(2)>2/3` supplies rational lower bounds on the claimed right sides, so these checks are exact inequalities, not floating-point sampling.

The full current identity fixture is `B=G+(G²−1)/4`. Its cumulants are nonzero, and the finite Stein recurrence is evaluated exactly in the Gaussian Hermite basis. This fixture is not globally Lipschitz; it tests the identities only and is not offered as a bound-admissibility example. The left side is computed by differentiation of the Gaussian density while keeping physical x fixed, independently of the fixed-root derivation on the right side.

None of these scalar or finite-rank checks proves the general theorem. The dimension-free proof is the induction, covariance factorization, cross-partition subtraction, and graph argument above.

## 7. Source and execution boundary

The imported centered Stein recurrence is consistent with the normalization used here: `E T_j=κ_(j+1)/j!`, `||T_j||_(L²;HS)<=L^j e`. The fourth-cumulant predecessor's five currents are recovered exactly by setting m=4 in the generator.

All conclusions are analytical. They provide a repeatable finite law-consumer target with computable constants. They do not realize K_r through native original-gradient VALUE queries, do not license differentiation of completed generated sources, and do not resolve the observer-boundary/cycle/current closure requirements. They imply no sublinear coefficient-cost claim.
