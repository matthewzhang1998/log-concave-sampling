# A typed repair algebra, its residual generator, and the exact closure boundary

2026-10-05. A formal construction from the existing audited source/current interfaces. No native sampler was executed. No new all-order unsmoothed theorem is claimed.

## 1. Answer in one paragraph

There is a useful algebraic structure: a graded algebra of typed marked graphs, equipped with derivative insertion, Gaussian contraction, connected-cluster extraction, bank-tree gluing, and heat-source replacement. The rank generators are its local instructions. The **meta-generator** is the full residual expansion of their common-endpoint compiler. At a fixed heat floor, if every emitted current has an admitted leading-current right inverse and **every old and new port that can be hit has a uniform positive reserve**, this residual operator strictly raises a filtration. It can then be inverted by finite formal substitution at every requested cutoff. That is a genuine repeatable construction. The failure to extend it to the original unsmoothed target is concentrated in explicit escape operations: heat traces or insufficiently improving shells, an exposed source-zero bank observer, and any residual without a certified source/current port. Merely calling the structure a group does not fill those ports.

## 2. The objects must remember more than rank

A connected generator is a labeled port graph with the following record:

- A finite set of original force occurrences. Occurrence v has primitive order k_v, hence k_v+2 tensor slots. Identical force labels at different occurrences are NOT identified.
- A matching of internal slots into edges, a list of physical external slots, and typed open bank/probe/current ports. A loop and a multiple edge are distinguished.
- Actual affine caller ancestry and coefficient-bank names, original structural and primitive clocks, private shields t_v, source/anchor versions, and all injection rows.
- A marked-spanning-tree certificate when one exists, a choice of genuine root, and the source-qualified selected-spine/current certificate. Graph cut admissibility alone is not native admissibility.
- Ownership: retained external record E, primitive coefficient bank B, owned dynamic tapes U, protected physical marks, source-zero carrier rows, one common physical reference H, and the untouched keep.
- Scalar coefficient, nominal amplitude a, source allocations beta_v, width exponents gamma_v, response normalizations, first-path certificates, covariance gaps, error floors, and a versioned evaluation-DAG cost record.

An unresolved term can instead have a declared escape type, such as HEAT_TRACE, RETAINED_OBSERVER, MISSING_CURRENT_PORT, NONPOSITIVE_RESERVE, or UNPAID_CLOCK. Keeping that term is mandatory; an escape is not zero.

For a completed coefficient graph, with N force occurrences, M internal edges and R physical legs,

    K = sum_v k_v,                 K+2N = 2M+R.

Open bank derivative ports add their own count on the right until contracted. Physical rank R is not the grading that makes repair terminate.

The same graph can have several analytic realizations with different source allocations, clocks, callers or ownership. Those are different typed generators. A change of allocation or heat representation requires the actual source program to be rebuilt and its affected calls replayed.

## 3. The algebra and the genuinely finite instruction set

Fix a finite original force-label dictionary and a finite source-version/clock list for a requested computation. Let A be the real vector space of all such formal graph expressions, including escape types. Addition means a labeled sum with exact coefficients; multiplication is disjoint union with names renamed apart except for explicitly shared external bank names. Quotient only by proved slot-permutation and multilinearity identities. Do not quotient away caller ancestry or Gaussian ownership.

The underlying commutative algebra is the symmetric algebra on connected typed graphs. It is infinite-dimensional. More strongly, it generally has infinitely many connected algebra generators: disjoint-union multiplication cannot manufacture a new connected graph. What is finite is the following collection of **operation schemes**, not a finite list of tensor ranks.

1. **Anchored VALUE source and signed difference.** Start with g(z+t x)-g(z). Apply a signed two-point difference in one fresh probe, affine rescaling, and known scalar normalization. Repetition constructs every C_k source. More precisely, set a_k=(k+1)^(-1/2) and Delta_p h(z)=[h(z+a_k t p)-h(z-a_k t p)]/2. Applying k commuting Delta operators to g(z+a_k t x)-g(z), then dividing by A t, is exactly the supplied active source. Choose the common a_k after k is known. This is a finite instruction schema with unbounded input k; it uses 2^(k+1) original VALUES before literal duplicate caching. The common radial pullback is another fixed source-preparation instruction.
2. **Relabel, sum, tensor product, and port permutation.** Their scalar and incidence effects are literal. Shared-bank equality is retained.
3. **Bank differentiation.** Apply the product rule at every occurrence and to every actual caller/normalization that depends on that bank. For a heated force coefficient a spatial hit promotes k to k+1 and adds one bank slot with the actual injection row. Hits on scalar geometry, covariance or heat increments are separate branches. A source-level derivative label is not a derivative oracle call.
4. **Gaussian pairing/tilt.** For every unprotected eligible Gaussian slot, choose its physical tilt leg or its Wick pairing, with labels and covariance retained. Each choice is made once in a complete finite expansion. A pair consumes two external slots and adds an internal edge. Pairing a protected mark or producing an unadmitted self-loop emits an escape rather than silently entering the native class.
5. **Connected extraction.** Convert moments to cumulants by the exact set-partition formula below. For the source-qualified conditional program, retain all root-copy, side-occurrence, public, auxiliary and finite-mean branches, not only full copies of the selected main graph.
6. **Whole-bank bridge.** Select a labeled tree on the coefficient arguments, differentiate once at each edge endpoint, enumerate every occurrence hit, and contract the two new bank ports. Retain the positive replica covariance, all original clocks and injection rows. Equality holds after the complete bank expectation, not pointwise in the retained old bank.
7. **Heat replacement.** A cold/warm source difference is one force occurrence with a branch label, the common private variables and the same ancestry. A sequential telescope replaces one occurrence at a time. A heat-variance derivative is a different rule: it inserts two same-occurrence derivative slots and contracts them. The latter creates a self-loop and is not admitted by the present graph compiler.
8. **Common-endpoint current composition.** Apply the finite chain/product rules, covariance-root identities, visible Hermite stripping, and the supplied native finite-current jet table at one endpoint. A positive simultaneous compiler uses one common reference H and charges all cross-group branches. This rule has an imported semantic contract; the graph algebra does not prove that contract.

An operation taking an arbitrary finite list or integer is still one operation scheme. Calling this a finite *group generating set* would conflate that fact with finite generation under multiplication.

### A modest Hopf structure, with no extra analytical claims

For each connected graph C define Delta(C)=C tensor 1+1 tensor C, and extend Delta multiplicatively to disjoint unions. Define epsilon(1)=1, epsilon(nonempty)=0, and S(C_1...C_m)=(-1)^m C_1...C_m. Direct expansion proves coassociativity and the antipode identity. This elementary connected-component Hopf algebra records disjoint-union factorization and its formal connected bookkeeping. Actual expectation of two graphs sharing an external coefficient bank is generally NOT multiplicative, so that expectation is not automatically a character of this algebra. Its physical cumulants must still be computed by the explicit partition/whole-bank rules below; the Hopf construction does not assert otherwise.

It is NOT a claim that the physically allowed graph rewrites are generated by that coproduct, nor a substitute for the bank-tree identity. In particular its antipode neither makes a signed density positive nor pays an inverse width. Spatial bank derivatives are commuting derivations; their square is generally nonzero. Thus calling this a differential graded algebra without introducing a separate square-zero differential would be inaccurate.

## 4. A complete residual meta-generator

The universal scalar connected instruction is

    kappa(X_1,...,X_m)
      = sum_partitions pi (-1)^(|pi|-1)(|pi|-1)!
                             product_(C in pi) E product_(i in C) X_i.

A finite conditional Gaussian moment is obtained by summing every eligible labeled pairing after choosing physical shifts. This is an algorithm with rational/integer coefficients, not an asymptotic analogy. Use literal repeated labels and repeated clocks; do not substitute independent histories.

For a fixed finite program, the meta-generator ExpandResidual takes:

- its normalized intended current;
- its supplied finite primitive current-jet identities;
- its declared retained coordinates, common endpoint and positive reference;
- a finite extraction order and exact allocation/width/precision records.

It performs the following deterministic steps.

1. Expand the entire requested family together, including all old/new products. Expand the finite source/filter/native jets to the supplied order. Do not Taylor expand an unsmoothed original C2 force to arbitrary order; higher force labels are the privately heated coefficient/source interfaces.
2. Enumerate all product-rule hit locations, all chain-rule set partitions, covariance-root/time/readout branches and all scalar normalization derivatives. Where a primitive current table is unavailable, emit MISSING_CURRENT_PORT with the literal unevaluated expression.
3. Enumerate complete conditional pairings and physical tilts, then perform the partition sum selecting connected coefficients. Keep lower physical ranks, finite means and covariance/normalization currents.
4. For every whole-bank cumulant, enumerate labeled replica trees and all repeated endpoint hits, with the original shared-bank dictionary. A new bank edge creates derivative slots, not new clock mass.
5. At the same physical endpoint, subtract only the declared intended main current and cancellations justified by exact labeled identities. If two terms have distinct caller/ownership records, they cannot be canceled merely because their averaged tensors look alike.
6. Type-check every surviving term. Return an admitted current with all ledgers, a bounded above-cutoff remainder with its analytical certificate, or an explicit escape term. Return the numerical/native floors separately even when the exact coefficient is zero.

### Formal soundness and completeness

Fix the primitive identities, finite jet order, finite clock list and finite source program. Assume every primitive identity is exact to that order with a separately certified remainder. Then the output of ExpandResidual is exactly the complete formal current residual through that order, plus those primitive remainders.

Proof: structural induction on the program syntax. The base cases are the supplied primitive identities. Addition, product and composition use the exhaustive labeled product and set-partition chain rules. Gaussian expectations use the exhaustive matching identity. Cumulants use the displayed finite partition sum. The exact whole-bank tree identity expands each connected coefficient, and its product rules enumerate every possible hit. Subtraction removes only exact matched expressions. Every possible branch occurs in the relevant finite choice list; an untypable branch is kept as an escape. These observations prove both equality and completeness by induction.

A primitive remainder known only by a norm bound can be kept as a boundary error when its actual bound meets the requested cutoff. If it does not, the generator emits MISSING_CURRENT_PORT; no formal term or source right inverse is inferred from that bound. This includes an intrinsic quadratic native-to-ideal floor, even when the ideal polynomial cross-jet is fully known.

The theorem concerns the **formal current expansion supplied by those primitive identities**. It is not a proof that the existing finite native comparison library already has a table for every future program. Nor is it a convergence theorem for an infinite moment-generating series. At finite cutoff, no such convergence is needed.

## 5. The filtration that actually works at a fixed floor

Take common gamma>0, beta>0 and fixed private width floor t_v>=c_v alpha^gamma. Define

    G=a-gamma K,
    d=a-beta N-gamma K,
    Psi=d-2gamma.

The source allocation is alpha^(beta+gamma k_v)t_v^(-k_v) at every occurrence and the extra alpha^d times the exact signed scalar coefficient at the genuine root. It reproduces the intended coefficient exactly. The effective root radius has grade beta+d.

The two essential rewrite identities are:

**A connected r-argument bank tree:**

    a'=sum a_i,  N'=sum N_i,
    K'=sum K_i+2(r-1),  M'=sum M_i+r-1,  R'=sum R_i,
    Psi'=sum Psi_i.

**A conditional nonmain offspring with h>=2 genuine root copies:**

    a'=beta N'+gamma K'+h d,
    d'=h d,                Psi'=h Psi+2gamma(h-1).

The second identity uses the actual retained side occurrences; it does not assume a full copy of every side. It imports the audited source-qualified conditional grammar. Physical Wick/tilt choices leave a,N,K and Psi unchanged and form a finite local enumeration. Hermite stripping decreases its finite visible degree q; it is not an outer repair step.

For heterogeneous widths fix Gamma>=sup_v gamma_v in advance, require beta_min=inf_v beta_v>0 (as in a fixed finite positive label dictionary), put H=sum gamma_v k_v and B_beta=sum beta_v, and use

    Psi=a-B_beta-H-2Gamma.

The derivative width charge of each inter-packet bridge is at most 2Gamma. More exactly, if an r-argument tree adds width charge DeltaH and a further net inverse-alpha injection/normalization loss ell, then a_child=sum_i a_i-ell and

    Psi_child=sum_i Psi_i+2Gamma(r-1)-DeltaH-ell.

The superadditive rule Psi_child>=sum_i Psi_i is admitted only when this FULLY charged right side is at least the parent sum, or the loss was already honestly prepaid in the parent grades. Without that check the branch emits a missing-gain/NONPOSITIVE_RESERVE escape even if its absolute reserve remains positive. No injection loss is hidden in a constant. For h identical root types the conditional formula becomes h Psi+2Gamma(h-1). For h potentially different root types, the exact allocation instead gives d_child=sum_i d_i and Psi_child=sum_i Psi_i+2Gamma(h-1); every retained occurrence keeps its own beta, gamma and source label. A mark receives no free smallness credit: a derivative of its increment row may remove the mark.

### Resource-tag extension

The same proof permits an occurrence charge q_v(k) instead of k if a **separately proved owned resource** pays the removed width factors and 0<=q_v(k+1)-q_v(k)<=1. Set H=sum gamma_v q_v(k_v). A bank hit then still costs at most gamma_v. The proposed owned-clock choice q(k)=(k-2)_+ fits the algebra once its node inequality and exact scalar-amplitude relocation are certified. It does not arise from relabeling k. Each clock weight is attached to one original occurrence, retained or copied only when that occurrence is legitimately retained or copied; derivative insertion never creates a fresh clock weight. The independently audited weighted-generator contract in induction-renewable-heat-budget-20261005 validates this ownership-qualified extension and a second finite heat reentry ending at grade 41/5. It does not validate indefinitely renewable clock credit. The present proof uses only the stated per-occurrence charge inequality.

## 6. A repeatable formal inverse theorem

This section is a formal coefficient algorithm. It is unrelated to time evolution, trajectory iteration, or any paused numerical programme.

Let J_ad be a typed space of admitted currents at a fixed heat/ownership boundary. Let V_ad be admitted source programs. Let L:V_ad->J_ad be their leading-current map. For **every emitted current type** used below assume an explicit right inverse h with Lh=Id; h includes exact selected/readout normalization and the required source/caller versions. The chosen h must preserve the declared filtration or enter every allocation, width and resource loss into the realized type BEFORE testing the reserve inequalities. This is not a reserve computed from the target tensor alone. If a type has no such h, the algorithm stops at that escape rather than assuming it away.

Freeze the finite old collection o, and normalize the full simultaneous common-endpoint formal compiler so that

    F_o(x)=x+R_o(x),                R_o(0)=0.

Any old residual independent of x is included in the target t to be canceled. R_o is the complete residual of compiling h(x) in the old context, not the sum of isolated per-packet errors. Its primitive current table must include all mixed old/new branches, moving covariance, normalization, auxiliary and retained-observer branches.

Assume:

A. Every new seed and **every old typed port actually hit in R_o** has Psi>=delta>0. An old O(alpha) first estimate alone does not establish this. A source-zero bank observer has to satisfy this condition or be excluded by ownership.

B. Every nonlinear residual branch is generated by the admitted conditional and bank-tree operations above, or another explicitly certified operation with the same positive-increment property. No linear uncharged bank/observer/heat branch is omitted.

C. The chosen formal compiler is locally finite in its reserve filtration and the right inverse h is available on all its emitted types. This can be supplied by a fixed finite primitive-jet/side-occurrence census with finitely branching chosen source realizations, or assumed and proved separately for a larger grammar. The effective-grade stopping argument below alone does NOT prove reserve-local finiteness of the universal graph algebra. All finite analytic/source floors are carried separately.

Then, whenever x-y has reserve at least p, every term of R_o(x)-R_o(y) has reserve at least p+delta. To prove this, expand one residual monomial and telescope its x occurrences: each difference term has one distinguished x-y occurrence. A mixed bank branch contains at least one additional argument, which may be old or new; by A it contributes at least delta. A conditional branch has h>=2 root copies, possibly of different types; its exact reserve is the sum of the root-type reserves plus 2Gamma(h-1). The other copies therefore contribute at least delta, and that extra term is nonnegative. Compositions add these nonnegative gains. All old occurrences remain in the telescope with their actual reserve. Finite Wick/tilt and Hermite extraction choices do not consume that gain. Assumption B is exactly what excludes an uncharged exceptional linear branch.

It follows that the equation

    F_o(x)=t,  equivalently x=t-R_o(x),

has a unique solution in the completed admitted filtration, determined by

    x_0=0,                  x_(m+1)=t-R_o(x_m).

Indeed successive differences gain delta at each substitution. Two solutions would have a first nonzero filtered difference, but the same equation would raise its valuation, a contradiction. At every requested finite cutoff the substitution stabilizes after finitely many steps. This is a precise meta-generator for the whole correction family, including feedback of previous repairs, rather than a prescription to add the next physical rank.

This existence statement is conditional on A-C and the supplied compiler identities. Current results give a substantial fixed-heated subgrammar satisfying the combinatorial recurrences; they do not establish A-C for every residual of the original unsmoothed task.

### Local finiteness and an implementable truncation

There are two distinct constructions. The inverse theorem solves an explicitly chosen reserve-locally-finite formal compiler. The operational correction queue below stops at analytically bounded high-G residuals and does not compile them into new sources. That queue is an analytically controlled restricted compiler, not automatically the G-truncation of the full formal inverse: G can decrease when a conditional offspring omits sides. In particular, arbitrarily many side occurrences may share one fixed root surplus/Psi in the universal grammar unless an actual finite side/jet census is specified. No universal pure-Psi local-finiteness claim is made.

For the operational queue retain G<P. Since G=beta N+Psi+2gamma,

    N<(P-2gamma)/beta,           Psi<P-2gamma.

For the conditional/bridge grammar just described, each repair lineage gains at least delta, so its length is bounded by a computable multiple of P/delta. Bank-tree arity and conditional root multiplicity are bounded for the same reason. Conditional copies do not increase k_v; one bridge stage increases a tracked k_v by at most its bounded number of incident edges. A conservative bound is

    k_max <= k_max,seed + ceil(P/delta)^2.

For this core grammar, valences, physical ranks, finite Gaussian polynomial degrees, labeled hit choices, source orders and graph sizes are therefore bounded. An extra operation allowed by hypothesis B must supply its own finite branching and bounded order-growth rule, or a separate local-finiteness certificate under hypothesis C; the displayed k_max formula is not asserted for arbitrary added operations. With a finite clock/source dictionary, every local expansion is finite, and the entire truncated expression tree is finite. Histories can grow enormously; finiteness is not a favorable complexity estimate.

A deterministic algorithm is: enumerate in increasing reserve; complete the finite visible-degree extraction for each term; merge only identical complete typed keys; solve its leading-current equation with h; emit all residuals; return escape terms verbatim; and bound discarded G>=P terms using their actual one-energy/native/keep certificates. A high-G cutoff is an analytical remainder operation, not deletion in the untruncated algebra or a rewrite-stable quotient. Boundary-current bounds must hold for the actual completed simultaneous program and endpoint; later changes require the full cross-current/first/replay ledger to be recomputed. Descendants of an unexecuted high-G correction source need not be executed, because that correction source was never added.

For heterogeneous allocations, use beta_min in the displayed N bound and Gamma in the bridge offset. For resource-tag charges q=(k-2)_+, k_max remains bounded by the same lineage argument. Bounding N alone would not bound K without that argument.

## 7. A minimal worked conditional example

Use a scalar selected-spine symbol

    V(theta)=lambda theta^2 (P+a theta)(Q+b theta),

where P,Q are independent standard Gaussians. The two permanent endpoint marks supply theta^2. The a and b coefficients retain their actual public readout/side-source grades; they need not have the same grade. This is a finite coefficient-identification example, not an executed unbounded polynomial sampler.

Its first two connected terms are

    E V = lambda a b theta^4,
    kappa_2(V)/2 = lambda^2 theta^4 [1+(a^2+b^2)theta^2]/2.

Thus repairing the rank-four main term immediately creates a rank-four and two labeled rank-six branches. A rule that adds only a higher physical rank misses a same-rank residual. Keeping a and b separate avoids an unjustified factor of two.

The complete formal logarithm is particularly transparent. For t=lambda theta^2, c=a theta, d=b theta, a direct two-Gaussian integral gives

    log E exp[t(P+c)(Q+d)]
      = -log(1-t^2)/2 + [t c d+t^2(c^2+d^2)/2]/(1-t^2).

Therefore, through fourth spine-amplitude order,

    log E exp V
      = lambda ab theta^4
        + lambda^2 theta^4/2
        + lambda^2(a^2+b^2)theta^6/2
        + lambda^3 ab theta^8
        + lambda^4 theta^8/4
        + lambda^4(a^2+b^2)theta^10/2 + O(lambda^5).

The set-partition/Wick meta-generator produces exactly this list, including the lower physical ranks, signs and multiplicities. One does not execute the logarithm as a density or assume its global MGF exists for arbitrary theta; the finite jet identity is sufficient.

For a concrete bridge grade, take two admitted three-vertex path graphs. Each has N=3,K=1,M=2,R=3, with all vertices physically marked. Let a=3, beta=1/4, gamma=1/8. Each has Psi=15/8 and G=23/8. A whole-bank bridge produces

    N'=6, K'=4, M'=5, R'=6,
    Psi'=15/4=2 Psi,           G'=11/2.

Its marked spanning tree is the two paths plus the new bridge. The two bank derivative hits raise local k, leave the original physical marks intact, and introduce no extra original clock. A generic physical pair between nonprotected ports may lower R by two while preserving Psi, but consumes a finite local choice; pairing a protected port is refused by the type checker.

## 8. Exact escape operations and counterexamples

### 8.1 Heat is a neutral rewrite plus a currently unadmitted trace

In the unweighted ledger, for a variance-smoothing derivative a factor h^2 with h=alpha^gamma is accompanied by two derivative slots:

    a -> a+2gamma,       K -> K+2,       N unchanged.

Hence G, d and Psi do not increase. The same neutrality holds for the owned-clock charge q(k)=(k-2)_+ once k>=2; any low-order charge credit is finite and exhausted after finitely many insertions. Contracting those new slots is a SAME-VERTEX trace, not a bridge between distinct coefficient arguments. The loopless marked-spanning-tree certificate is lost. Iterating the formal heat operation gives arbitrarily high k at unchanged reserve; reserve truncation alone is not locally finite.

This demonstrates failure of the present certificate, not impossibility of a better heat theorem. Even the simplest generic cap warns against ignoring closed components: the matrix I has operator norm 1, but its full self-trace is D. A proper-cut estimate alone is not a dimension-free trace estimate.

The executable cold/warm VALUE shell avoids the literal loop. Its present improvement comes from a paid allocation/clock resource, however. With fixed N,K and no new resource the inequality

    N(beta_old-beta_new) > K(gamma_new-gamma_old)

spends finite positive beta budget. A sequence of positive gains tending to zero is not a well-founded proof. The original-target heat remainder remains a separate emitted type until a finite target-qualified source realization and bound are supplied.

### 8.2 A minimal zero-reserve old observer destroys grade local finiteness

Let B be standard scalar Gaussian, choose 0<epsilon<c and c+epsilon<=1, and use the perfectly admissible original gradient

    g(x)=c x+epsilon sin x,       g'(x)=c+epsilon cos x in (0,1].

Set W=B+alpha g(B). Here the old B itself has an unattenuated physical readout and is also the coefficient caller. Gaussian tilting gives

    log E exp(theta W)
      = theta^2/2
        + alpha [c theta^2+epsilon exp(-1/2) theta sin theta]
        + O(alpha^2).

The first-order term contains every even physical rank 2,4,6,... at the SAME alpha grade. Its rank-2m cumulant contribution from the sine is

    alpha epsilon exp(-1/2) (-1)^(m-1) (2m),       m>=1.

The retained observer contributed no positive reserve. This is an exact smooth, monotone-gradient counterexample to inferring rank local finiteness or strict residual gain merely from an O(alpha) perturbation. It does not disprove a redesigned representation; it identifies why the source-zero/ownership hypothesis is essential.

At the bare filtered-algebra level the same issue is R_o(x)=c_0 x with a degree-zero old label c_0. Then R_o(x)-R_o(y) has exactly the same degree as x-y. The formal iteration need not stabilize at any cutoff. If c_0=-1, F_o(x)=0 has no right inverse for a nonzero target. Even when 1+c_0 is invertible, solving that linear block requires a new leading map, not the claimed positive-gain proof.

### 8.3 Equality levels and nonlinear compilation cannot be erased

A bank-tree replacement has the correct whole-bank mean but generally leaves a random conditional coefficient. It emits another bank-dependent current. Likewise a private-source gradient does not certify the completed terminal sampler's curl, and an analytical all-cut tensor does not furnish an executable positive source. These are distinct types in the algebra, not interchangeable representations.

## 9. What group and field language can legitimately mean

After clearing rational grade denominators and imposing a genuinely finite positive filtration cutoff, the augmentation ideal I of the resulting finite jet algebra is nilpotent. Then

    (1+u)^(-1)=1-u+u^2-...,
    log(1+u)=u-u^2/2+...,
    exp(u)=1+u+u^2/2!+...

are finite expressions. The elements 1+I form a nilpotent filtered group (commutative for the scalar graph-product algebra). Tangent-to-identity formal substitutions give a generally noncommutative filtered composition group. The compiler inverse in Section 6 belongs to this triangular formal world when its hypotheses hold.

There is no natural finite field: coefficients are real/rational, positivity and factorial denominators matter, and a nonzero nilpotent element cannot lie in a field. A finite-dimensional jet quotient is still an infinite set over the reals. There is no useful single cyclic generator that manufactures all labeled shared-bank/ownership structures.

The algebraic analogy becomes useful only when written as explicit operations plus a terminating filtration. Without positive reserve, the finite quotient required by this paragraph may itself fail to be finite-dimensional, as the observer example demonstrates.

## 10. Formal closure, positive realizability and cost are three separate theorems

**Formal closure:** Enlarge the syntax to include loops, observer ports, heat marks, covariance/time branches and unresolved primitives. ExpandResidual is complete at a prescribed finite expression/jet cutoff. The admitted fixed-heated subgrammar has a reserve-qualified formal inverse under its stated contracts.

**Positive finite-VALUE realizability:** Each emitted type additionally needs an actual anchored VALUE source, all-cut/one-Hilbert and derivative-frame bounds, exact ownership, gapped positive covariance, same-endpoint native/current extraction, common-reference joining, untouched keep, actual complete first and terminal curl, and propagated finite floors. The local active C_k construction supplies one layer, not all these layers. The lead-current right inverse h is where these requirements enter the formal theorem.

**Cost:** Evaluate the actual versioned DAG. For each original force occurrence charge 2^(k_v+1) times its finite native/filter repetition count, original caller cost, exact clock multiplicity and every affected replay. Add scalar roots, readouts, captured tapes and precision work. A heat mark uses both cold and warm source branches. Exact caching requires identical complete source versions. Fixed-order native serial compilation may have zero inverse-alpha exponent while a whole-bank tree cubature or a resource-creating node split has a nonzero one. On the expanded DAG, sum parallel bills and multiply nested bills; inverse-alpha exponents follow the corresponding max-plus path recurrence.

Neither a finite number of operation schemes nor finite generation at each P implies c(P)=o(P). A repaired heat rule must return its resource balance and its actual cubature/query exponent together. The current framework gives an exact checklist and an explicit obstruction location, not that missing asymptotic theorem.

## 11. Operational next theorem to prove

The useful next object is not another list of ranks. It is a **certified residual dictionary** for the remaining escape generators. For each such generator, return:

1. its exact common-endpoint current and retained-bank equality level;
2. a leading-current right inverse using finite original VALUES;
3. every descendant type, including lost heat marks and exposed-observer branches;
4. a reserve increment bounded below by a target-qualified positive epsilon, or an alternative locally finite triangular block solver;
5. actual first/frame/curl/positivity and error-floor certificates;
6. a computable cost-and-replay bound.

If this dictionary is filled for all reachable escape types, Section 6 immediately supplies the repeatable formal architecture. If it is not, the generator returns the first unfilled type rather than relabeling a norm remainder as a source. That makes the algebra a practical way to discover and test the missing theorem.

## 12. A precise algebraic dualism: currents versus observables

There is an actual pairing, rather than only a metaphor. At a specified common endpoint X, let a finite current be J=(C_0,C_1,...,C_m), with all coefficient/tape ownership recorded, and let phi be a smooth compactly supported observable. Define

    <J,phi>_X = sum_(r=0)^m E[C_r : D^r phi(X)].

The r=0 term is C_0 phi(X); it keeps normalization or weighting branches visible until their cancellation is proved. Genuine probability-path currents annihilate the constant observable. That mass-zero fact does not permit deleting an arbitrary random C_0 without its exact relation.

Two current expressions are equivalent when this pairing agrees for every such phi. Let N_X be the annihilator of all observables; the observable current space is the raw current space modulo N_X. This is a linear quotient at a **fixed endpoint and retained-record contract**. Changing endpoint or exposing another tape changes the pairing and can invalidate a relation. Conditional equality requires the identity after multiplying by every bounded test of the retained record, not just after its unconditional average.

For example suppose X=X_0+eta Z, eta>0, and the keep Z is standard Gaussian independent of X_0 and every coefficient C_r. For r>=1 the exact Gaussian integration-by-parts relation is

    E[C_r:D^r phi(X)]
      = eta^(-(r-1)) E[(C_r:H_(r-1)(Z)) . grad phi(X)].

Here one coefficient output slot remains free. The associated rank-one vector is an analytic representative of the same current. Thus a high-rank current and a rank-one velocity can be dual-equivalent while having entirely different executability and cost. Conditional expectation of that vector given X gives the physical continuity velocity when integrable. The keep factor, Gaussian polynomial norm and coefficient one-energy bound are real costs, not discarded symbols.

If C_r depends on Z, product-rule derivative terms of C_r appear. The relation above is then invalid as stated. This is a precise algebraic form of the ownership requirement. It also explains why a same-endpoint current theorem is stronger than matching a log polynomial at different endpoints.

Let L be the leading-current map from typed executable source programs into this observable-current quotient. Then:

- Im L is the linearly realizable leading-current class under the specified source and ownership contracts.
- ker L comprises source combinations whose leading current vanishes; their finite native residuals need not vanish.
- J/Im L records obstructions to leading-current realization. One may call it an obstruction quotient; it is not automatically a cohomology group or a module over all graph operations.
- A right inverse h chooses a source representative for each reachable current class. Different choices differ by ker L and may have different feedback, resource use and query cost. The reserve filtration is on these certified typed representatives. A quotient identity that changes allocation or consumes a resource is not filtration-preserving until its full ledger is supplied.

The nonlinear positive compiler is not a linear duality isomorphism. Positivity belongs to the resulting probability law; the signed graph coefficients live in its tangent/current calculus. The meaningful problem is to enlarge Im L, understand its relations, and choose representatives whose residual map is filtration-improving.

### A conservation relation for apparently fresh analytic resources

An inserted positive resolvent R_l=(N+l)^(-1) does not create free new clock mass. Restoring the original target requires the full identity R_l(N+l)=Id. For the Gaussian Ornstein-Uhlenbeck number operator, N=x dot grad-Delta. Its compensator therefore emits both multiplication/derivative and Laplacian branches, including the same-vertex trace currently outside the marked-tree compiler. A dummy clock without this compensator changes the target; a clock unrelated to a private shield cannot claim that shield's node-mass inequality. This is an example where the algebra exposes the exact conservation law an attractive analytic rewrite must obey.

## 13. Three independent next construction tasks

### Task A: Build the current-image and residual table

Choose the smallest nontrivial family containing C0, C1, one connected bank bridge, the scalar rank-four example above, one covariance/time branch, and the proposed carrier-stable superposition. For each generator, compute its leading current in the fixed-endpoint observable quotient, an explicit VALUE-source representative, and its entire first nonmain residual. Record exact equality scope and retained tape names. Check all critical overlaps: derivative then join versus join then derivative; physical tilt before versus after connected extraction; grouped compilation versus isolated compilation. Differences belong either to a proved relation or to a named missing generator. Deliver a machine-readable table and a complete finite test, not a statement that all graphs are admissible.

Success condition: the source/current image and every escape at this small generating boundary are explicitly enumerated, so the filtered right-inverse theorem has verifiable inputs. This is independent of improving any heat bound.

### Task B: Resolve the minimal trace escape in the observable quotient

Start with one privately heated force occurrence carrying two new same-vertex derivative slots inside a marked tree. Retain the whole graph and caller/clock dictionary. Expand the heat derivative, all clock/resolvent compensators, and the keep integration-by-parts relation before bounding individual pieces. Ask whether the complete trace family is equivalent to a finite sum of already realized distinct-vertex/marked-difference currents, possibly plus a new primitive with an explicit source.

Success condition: either a concrete relation and finite VALUE right inverse with positive reserve/cost, or an explicit nonzero obstruction class under the current dictionary. A proof of nonmembership must refer to the defined image; failure of one norm estimate does not prove nonmembership. Cancellation of complete trace families is allowed even if isolated trace terms have poor bounds.

### Task C: Replace a zero-reserve tower by an invertible leading block

Use W=B+alpha[cB+epsilon sin B] as the required test. Its infinite same-grade observer tower is summarized by one exact Gaussian-tilt expression, theta E[g(B+theta)]. Treat that whole family as one typed functional generator instead of truncating physical rank. Determine whether its leading interaction with admitted currents has a computable block inverse at the retained-carrier boundary. If it does, redefine L to include that block and test whether the remaining residual gains positive reserve. If it does not, return the exact missing source or ownership relation.

Success condition: an explicit finite original-VALUE realization and complete first/ownership/current identity for the resummed block, or a precise obstruction. A symbolic inverse of a generating function alone is not success. The proposed carrier-stable native port can supply superposition inside its stated coefficient/carrier independence boundary, but must not silently absorb the counterexample's shared coefficient/carrier B. Its native-to-ideal floor is quadratic in the summed root envelopes, so it transfers a formal ideal jet only BELOW that floor unless those quadratic native coefficients are themselves identified. At higher cutoffs the floor remains an explicit unresolved row, not an exact all-orders native jet.

These tasks separate three useful advances: make the current dictionary complete, enlarge it across the heat trace, and change the leading algebra where reserve zero makes ordinary filtered inversion inapplicable. Each can progress without assuming the other two are solved.

## 14. A concrete trace-to-bridge relation

This is an additional exact algebraic construction, not merely a proposed research task. Freeze every scalar clock, width, contraction tensor and normalization. Give the N original force occurrences independent formal center variables z_1,...,z_N in R^D, and write F(z_1,...,z_N) for their contracted tensor product. Define constant-coefficient operators

    H_ind = (1/2) sum_v Delta_(z_v),
    H_com = (1/2) (sum_v D_(z_v)) dot (sum_v D_(z_v)),
    B_cross = sum_(u<v) D_(z_u) dot D_(z_v).

Then the exact operator relation is

    H_ind = H_com - B_cross.

Proof: expand the square in H_com. The diagonal terms are H_ind. Each unordered distinct pair occurs twice, canceling the factor one-half. Tensor contractions commute with these differentiations because their coefficients are frozen. If an actual coefficient is not frozen, its product-rule derivatives must be included as additional typed branches.

H_ind inserts the sum of the same-vertex heat traces. H_com differentiates the entire graph under one common translation z_v -> z_v+y. Each term of B_cross inserts a derivative slot at two DISTINCT original occurrences and pairs those slots. Thus its graph has one new edge, possibly parallel to an old edge, and never a new self-loop. The old marked spanning tree and every old physical mark survive. Repeating B_cross still only joins distinct original occurrences. This is coefficient-graph cut admissibility; native/current realizability is a separate field of the record.

For the smallest nontrivial example F=f(z_1)g(z_2) in one dimension,

    (f''g+fg'')/2 = H_com(fg)-f'g'.

The complete two-vertex trace family is therefore a common-translation derivative minus one loopless edge insertion. For N=1 there is no bridge, and the remaining common-translation generator is exactly the original local heat generator. Nothing has been canceled for free.

All three constant-coefficient operators commute. Accordingly their finite formal jets obey, for every m,

    exp(s H_ind)
      = exp(s H_com) sum_(j=0)^m (-s B_cross)^j/j!
                                                   modulo s^(m+1).

This is an identity of finite differential-operator jets. It is NOT permission to execute exp(-s B_cross) as a probability kernel; that factor need not be positive. It is NOT an assertion of convergence of a heat Taylor expansion down to zero smoothing for general C2 sources. At a positive shield its finite derivative coefficients exist; a useful finite remainder still needs its own actual heat/clock and source certificate.

There is also a useful single-operator interpretation. For a common Gaussian-bank derivative D and H=(1/2)Delta, define the failure of H to obey the product rule by

    {F,G}_H = H(FG)-(HF)G-F(HG).

The ordinary coordinate product rule gives {F,G}_H=DF dot DG, exactly the bank-bridge instruction. The bracket obeys the ordinary product rule in each argument. Thus diagonal heat and distinct-argument bridges arise from the SAME second-order operator; they are not unrelated combinatorial gadgets. No square-zero differential or extra cohomology structure is assumed.

### What this relation gains and does not gain

It shrinks the unresolved generator question: modulo a certified common-translation current relation, the sum of local heat traces is expressible through loopless bridge insertions. It offers a concrete first calculation for Task B: determine the image of H_com under the fixed-endpoint observable pairing and whether the induced boundary/covariance/observer terms already lie in Im L.

The relation does not by itself improve reserve. An s of grade 2gamma pays exactly the two extra inverse-width derivatives in the saturated charge regime. Nor does it establish a positive source for every added same-graph cyclic edge. If H_com is moved by integration by parts in a retained caller, all Gaussian-density, endpoint-chain and coefficient product-rule terms remain; the old observer counterexample explains why those terms can matter. If it is transferred to an untouched keep, its precise coefficient/endpoint dependence must first meet the keep relation's hypotheses.

A successful enlargement can therefore target a sharply specified pair of generators: the common-translation current and distinct-vertex same-graph edge insertion. The relation proves formal sufficiency for the full local-trace sum once both are admitted, while leaving their reserve, finite remainder, ownership and cost requirements visible.
