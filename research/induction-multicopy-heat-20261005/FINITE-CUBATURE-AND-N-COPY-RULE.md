# A computable n-copy rule and a conservative positive finite cubature

2026-10-05. This closes finite executability of the coefficient census, with an honest potentially expensive numerical bill. It is not a public-log cubature theorem, an unqualified native smallness theorem, or an eventual-sublinear high-accuracy claim.

## 1. Exact finite input and census

Input consists of n>=2 finite original-VALUE coefficient graphs, their marked spanning trees, original clock labels and positive weights, affine scalar query rows on the COMPLETE common standard bank B, positive original private widths, physical marks/readouts, source versions and numerical floors. All retained observer and source-zero physical bank rows must meet the stated return boundary; otherwise this is only a tensor identity.

For each labelled replica tree T, enumerate one original force hit at each end of every bridge. An additional q caller derivative is a finite ordered force-hit list with its actual injection rows. Let N_i be the number of forces in argument i and N=sum N_i. For fixed input histories the exact number of tree/bridge-hit combinations is

    H_n=(product_i N_i) N^(n-2).                         (1.1)

Indeed one tree contributes product_i N_i^(degree_i), and weighted Prüfer gives the displayed sum. This is 243 for 3/3/3 and 648 for 3/3/6. Each of q additional ordered caller hits has at most N choices, giving at most H_n N^q histories before merging. Signs, factorials, physical permutations and any identical-input multiplicities are inserted once under a fixed convention.

Each history has

    K_new=sum_i K_i+2(n-1),
    k_v=k_v,old+number of bridge hits at v.

The union of original spanning trees and the replica tree bridges is a marked spanning tree. All additional Price edges used below connect distinct force occurrences, so they create no self-loop and preserve that certificate.

## 2. Exact redistribution algorithm

Compute the known scalar K_T(t), c=min t, delta_i=1-max incident t, and H=K-c11^T-diag(delta)>=0. Retain the actual B row sqrt(c) and add the independent H-root and delta-copy innovations.

Within each delta-copy innovation, inspect the recorded independent scalar query blocks. From any known query covariance Gamma choose a diagonal D>=0 with Gamma-D>=0 and move delta D into independent force shields. A known small scalar eigenvalue is a conservative choice; the simultaneous cubic maximum floor m_star=max_j m_j is a sharper choice. One executed choice must support all its claimed caller derivatives; one may not change its heat allocation merely to bound a different caller. For a previously constructed multi-copy packet, its recorded independent residual blocks remain available INSIDE the new delta innovation: extracting heat there preserves the packet's exact original conditional coefficient. Never alter an already recorded old disintegration and assume its joint mixed currents stay unchanged.

For each node/history compute its literal scalar factor b_h from original weights, new positive quadrature weights, actual inverse widths, injection scalars, coefficient signs, fixed source-response factors and actual readout products. Compute separately:

    S0=sum |b_h|,
    S1=sum |b_h| sum_v a_v/t_v,
    S2=sum |b_h|^2,
    S2,1=sum |b_h|^2 (sum_v a_v/t_v)^2,

and their nodewise maxima. Here a_v bounds the ACTUAL complete captured-caller row, not a guessed unit row. All of these are finite scalar sums, computable without an original derivative or tensor oracle. Repeated clock labels use their actual combined powers. Exact dyadic-cell interval bounds can certify the sums without enumerating each tensor entry.

This is the general n-copy heat/importance allocation rule: allocate available independent covariance to the forces with the largest unpaid inverse-width powers, compare the resulting scalar envelopes, and use the smallest certified envelope. The selection is deterministic from known clocks/history. Freeze ONE executable allocation for each main history and evaluate S0, S1, S2 and S2,1 at that SAME allocation. Choose among candidates against the entire required derivative budget; never reselect D separately for different caller proofs. The simultaneous m_star construction does this automatically for cubic copies. This is not importance sampling a random tensor or changing the target measure.

A safe fallback uses no extra extraction at all. If all old private widths are at least sqrt(eta), and actual row norms are at most A_row, the q-hit coefficient envelope is at most

    H_n N^q A_row^(2(n-1)+q)
       (product of literal input masses and fixed source constants)
       eta^(-(K_new+q)/2).                              (2.1)

Use the actual local constant at each k_v+hit, rather than treating fixed-rank factorials as one. This rank-explicit fallback proves finite bounds; it need not give useful high-accuracy scaling. The three-copy notes prove substantially better exact integrated envelopes in their stated families.

## 3. Finite positive cubature without a tensor oracle

Fix one tree and bridge-hit history, with the original finite clocks frozen. Let F(t) be its tensor integrand after Gaussian expectation. Evaluate no tensor entries in execution. F is only an analytical object used to certify a scalar quadrature rule.

For every normalized local original-gradient tensor of order k choose a known constant C_k bounding all proper cuts and HS/sqrt(D). The active C_k source theorem supplies explicit constants; a conservative larger factorial constant is allowed. Let

    B0=|w| product_v C_(k_v) sigma_v^(-k_v)

include all literal known old coefficients and injection products. For an edge e, order-sector differentiation of K_T gives derivatives of its pair entries in {0,1}. Gaussian Price differentiation therefore sums over at most n(n-1)/2 distinct replica pairs. For each such pair (i,j), expand its N_i N_j force hits. Compute B1,e by summing the same displayed product with k_v incremented once at each selected force, multiplying the two new actual row factors. This is a finite KNOWN scalar expression.

Each derivative adds one edge between distinct replicas to the already marked graph. The marked-spanning-tree theorem yields

    ||partial_e F||_HS <= B1,e sqrt(D),
    every proper cut <=B1,e,                            (3.1)

where derivatives exist on each open order sector. Original private heat is positive; Gaussian covariance differentiation is justified even when a bank covariance is singular by positive regularization and its limit. F is continuous across sector walls. Integrating piecewise along a line gives the global bound

    ||F(t)-F(t')||_HS <=sqrt(D) sum_e B1,e |t_e-t'_e|.     (3.2)

No hidden dimension-sized density derivative or numerical tensor entry is used. The extra Price edges preserve the existing spanning tree, so (3.1) still spends only one Hilbert factor.

Partition each edge-gap interval [0,1] dyadically down to the original eta, with one last interval [0,eta]. Further divide each panel into equal midpoint subintervals of length at most h. Tensor the midpoint rules over the n-1 edges. All nodes are interior, all weights are positive, and they sum to one. Each node weight u<=2s at its gap s; each original dyadic panel retains its exact total mass. Thus the bridge mass envelopes used in the three-copy proofs remain valid.

By (3.2), a sufficient mesh condition for absolute coefficient error <=epsilon sqrt(D) is

    h <= min(1, 2 epsilon/max(1,sum_e B1,e)).             (3.3)

For all histories together, replace sum_e B1,e by the literal sum over histories, or allocate epsilon_h with sum epsilon_h<=epsilon. The number of nodes on one edge is at most

    2+ceil(1/h)+ceil(log2(1/eta)),

and the tensor rule has at most that quantity to the power n-1. This is an explicit finite certificate. The possibly large inverse-width and inverse-tolerance dependence is part of the full work bill. A higher-order positive rule may improve it, but no improvement is assumed here.

The same rule can be used for all derivative histories by taking the maximum required known Lipschitz bound. The conservative majorants can be computed at original private widths, even if the executed source uses the improved exact redistributed widths. Thus there is no circular dependence of the quadrature proof on sampled source outputs.

## 4. Native positive-law admission and cost are separate checks

At each node compile the marked force graph from complete original-gradient VALUE sources. All numerical covariance roots, original private shields, source-origin anchors, pair/filter/calibration tapes, active probe rows, physical readouts and untouched keep are retained exactly as their consumer requires. Known scalar signs live in real root amplitudes. There is no tensor coefficient, high original derivative or Monte Carlo expectation oracle.

For a graph with N forces and target alpha^a b_h, put nonroot amplitudes alpha^beta (or a declared weight-aware version) and root amplitude b_h alpha^(a-beta(N-1)). Before admission check the ACTUAL complete radius sums, caller paths, covariance gaps, readout shares, finite native priors and numerical floor propagation. The width-zero fixed-graph conditional theorem contributes C_h rho_root,h^2 sqrt(D), but it does not erase b_h, S1, any old endpoint first, or the C_h dependence on chosen readout/gap geometry.

A conservative midpoint count can increase readout and source sums by inverse-alpha powers when epsilon or old cutoffs scale with alpha. Those powers must be entered into the grade inequalities. This paper supplies a finite original-VALUE coefficient construction and an executable sufficient admission test, not a proof that the test succeeds at every intended P. The proposed public-log/high-order quadrature route is an outstanding improvement, rather than a silently substituted theorem.

For actual node count Q_T, the force work is

    sum_(T,h,node) sum_v 2^(k_v+1)
       M_native(k_v,b_v,epsilon_v)
    + complete original captures + old complete requested versions
    + known scalar replica/residual roots + physical rows and keep
    + all affected ancestors/anchors/replays.

Add live and retained D-root storage, scalar arithmetic precision, finite means and calibration, and every source Sobolev floor. Rank constants, graph count, node count, wrapper order, inverse-width VALUE precision and old complete replays remain visible. No promise of c(P)=o(P) follows from this finite certificate.

## 5. A stronger integrated rule for any number of old cubic copies

N-CUBIC-RANK-EXPLICIT-BOUNDS.md proves an explicit uniform rank formula, beyond the crude private-width fallback. For n old cubic copies and q additional analytical bank/caller hits, put

    e(n,q)=max(n+q-4,0).

With the same executed m_star allocation, the complete integrated all-cut and one-Hilbert allowance is

    C(n,q) L^(6n-1) alpha^(3n) eta^(-e(n,q)),

where eta bounds the ORIGINAL private and structural endpoint variances. The leading node factors and their q-hit envelopes have nodewise bound C(n,q) eta^(-e(n,q)); literal squares have integrated bound C(n,q)L^(6n-1) eta^(-2e(n,q)). The q=0/1 cases give the actual coefficient and captured-caller budgets. Higher q describes analytical heated-tensor derivatives, not a promise to execute higher derivatives of g.

This explains the special three-copy closure, and adds a fourth-copy main coefficient with only logarithmic costs. From five copies onward this safe rule exposes genuine positive inverse-cutoff powers in its upper bound. It is not a proof that those powers are necessary.
