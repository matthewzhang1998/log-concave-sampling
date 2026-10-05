# Marked spanning trees: the broader cut lemma is valid

2026-10-05. Independent follow-up to the sealed decorated audit. The earlier audit directory remains unchanged. This addendum corrects the earlier suggestion that an interchunk edge at a marked C1 leaf is itself a dimension-safe graph/topology obstruction.

## Verdict

The proposed marked-spanning-tree cut lemma is true under its literal hypotheses. Its argument by balanced terminal flows and an acyclic orientation is sound. The same certificate survives connected conditional Wick offspring and bank-hit bridges joining separate copies, including bridges incident to previously protected C0 leaves after their promotion to marked C1 vertices.

Consequently a marked-leaf bank bridge is not, by itself, a missing native graph topology or all-proper-cut resource. It falls outside the narrow center-plus-C0-leaf invariant, but admits the broader marked-spanning-tree construction below. Uniform source/clock/frame bounds, actual observer ownership and a well-founded mixed correction queue still require their own proof. No full arbitrary-order sampler follows from this cut lemma alone.

## 1. Statement

Let G be a finite connected loopless force multigraph. At vertex v let T_v be a tensor with one finite-dimensional Hilbert-space slot for every incident edge and every physical external port. Pair internal edge slots by their Hilbert inner product. Assume that for every nontrivial bipartition of its slots, the corresponding local operator norm is at most M_v.

Suppose G has a spanning tree T such that every graph-theoretic leaf of T has at least one physical external port. For every nontrivial partition of all external ports into global input and output ports, the contracted network has operator norm at most product_v M_v.

For a network with at least two external D-dimensional physical ports, isolating one physical port gives full Hilbert norm at most sqrt(D) product_v M_v. Different port dimensions use the square root of the chosen isolated port dimension. The one-vertex/vector exceptions require the corresponding local energy statement and are not hidden in the term “proper cut.”

The result does not need tensors to be symmetric, or G to be simple. It does require ordinary pairwise edge contractions, local control of every proper cut, and the stated marked-tree certificate.

## 2. Every tree edge separates terminal-containing components

Delete an edge of T. Each of its two components contains a leaf of the original T: start at the cut endpoint and follow a path away from the deleted edge until it terminates; that terminal vertex is an original tree leaf unless the component is a singleton, in which case its unique vertex was itself an original leaf. Therefore both components contain an external port.

This fact is essential to the generic-flow argument. An unmarked branch attached to a marked core does not meet it.

## 3. Balanced positive terminal flow

Fix any nonempty proper global input/output partition. Put a positive mass at each input port and a positive demand at each output port, with total supply equal to total demand. The positive balanced assignments have nonempty relative interior.

For each tree edge e, choose one component A of T-e. Its net supply is a linear functional on this balanced assignment space. It is not identically zero: A contains a nonempty proper subset of all terminals. A net-supply functional proportional to the global balance functional would require either no terminal in A or every terminal in A, both excluded by Section 2. Consequently the zero-flow assignments for e lie in a proper hyperplane.

There are finitely many tree edges. Choose a positive balanced assignment outside the union of those hyperplanes. Each tree edge then carries its unique nonzero flow from the positive-net component toward the negative-net component. Treat every global input port as an incoming edge supplying its mass, and every output port as an outgoing edge demanding its mass.

Conservation at each force vertex implies that its tree-plus-terminal ports contain at least one incoming and one outgoing port. Indeed every tree edge has nonzero magnitude, all external magnitudes are positive, and an all-incoming or all-outgoing vertex would violate conservation. This includes vertices having no external ports and vertices carrying several external ports assigned to opposite global sides.

## 4. Add every cycle edge without producing a directed cycle

Any orientation of the tree is acyclic. Take a topological order of its vertices. Orient every non-tree edge from its earlier endpoint to its later endpoint. Parallel edges are treated separately. G is loopless, so their endpoints differ and each can be oriented strictly forward.

The entire directed graph is now acyclic. Every vertex still has its incoming and outgoing tree/terminal ports from Section 3. Adding more directed ports cannot make its local cut improper.

View T_v as an operator from the tensor product of its incoming slots to the tensor product of its outgoing slots. Its norm is at most M_v by the local all-proper-cut hypothesis. Process vertices in topological order. Each step is this local operator tensored with the identity on all other live wires, followed by a permutation of tensor factors. All wires are ordinary Hilbert-space identity maps. No cap, cup, diagonal copying map or normalized/un-normalized trace is introduced. Thus each step has norm at most M_v, and composition has norm at most their product.

At the end the live wires are precisely the global output ports; initially they were precisely the input ports. The circuit contraction is exactly the original tensor-network contraction. This proves the global cut statement. Isolating one D-dimensional physical index then gives the single-sqrt(D) Hilbert estimate.

For complete index bookkeeping, let the live-wire set immediately before a vertex step contain: every unconsumed global input, each internal edge whose tail has been processed but whose head has not, and every already produced global output. Each internal edge has one tail and one head, so its index is produced exactly once and consumed exactly once. Each external input is consumed once; each external output is produced once and carried to the final boundary. At a vertex, all its incoming wires are already live and its outgoing wires are fresh, because every internal edge points strictly forward. Apply the local matrix to precisely these incoming tensor factors and leave all other live factors unchanged. Distinct output slots are distinct factors, even if their dimensions agree. This never copies an index, inserts a diagonal map, or contracts two live outgoing wires. In the component formula, multiplying the matrices sums exactly once over each internal edge index, which is precisely the defining network contraction. Every vertex with no incoming INTERNAL edge has an incoming external leg, and every vertex with no outgoing INTERNAL edge has an outgoing external leg. Thus there is no hidden vector-source or scalar-sink factor lacking the required physical boundary.

This proof is stronger and cleaner than attaching chunks and attempting to prove global cut inheritance by an informal repeated Cauchy–Schwarz argument. It explicitly eliminates the feared self-trace by constructing an acyclic orientation for every desired global cut.

## 5. Native selected-spine realization

The same marked spanning tree supplies a legal selected-spine decomposition when each local original-VALUE C_k adapter admits two selected ports and k remaining external ports, matching the literal vertex valence k+2.

Choose a root at a vertex with a physical mark and use one such mark as the root physical endpoint. Direct the tree away from this vertex. At each unmarked nonroot vertex, select one child tree edge as continuation of the selected spine. Such a child always exists: an unmarked nonroot cannot be a tree leaf. At a marked nonroot vertex, terminate its incoming selected spine at one designated physical mark. Its child tree edges initiate side spines through the other external adapter ports. At the root, select one child edge as the other selected port; its remaining children likewise initiate side spines. Every spine eventually terminates at a designated physical mark.

Each local tensor uses exactly two selected ports. All unused tree-child edges are side-output probes; all non-tree edges are cut and supplied by fresh shared Gaussian colors; every remaining physical port is a direct-public probe. Valence k+2 gives exactly k such external probes. A marked internal C1 vertex with two tree incidences and one physical mark is legal: use its parent edge and physical mark as selected ports, and its other tree edge as the single side probe. Nothing requires it to remain a two-force center-C0-leaf chunk.

This establishes the static source grammar relative to the admitted fixed-order C_k/selected-pair compiler. It does not by itself prove new coefficient normalizations, width/clock guards, costs or same-endpoint remainder estimates for an unbounded collection of adapter orders.

## 6. Conditional offspring preserve a marked spanning tree

A partially inserted reverse root-spine monomial retains complete selected spines. An omitted side contributes a centered boundary innovation at its parent; it does not delete a selected edge inside an already retained spine. Every leaf of the retained force tree is therefore a designated physical endpoint of one selected spine. The monomial has its own marked-leaf spanning tree before Gaussian averaging.

In a product of monomials, take the disjoint union of these trees, keeping repeated force labels as distinct occurrences. Wick contractions add edges; they do not remove existing selected edges or consume the designated selected physical endpoints. Direct-public marks may be consumed, but they are not needed for this certificate. Under the usual once-per-local-color condition, the added edges are loopless.

A connected cumulant diagram connects these monomial trees. Contract each monomial tree to a temporary supervertex and choose a spanning tree of the connected supergraph, using one actual Wick edge for each chosen superedge. The union of the original monomial trees and these connecting edges is a spanning tree of the whole force graph. Adding the connecting edges can only increase old degrees, so every leaf of the new spanning tree was already a leaf of one monomial tree and remains physically marked. Unused Wick edges are merely non-tree edges.

This argument does not require each center to retain its own C0 leaf, does not require a short root spine, and does not require cycles to stay within a center-only graph.

## 7. Bank-hit bridges and the exact correction to the earlier boundary

Differentiate a tensor at a physically marked vertex: a C0 leaf becomes a C1 tensor with the same physical mark and an additional bank slot. Its original tree edge and physical slot remain. If two differentiated tensor copies are joined through their bank slots, that joining edge connects two otherwise separate marked-tree copies. The union of both trees and this bridge is again a spanning tree; none of its leaves has lost its physical mark. The all-proper-cut theorem applies to the resulting local higher-order tensors, assuming their actual source estimates and bank-injection rows are supplied.

Thus the prior statement that such a marked-leaf bridge lacks a dimension-safe native topology is too strong. What failed was the narrower center-spanning-tree/two-force-side-spine proof, not the existence of a legal broader proof. This addendum supplies the broader topology and cut invariant.

Important boundaries remain:

- A self-loop contraction at one force occurrence is outside the loopless theorem. Every proposed whole-bank history must verify its occurrence-level edge structure; calling two derivative slots different labels is not sufficient if they belong to the same tensor occurrence.
- Retaining a physical root as an outside observer is an ownership question, not a tensor-cut question. The retained-public obstruction is unchanged.
- Higher C_k orders need their actual common-heat and one-hit/multiple-hit clock envelopes. Static cut inheritance does not establish those bounds.
- A mixed old/new current can contain only one new root. The pure conditional recurrence d'=h*d with h>=2 is not automatically a well-founded measure for that mixed queue.
- Complete dynamic-bank Jacobian/frame estimates and the positive same-endpoint consumer must be extended to the chosen finite larger graph family; a static graph norm alone is not a transport theorem.

The correct next question is therefore quantitative/recursive closure of this larger grammar, rather than existence of a native topology for marked C1 leaf bridges.

## 8. Why the marked-tree hypothesis is substantive

Without it, local all-proper-cut bounds do not bound every global proper cut. Here is a fully symmetric tensor example.

Take dimension D>1 and unit vector u=e_1. Form a triangle with tensors A,I,I at its vertices, and attach A by one extra bridge to a fourth tensor B carrying two external ports. Put

    A_ijk = (u_i delta_jk + u_j delta_ik + u_k delta_ij)/sqrt(D+8),
    B_ijk = u_i u_j u_k.

All local proper-cut norms of A are at most one: its 1-versus-2 flattening Gram has eigenvalues 1 in direction u and 2/(D+8) on u-perpendicular. Symmetry covers every cut. B and the two identity tensors also have every local proper-cut norm one.

Closing the triangle traces two slots of A and produces

    ((D+2)/sqrt(D+8)) u.

After the bridge contraction with B the global input-to-output operator is ((D+2)/sqrt(D+8)) u u*, whose norm exceeds one for D>1 and grows like sqrt(D). The graph is connected and loopless, but every spanning tree has an unmarked leaf inside the triangle. It therefore violates exactly the marked-spanning-tree hypothesis.

The proposed certificate avoids this genuine self-trace obstruction; it does not prove a false all-graphs theorem.
