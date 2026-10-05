# Independent addendum: all proper cuts of marked-spanning multigraphs

2026-10-05. This addendum audits the separately published `MARKED-GRAPH-COROLLARY.md`. The previously audited source note and original audit report remain unchanged. The present result strengthens their restricted multigraph discussion.

## Verdict

**PASS.** A connected loopless multigraph with a specified spanning tree whose leaves retain external physical slots has the claimed product bound for every proper external tensor flattening. The earlier one-HS restriction is not the strongest available coefficient theorem: the topological-order extension supplies all proper cuts too.

This is a theorem for the fully contracted coefficient tensor. It does not establish the missing native marked-tree comparison port, nonlinear derivative-frame propagation, positive sampler return, observer-current closure, or sublinear cost. It does not authorize a local self-trace inside an arbitrary selected-spine chunk.

## Proof check

Fix a nontrivial split of external slots into physical inputs and outputs. Apply the already audited nonzero balanced-flow construction to the specified spanning tree. Every vertex has at least one incoming and one outgoing leg among its tree and external legs. The underlying tree has no cycle, so this orientation has a topological vertex order.

For each remaining graph edge, orient it from the earlier endpoint to the later endpoint in that same order. Looplessness is essential: the endpoints must be distinct. Parallel edges cause no problem. Every added edge respects the strict order, so the complete directed graph is still acyclic. Adding edges cannot erase any pre-existing incoming or outgoing leg. The resulting input/output split of every local tensor remains proper.

Initially the live wires are the physical input factors. Process vertices in the common topological order. All incoming edges at the current vertex have already been produced, and its external inputs were present initially. Permute those factors together, apply the local output-by-input matrix tensored with identities on the other live factors, and replace consumed wires by its outgoing wires. Every step has norm at most the local all-cut constant; tensor-factor permutations have norm one. After the final vertex, only physical output wires remain.

The matrix entries of this circuit sum precisely once over every internal edge index, so the circuit equals the original full tensor contraction. Consequently its operator norm is at most the product of local proper-cut bounds. No independence, positivity, planarity, or one-output restriction is used.

The one-HS statement also remains valid by the original marked-leaf attachment argument. Adding an external bank-derivative slot preserves the marked spanning tree. Joining two certified connected graphs by a bridge preserves a certificate: join their old spanning trees with that bridge; any leaf of the new tree was already a leaf of an old tree and remains marked. These facts do not themselves verify the native program or its derivative-frame comparison.

## Why the marking hypothesis is substantive

There is a simple loopless triangle outside this class. Let its vertices be A, B, C, with both external physical slots at A. For unit vectors u,v in R^D, take

    A_(a,b,i,j) = u_a v_b delta_(i,j) / sqrt(D),
    B = I_D,   C = I_D.

Every proper flattening of each local tensor has norm at most one. For A, cuts keeping i,j together have norm one; cuts separating them have norm 1/sqrt(D). Yet contracting the triangle yields

    sqrt(D) u tensor v,

whose physical operator norm is sqrt(D)>1. No spanning tree of this triangle has all its leaves marked: only A has external marks, while every nontrivial spanning tree has at least two leaves.

Thus looplessness alone does not suffice. The marked-leaf spanning-tree hypothesis supplies real content, and the corollary cannot be promoted to arbitrary cyclic networks or dimension-free traces. This is an abstract tensor counterexample, not a claim about a realizable original-gradient source outside the audited source assumptions.

## Independent executable checks

`check_independent_graph_corollary.py` imports no author graph generator. It constructs loopless multigraphs from recursively grown spanning trees, retains physical marks at every tree leaf, and allows parallel extra edges. Internal edge dimensions alternate between two and three.

For every tested physical input/output split it independently:

1. Constructs the exact integer balanced tree flow.
2. Topologically orders that tree and orients each extra edge forward.
3. Verifies every local cut remains proper.
4. Constructs the actual live-wire matrix circuit.
5. Verifies its operator norm after every vertex.
6. Compares the resulting matrix with a separate full `einsum` tensor contraction.
7. Checks the product bound for the completed cyclic physical flattening.

Final results: **3,667 assertions passed**, across 20 loopless multigraphs, 112 physical-cut circuits, and 484 live-wire steps. The maximum circuit/full-contraction discrepancy was 6.25e−17. The triangle counterexample was verified at D=2,3,7.

These numerical checks diagnose the implementation of the algebra; the general result rests on the acyclic live-wire proof above. The audited corollary's hash, diagnostic hashes, and unchanged-old-input verification are recorded separately in `GRAPH-COROLLARY-MANIFEST.json`.
