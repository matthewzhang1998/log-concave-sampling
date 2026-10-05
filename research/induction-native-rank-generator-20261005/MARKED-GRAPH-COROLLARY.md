# Full marked-spanning-graph cuts from an acyclic circuit

2026-10-05. This strengthens the coefficient-tensor scope in Section 8 of ACTIVE-PROBE-GENERATOR.md without changing those independently audited bytes. The same result was independently derived in the shared-variable work's MARKED-SPANNING-TREE-EXTENSION.md. It is not by itself a native sampler-return theorem.

## Theorem

Let G be a finite connected loopless tensor multigraph with a specified spanning tree T. Every graph-theoretic leaf of T must retain at least one external physical slot. Suppose every proper flattening of each local tensor A_v has norm <=c_v. Then **every proper external flattening of the fully contracted graph** has norm <=product_v c_v.

The one-Hilbert bound is also preserved: root T at a marked leaf whose tensor has HS norm h_root; the completed graph has HS norm <=h_root product_(v!=root)c_v. For the active original-gradient response tensors, this is one sqrt(D) with explicit local rank/width constants.

## Proof of every proper cut

Fix a nonempty physical input set and its nonempty output complement. Orient the spanning tree using a balanced boundary flow, exactly as in Section 8 of the source note. Its explicit masses w_j=2^(2^j), multiplied by the opposite-side sum, ensure every tree edge has nonzero flow. Every vertex already has at least one incoming and one outgoing tree or boundary leg.

Choose a topological vertex order of that oriented tree. Direct every remaining edge from its earlier endpoint to its later endpoint. Parallel edges are allowed. Since no edge has equal endpoints, all additions respect the same strict order. The enlarged graph is a DAG. Adding these edges does not remove the incoming and outgoing legs already established at a vertex. Its local tensor is therefore still used only through a proper flattening.

Process vertices in topological order. Regard the still-unconsumed physical inputs and internal edges as live tensor-factor wires. At a vertex, permute its input wires together, apply its local operator tensored with identities on other live wires, then carry its output wires forward. Every step has norm at most c_v, and wire permutations have norm one. After the last vertex the live wires are exactly the physical outputs. This circuit is algebraically the original tensor contraction, so its norm is at most product_v c_v.

The one-Hilbert proof attaches vertices in a tree-parent-before-child order, contracting every edge to the built group at once. A new vertex retains a future tree child or a physical mark, keeping its local attachment cut proper. Hence only the initial root contributes an HS norm.

## Source and induction consequence

At each force occurrence, the active original-gradient source supplies the normalized local tensor and the numerical all-cut constants from equations (2)-(4), or their explicitly paid bounded-preparation floors (8). A coefficient-bank derivative may add an external slot at a vertex without destroying the marked spanning tree. A bridge between separate certified graphs preserves a certificate by joining their old spanning trees with that bridge. Therefore C1 hits on formerly C0 marked leaves are not an intrinsic tensor-cut obstruction.

This theorem controls the **fully grouped analytical coefficient graph**. A native execution still needs an actual finite source/pair/filter program, positive covariance gaps, complete first paths, finite same-endpoint return census, extended derivative frames, retained-observer ownership, numerical Sobolev/curl floors and replay. In particular this theorem alone does not authorize estimating a self-trace inside an arbitrary selected-spine chunk. The global circuit proof and the native frame comparison are distinct certificates.

Self-loops are excluded. Without a marked-leaf spanning tree, neither the orientation argument nor the resulting product-cut claim is generally valid. A component with no physical terminals cannot simply be erased or assigned a dimension-free trace.

## Validation

The independent audit addendum supplies direct full-network-versus-live-wire-circuit equality tests with heterogeneous edge dimensions, proper-cut inequalities, and a counterexample outside the marked-spanning-tree class. This is separate from the original source and tree-only audit and does not modify their pinned inputs.
