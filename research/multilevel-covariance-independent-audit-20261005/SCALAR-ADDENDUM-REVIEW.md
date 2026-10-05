# Review of the separately supplied scalar setup addendum

2026-10-05.

**PASS.** The inspected author addendum is `SCALAR-CLOCK-IMPLEMENTATION-ADDENDUM.md`, SHA-256 `cc4404fe4becbb1c5f4768e043ca4f1a6d0d531bfc10f15d3b5a1f7afd06f89e`; its immutable copy here is `SCALAR-ADDENDUM-SNAPSHOT.md`. The main draft remains unchanged at SHA-256 `9bc692e910de080b46690bd698e460e386ce6ee1a4de6c6c7d39cddb293ff40e`.

The normalized geometric moment generating function is correct, including the index scaling by n-1. Its Taylor coefficients multiplied by factorial give the moments. The O(m²) jet computation does not enumerate n atoms.

For normalized moments, the Cholesky construction J=L^(-1)H1 L^(-T) is multiplication by x compressed to the orthonormal-polynomial span. The first polynomial is the constant one. Thus the eigenvalues are Gauss nodes, and the squared first eigenvector coordinates are their weights. The rule is positive and exact to degree 2m-1. If n<=m, keeping the original atoms avoids unnecessary endpoint or rank degeneracy; the supplied code does this.

The addendum's spaced-subset Vandermonde lower bound is valid and improves the more conservative bound in the main independent audit. For n>=m, m lattice points can be chosen with separation at least 1/(2m), while their normalized weights are at least 1/(4n). The resulting log inverse condition bound is O(m² log(2m)+log n). Certified moment, Cholesky, eigensystem, and positive weight computations consequently need only polynomial bit precision. The explicit root/weight lower-bound argument in Section 3 of the independent audit can be used with this stronger Hankel estimate.

The finite-precision instructions correctly distinguish analytical moment exactness from encoded coefficient floors. Positive weights can be encoded and renormalized; node encodings remain inside their panels. Absolute numerical budgets retain all downstream factors, including native normalization, response widths, and small reserve h. For A=0 the exact zero source is returned directly rather than normalizing by zero; for positive A the inverse-radius arithmetic affects precision and remains explicitly budgeted.

The author checker, SHA-256 `957f9221dcce16dc66d3322414c22a61d198ef9cad5ddfcf39708af14793ad57`, was copied to `author-check-rerun/` and rerun there without modifying any author files. All 2,004 assertions passed. It constructs an eight-node positive Gaussian rule for a panel with 524,288 atoms using 200-digit moment jets, without enumeration. The rerun results match the supplied result hash `b3df42c185445c9d677a5a21f7fc5e4608501de49b2c12ff8362cb030e6c8622`.

This is a working arbitrary-precision numerical prototype plus a constructive precision proof, not a production interval-certified solver. No additional mathematical blocker to the scalar public-log setup was found. The original audit's separate mean-LAW/caller/final-compiler limits remain unchanged.
