# Scalar positive lattice-clock setup: finite cost and conditioning

2026-10-05. This addendum makes the scalar construction in Section 3 of the main note explicit. It concerns known public coefficients, not an oracle for g.

A panel has n consecutive integer sites, positive weights exp(-2 Delta j), and span at most delta<=log 2. The ratio of its largest and smallest atom weights is at most 4. Rescale its support to x_r=r/(n-1), r=0,...,n-1; the n=1 case is exact. Its normalized moment generating function is

    G(z) = [(1-q^n exp(n z/(n-1)))/(1-q exp(z/(n-1)))]
           /[(1-q^n)/(1-q)],     q=exp(-2 Delta).

Compute its Taylor jet through order 2m by finite exponential-series multiplication and division, then multiply coefficient ell by ell! to get normalized moments mu_ell. This uses O(m^2) scalar arithmetic operations and evaluation of the known exponential coefficients, irrespective of n. Large exponents n are handled by multiplication of public logarithms/exponents; they are not loops over n sites. Tiny differences use expm1/series with certified error.

Let H=(mu_(i+j))_(0<=i,j<m) and H1=(mu_(i+j+1)). For n>=m, H is positive definite. Its Cholesky factor L defines the symmetric Jacobi representation J=L^(-1) H1 L^(-T). Compute its m real eigenvalues x_i and weights w_i=(v_i)_0^2 for orthonormal eigenvectors v_i. Then 0<=x_i<=1, w_i>0, sum w_i=1, and the rule is exact through degree 2m-1. Multiplying by the explicit panel mass and undoing the scalar rescaling gives the desired public positive rule. Cases n<m retain the n atoms exactly.

A crude conditioning bound suffices for the claimed polylogarithmic cost. For n>=m choose m grid sites separated by at least 1/(2m). Their Vandermonde matrix V on monomials of degrees 0,...,m-1 satisfies

    |det V| >= (2m)^(-m(m-1)/2),
    ||V||op <= m,
    lambda_min(H) >= (1/(4n)) sigma_min(V)^2
                  >= exp[-C m^2 log(2m)-log(4n)].

Thus certified Cholesky/moment arithmetic requires polynomially many bits in m, log n, log(1/Delta), and the requested absolute precision. The resulting irreducible Jacobi off-diagonal entries have lower bounds of the same exponential-polynomial form. Standard finite symmetric eigenvalue isolation or characteristic-polynomial interval isolation then certifies positive weights and separated roots using polynomially many bits and operations at this fixed finite degree. Equivalently one can use exact algebraic moment data with interval refinements for their exponential coefficients. No inverse-n precision or operation COUNT appears: log n bits do.

At finite numerical precision, store positive rational weights certified within the allowed intervals and normalize them by their positive total to the encoded exact panel mass. Store nodes certified inside the panel. Moment exactness is analytical; the total encoded operator discrepancy gets its own absolute floor. The sum of absolute weights is the known finite mass. Apply the preallocated downstream path sensitivity budget to all moment, node, weight, bridge-row and response-width errors.

The accompanying checker implements the Taylor-jet moments and Hankel/Jacobi rule with arbitrary-precision mpmath. It also runs a panel containing over 2^19 sites without enumeration; this tests the stated algorithm rather than substituting an atom-by-atom setup.
