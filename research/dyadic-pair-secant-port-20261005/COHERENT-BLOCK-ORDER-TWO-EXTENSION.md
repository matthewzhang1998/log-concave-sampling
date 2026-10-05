# Order-two interaction truncation for the full coherent history block

2026-10-05. Analytical extension of FINITE-ORDER-TWO-C2-TAIL.md. This uses the actual pair (F2,F3-F2), with a single coherent future. It supplies no conditional coefficient oracle or executed history grid.

## Theorem

Under the source, physical block and dyadic hypotheses of the F2 tail theorem, let

    W=(F2,Delta3),  Delta3=F3-F2,
    H_k^W=E[W|S_(k+1)],
    r_k^W=P_(k,>2)H_k^W,
    Q_k^W(Y)=E[r_k^W(r_k^W)*|Y].

The physical output is R^(2D), with its two slots coherent. There is an absolute constant C such that

    ||Q_k^W||_(L2(Y;HS))
      <=C sqrt(D) min{A^2 d_k^2,
                      A^5 b^(3/2)(1+d_k^(-1/2))}.                 (1)

If A^3 b^(3/2)<=1, then

    ||sum_k Q_k^W||_(L2(Y;HS))
      <=C sqrt(D)[A^(22/5)b^(6/5)
                   +A^5 b^(3/2)(1+log_+T)].                     (2)

As always Q_0^W=Q_1^W=0. Every Q_k^W is PSD and includes its F2/Delta3 cross block. Retaining the base covariance and every singleton/pair block at all retained levels therefore gives a positive approximation to the COMPLETE block covariance, with the paid error (2) plus the terminal dyadic/tail residual. For fixed b and logarithmic horizon, this is o(A^4 sqrt(D)). No such dimension-uniform conclusion is asserted for b=D growing arbitrarily with A.

## Proof

Work first in one physical block of dimension n. The true history difference obeys, by the original A-Lipschitz recursion and Minkowski,

    ||Delta3||_Lp <= A^3 ||X_0||_Lp,  p>=1.                        (3)

Indeed F1 has Lp energy at most A||X0||p, Delta2=F2-F1 at most A||F1||p, and Delta3 at most A||Delta2||p. These are energy estimates only; the full private/caller first of Delta3 is still the O(A) bound supplied by the complete sensitivity kernel.

The F2 proof constructs an actual comparison V_k satisfying

    P_(k,>2)H_k^(F2)=P_(k,>2)V_k,
    ||V_k||_L4^2 <=C A^5 n^2 J_k,
    J_k<=C(1+d_k^(-1/2)),  J_k>=1-e^-T.

Use the stacked comparison

    V_k^W=(V_k,E[Delta3|S_(k+1)]).

Its high-order projection is exactly r_k^W. Conditional Jensen, (3), the Gaussian eighth/fourth moments and T>=1 give

    ||V_k^W||_L4^2
       <=C[A^5 n^2 J_k+A^6 n]
       <=C A^5 n^2 J_k.                                      (4)

There is no assertion that P_(k,>2) contracts L4. For every fixed test vector in R^(2n), conditional L2 orthogonality instead gives

    E[r_k^W(r_k^W)*|S_k]
           <=E[V_k^W(V_k^W)*|S_k].                           (5)

Taking conditional expectations onto Y, using HS<=trace for the positive right-hand side, then conditional Jensen, proves the one-block coarse bound in (1), C A^5 n^2 J_k. Thus the true mark and every high-order cross term are bounded together, without separately multiplying dimension-sized energies.

For the fine-scale cap use the full path sensitivity already proved for W:

    beta_W(s)=A e^-s[3+3As+(As)^2/2].

The standardized midpoint roots have supported OU regression hats with amplitude sqrt(tanh(d_k/2)). Conditional Gaussian Poincare therefore gives

    Cov(H_k^W|S_k)<=C A^2 d_k^2 I_(2n).

Hoeffding orthogonality makes its high-order covariance smaller in Loewner order, so its one-block HS norm is at most C A^2 d_k^2 sqrt(n), absorbing sqrt(2) in C.

With a physically block-separable g, all slots inside a physical block share that block's actual future, and distinct physical blocks are conditionally independent with centered high-order terms. The full covariance is block diagonal after grouping (F2^(j),Delta3^(j)). Squaring the per-block HS estimates and using sum n_j=D, sum n_j^4<=b^3D proves (1). Summing the exact same dyadic crossover d_*=A^(6/5)b^(3/5) proves (2).

## Finite level, readout and positivity

Let S_K be the last retained skeleton. The exact positive covariance identity is

    Cov(W|Y)=Cov(E[W|S_0]|Y)
       +sum_(k<K) sum_(nonempty S) E[d_(k,S)d_(k,S)*|Y]
       +E Cov(W|S_K).

Retain all |S|<=2 terms. The omitted high-order terms and residual are both PSD. Their sum bounds the full positive block error. The residual remains

    C A^2[d_K^2+(1+T^4)e^(-2T)]sqrt(2D)

in L2-HS. Applying [I I] yields a positive approximation to Cov(F3|Y); its covariance-error norm costs ||[I I]||^2=2. A buffered LAW readout, if separately supplied, costs sqrt(2) in W2 and changes its buffer exactly as specified by that actual law.

No covariance-increment positivity is inferred. The approximation retains the original coherent cross block rather than adding an independent Delta3-only covariance. Its analytic order cutoff is fixed at TWO, but the literal retained coefficient count remains order 4^K and the full coarse-skeleton/environment source is not supplied here.

## Independent check of the F2 input proof

The extension's author independently reviewed the companion proof. The substantive checks are: the entire future is decomposed relative to the single Gaussian X_t while all residuals are fixed; its regression coefficient on each later site is in [0,1] and zero after the fine-cell endpoint; entropy transport uses only ||D F1||<=A; mollification changes every ancestral source coherently; the nonlinear smooth first-order term is exactly at most pairwise in independent coarse BRIDGE objects; and the high-order projection is used only in conditional L2. These checks avoid both an illicit expectation-through-g step and a false L4 contraction of Hoeffding projections.
