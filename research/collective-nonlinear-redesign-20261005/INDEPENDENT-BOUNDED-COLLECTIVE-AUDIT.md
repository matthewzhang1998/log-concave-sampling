# Independent audit of the bounded two-level collective source

2026-10-05. Scope: the finite two-level construction and strong-error proof in `STRATIFIED-COLLECTIVE-HISTORY.md`, and its conditional interface in `COLLECTIVE-HISTORY-COMPILER-PORTS.md`. No arbitrary-depth extension, executed native service, or full compiler recertification is covered.

## Verdict and fixes

The strong approximation, linear source-occurrence count, whole-bank operator bounds, positive-mixture transfer, conditional energy envelope, and numerical VALUE floor are valid. Two small corrections are needed:

- Explicitly require t_l∈[0,1] and nonnegative normalized outer weights, in addition to the exact first moment. The pinned rule may already supply these assumptions, but the interface should state them.
- The imported bracket ell_E(a_seed+mu)+ell_E³(1+mu^(-1/2)) evaluates to 9A²+27A³+27A^(5/2) at (ell_E,a_seed,mu)=(3A,2A,A). The printed 16A² expression is a safe upper allowance rather than the exact substitution.

## Proof checks

Conditional on the entire realized OU path, clocks in distinct strata are independent. This makes the stratum errors centered and independent even though the values on the OU path are strongly correlated. Conditioning further on U_i retains that fact for its later-bin suffix; the multiplier exp(U_i) can then be bounded uniformly inside bin i. These are valid uses of variance cancellation.

In reverse bin index k, (w^(k))² eta_k≤9k/M⁴. The weighted same-bin remainder is at most (32/3)A√D/M^(3/2), and the weighted suffix remainder is at most 9A√D/M. Therefore the derived terminal bound

    A²√D/M [3+12A+32A/(3√M)]

is valid and is bounded by the advertised 24A²(1+A)√D/M. The proof needs no second derivative and no pointwise bounded random force. Each suffix is evaluated in one backward pass, so the graph has 2M+1 original VALUE occurrences and MD private Gaussian coordinates before the outer rule.

For fixed clocks, each sampled OU row has whole-bank operator norm at most one. Both collective averages have absolute weight sum one. This proves the bank Lipschitz constants without a hidden factor M. The G-block leading residual derivative is symmetric; the remaining diagonal skew and off-diagonal blocks give (2 beta+1)A²(1+A). At A=1/2 the half-variance-normalized first coefficient is approximately 2.3914 and curl coefficient 5.7956, safely below declarations 3 and 6.

The residual conditional envelope A²(1+A)(|z|/2+kappa_p√D) and caller-origin bound follow from the actual row coefficients and exact outer first moment. The terminal VALUE floor is (1+A+A²)nu; adding the residual baseline gives (2+A+A²)nu. Captures, replays and numerical clock/row errors remain separately charged.

The deterministic outer-rule error splits as R1(mu_U-psi)+(Q_out-R1)mu_U and is correctly bounded using the strong error and the uniform source norm. A single frozen random outer clock would not have this property, and the packet correctly excludes that shortcut. Mixing positive conditional service laws over the clock labels is legitimate, assuming measurable conditional services and the displayed integrated conditional allowance.

The complete interface remains conditional on all native guards and fully expanded occurrence counts. Increasing M only reduces the history approximation; it does not remove a fixed-order native completion floor. The file states this limitation correctly.
