# Finite enumeration, radius guards and countertests

This supplements the fixed-P construction with an intentionally loose numerical census bound. It is not a favorable complexity estimate.

## A computable graph-list bound

For beta=1/2,gamma=1/3, set

    N_P = ceil(3P/2),
    H_P = floor(P/4),
    D_P = max(0,floor(log_2(P/4))).

A retained graph has at most N_P centers. A reverse one-cluster polynomial on n centers has at most 2^(n-2) side-insertion patterns. Its h-th finite cumulant uses at most h^h set partitions, at most 2^[h(n-2)] insertion assignments, and at most (4hn-1)!! Gaussian pairings. These are deliberately overcounts, with every pairing treated as a separate labelled term. Physical and force labels follow their actual occurrences, so no unknown tensor value is needed to choose a term.

Thus an explicit finite branching overcount is

    B_P = sum_(h=2)^H_P h^h 2^(h N_P) (4h N_P-1)!!.

For one initial graph, the number of executed correction packets is at most

    T_P = 1+B_P+...+B_P^D_P.

If H_P<2, take B_P=0. Round noninteger target cutoffs up before using this bound. Actual connectedness, parity, term cancellations and the surplus cutoff remove many terms. The checker enumerates the much smaller numerical-state overcount (n,d), not this full labelled graph list.

Allocate variance using this known overcount before the queue is executed. For example reserve a total nu, assign at most nu/T_P to each possible packet, and split each packet share equally among its n physical readout rows and a keep. Any unused shares remain Gaussian reserve rows. The inverse physical readout product is a known finite factor bounded by a power of T_P(N_P+1)/nu. When the initial list includes positive quadrature nodes, include their known finite count in T_P; its growth is a charged public-log factor.

Every numerical cumulant coefficient and original normalization constant can be enumerated with exact arithmetic and then bounded absolutely. The finite lists of C0/C3 source/filter calls, pair orders, covariance gaps, old-row injection norms and path factors determine the complete small-alpha guard. This makes the guard computable without pretending that the constants are small or that c(P)=o(P).

## Exact scalar fixture and orientation dependence

Choose the center spanning path a-b-d-c, root at a's leaf and select b's leaf as the root-spine terminal. The other selected spines are d-leaf(d) and c-leaf(c). Label the five cut edges by G0 for the other a-b edge, G1 for the other b-d edge, G2 for the other d-c edge, and G3,G4 for the two c-a edges.

With scalar unit coefficient targets, and after suppressing readout constants, the reverse polynomial is

    V/H = alpha^(d+2) z^2 theta^2 G0^2 G1 G3 G4 Yd
        + alpha^(d+3) z^3 theta^3 G0^2 G1^2 G2 G3 G4 Yc
        + alpha^(d+4) z^4 theta^4 G0^2 G1^2 G2^2 G3^2 G4^2.

All seven Gaussians G0,...,G4,Yd,Yc are independent in this reverse representation. At d=4 its first expectation is alpha^8 H z^4 theta^4. Its complete second log coefficient is

    H^2 [(3/2) alpha^12 z^4 theta^4
       + (9/2) alpha^14 z^6 theta^6
       + 121 alpha^16 z^8 theta^8].

The middle coefficient is orientation dependent. The independent auditor's exhaustive 192-orientation census finds numerator triples (3,9,242) for 128 choices and (3,18,242) for 64 choices before dividing by 2. The universal first feedback coefficient is 3/2. The Python fixture deliberately uses the one displayed orientation; it does not assert the 9/2 coefficient for every orientation.

## Countertests retained in the record

1. Equal alpha-per-force amplitudes reproduce the old same-grade alpha^8 z^4 variance. The new normalization changes the actual root-spine amplitude, not a label on that old variance.
2. Surplus, not effective exponent alone, is the queue order. At n=100,d=4 the input exponent is 212/3, while a four-center h=2 return has exponent 32/3. Executing a high-grade correction can expose a lower-grade current; one must not repair discarded terms and then drop their offspring.
3. Each first-variance diagram remains loopless 4-regular after all three G0 Wick pairings. Two-force chunks retain physical leaves in all 192 short-spine orientations. No long-spine mark-restoration claim is used.
4. Finite characteristic jets are not an executable MGF oracle and do not alone give W2 control. The separate frame-and-remainder addendum is required.
5. Arbitrary old-row norms are not silently bounded. All O(alpha) first assertions require the recorded row/path guard, including alpha dependence.
6. A finite retained nonzero polynomial log generator is a current target, not a positive probability law. A law comparison requires an eligible positive joined target and its complete coefficient cancellation conditions.
7. The original alpha^5 tau smoothing debt and the whole-old-bank/tilt/terminal-curl gates remain open outside this sector.
