# Aggregate source ports with smaller endpoint cutoffs

2026-10-05. Companion to SMOOTHED-BRIDGE-PREFIX.md. This strengthens a sufficient source-port condition; it does not supply a new skew or cumulant service.

The standalone 17/5 theorem conservatively sets eta=A to keep every node's carrier-subtracted first O(Lambda A). This is sufficient, but not necessary for the FINAL positive-sum source. In fact the source ports remain valid whenever w>=A and every individual native service remains admitted at its actual local buffer.

The imported mean/Gram graph theorem gives each node an actual residual first, including caller first,

    L_t<=Lambda[A+A^(3/2)/sqrt(v_t)].

The K graph has first <=Lambda A and is included in this envelope. Let the finite positive outer rule have normalized weights omega_t. The exact reciprocal formula and positive dyadic panel bounds give

    sum_LAW omega_t v_t^(-1/2)
        <=C[w^(-1/2)+1],

independently of the lower cutoff eta>0. Thus

    sum_LAW omega_t L_t
       <=C Lambda[A+A^(3/2)/sqrt(w)]<=C Lambda A        (1)

whenever w>=A. The RAW endpoint nodes already have first O(A).

A known orthogonal rotation aligns each node's entire terminal source-zero Gaussian row rX_t+hG_h to one common D-root. Orthogonal rotations preserve the operator first bound of each complete node map. After embedding all node tapes into one common-G-plus-complements tape, weighted triangle inequality gives (1) for the actual complete private first, and separately for caller first. There is no max_t L_t factor: all node weights remain in the literal executing graph.

The displacement from the aligned terminal carrier has actual first <=C[L_t+Aw(1+L_t)] and vanishes at the whole-source origin. Its Lp energy at fixed caller z is therefore bounded by this first times |z|+sqrt(d_t), with d_t the complete active node dimension. Since d_t/D is a fixed public-log polynomial, weighted triangle inequality and (1) give aggregate displacement energy Lambda A(|z|+sqrt(D)). Multiplying by the outer terminal A-Lipschitz g gives the required remainder energy Lambda A²(|z|+sqrt(D)). This is a direct source-graph estimate; no LAW-to-RAW inference is used.

For curl, the scalar aligned-carrier coefficient multiplies a symmetric Hessian difference at every node. The only nonsymmetric products contain terminal Dg and a displacement first. Summation with the actual positive weights gives curl <=C A sum omega_t(L_t+Aw(1+L_t))<=Lambda A². Off-common-root paths obey the same estimate. The residual's unrestricted first is Lambda A because the symmetric Hessian-difference term is O(A). The final own-mean compiler therefore sees the same genuine-gradient baseline, remainder energy and curl as in the standalone theorem, although some individual node displacement first radii exceed A.

Every reentry of this FINAL source evaluates its full positive sum before consumption and recreates its complete bank. Thus a replay does not select an exceptional near node or erase its weight. An interface that exposes a single node as an O(A)-first RAW service would need a separate bound and is NOT authorized by this aggregate argument.

## Conditional extension to a deeper cutoff

The actual native mean radius at buffer u is Lambda A/sqrt(u), while its raw self-reserve first is a public-log factor times A^(3/2)/u. Hence the cutoff

    eta=C_guard Lambda_guard A^(3/2)

is eligible only when the declared C_guard Lambda_guard is chosen to satisfy every imported numerical self-reserve/radius guard, eta<=w, and all other explicit native conditions are checked. No asymptotic equality alone certifies admission. Alternatively, eta=A^beta with any fixed 1<beta<3/2 gives vanishing power margin in these guards, subject to the usual public-log small-parameter window.

At this endpoint scale and w=A^(5/6), h=A^(7/4), the new prefix bound is

    C[A^(15/4)+A^(23/6)+A^(9/2)]sqrt(D)
        +epsilon_bridge A^(17/6)sqrt(D).

The fresh-RAW endpoint debt is a public-log factor times A^(15/4)sqrt(D). This makes the prefix compatible with a separately admitted skew service whose integrated main remainder is A⁵w^(-3/2), also grade15/4. The present companion supplies the prefix and aggregate ports only. It does not certify that skew service, its Gaussian regression, or its full joined cost; those must come from its own audited packet.
