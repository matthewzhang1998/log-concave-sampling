# Exact joint terminal/twin gradient feedback lift

New source lemma for the third-iterate native port c8b20bf4. It retains BOTH
terminal queries on their exact joint input. It supplies a genuine-gradient
coefficient source, not a law allowing the outer ancestor to be independently
refreshed.

Let S,Z in R^d be independent standard, let g_t=grad V_t have ||Dg_t||<=1
and g_t(0)=0, and fix 0<epsilon<=1, a>0. Define

    G_Delta(S,Z)
       =a( g_t(S)-g_t(S+epsilon Z), -epsilon g_t(S+epsilon Z) ). (1)

This is the exact gradient in the COMPLETE 2d input of

    a[V_t(S)-V_t(S+epsilon Z)].

The potential is analytical only: execution uses the two displayed original
gradient VALUES, once each, plus known scalar/readout operations. They use the
SAME S,Z and finite original source version. With Delta=g_t(S)-g_t(S+epsilon Z),
its first block is a Delta, so the K3 ancestor feedback is M* times that block.

For the original master record (W,Z), S=PW, the exact gradient pullback is

    Ghat_Delta(W,Z)
       =a( P*[g_t(PW)-g_t(PW+epsilon Z)],
                                      -epsilon g_t(PW+epsilon Z) ). (2)

Since PP*=I, the known contraction [M*P,0] reads the same feedback from(2).
Known noncoisometric output maps require their ordinary fill/row certificate
when used inside a stationary channel.

## Actual norm and caller rows

Put L1=[I,0], L2=[I,epsilon I]. Then

    D G_Delta
      =a[L1* Ht(S)L1-L2* Ht(S+epsilon Z)L2],             (3)

which is symmetric. Hence its full first is at most a(2+epsilon^2)<=3a. Its
Z first is at most a epsilon sqrt(1+epsilon^2), and its S first is O(a).
These are exact finite first paths using original HVPs; no Hessian difference
is assumed small merely because epsilon is small.

Its full VALUE energy has the intact small terminal-difference scale:

    |G_Delta(S,Z)|
       <=a epsilon[|Z|+|S+epsilon Z|],
    ||G_Delta||_Lp<=C_p a epsilon sqrt(d).               (4)

The same bound holds for the master pullback. If S is captured rather than
integrated, keep its actual additional a epsilon |g_t(S)| profile. A full
conditional energy independent of S is not asserted. For a varying original
g_t label with direct caller bound L_t, the complete caller is O(a L_t),
including both terminal calls. Any actual S caller is composed with the
O(a) S first. The small epsilon VALUE bound is not differentiated.

The lower companion in(1) is essential for symmetry. At Z=0 it is
-epsilon a g_t(S), generally nonzero. Removing it to force conditional zero
would change the source and destroy the displayed full-gradient identity.
At the COMPLETE origin S=Z=0 the source is zero. Original anchor values and
their caller paths remain literal.

## Whole-input heat preserves the joint pair exactly

At a common-input OU clock substitute

    (S,Z)=c(S_root,Z_root)+s(U,V)

with U,V independent standard. The two original queries in(1) become

    c S_root+sU,
    c(S_root+epsilon Z_root)+s(U+epsilon V).

Their conditional Gaussian Gram is exactly

    s^2 [ I,I; I,(1+epsilon^2)I ],

with the corresponding captured means. In particular the second query does
not receive a newly independent additive terminal heat. The conditional mean
of(3) is a single symmetric selected matrix of an admitted gradient source.
A finite gradient-pair compiler on the full source(1) or pullback(2) therefore
extracts that WHOLE joint heat coefficient with its actual a radius, source
floors, known rows and captured-root caller bounds. Both original queries are
paid at every source call. No sampled matrix or inverse unknown Jacobian is
needed.

This identification supplies all joint feedback derivative blocks together.
It does not identify independently heated Ht(S) and Ht(S+epsilon Z) with the
original common-input coefficient under an additional nonlinear observer.

## Separate conditional-gradient port and the remaining gate

At fixed S there is also the smaller d-square gradient source

    k_S(z)=a[g_t(S)-g_t(S+epsilon z)],
    D_z k_S=-a epsilon Ht(S+epsilon z).

It has active first a epsilon and root-S first O(a). This conditional source
extracts the Z derivative only; it is not a substitute for the full S/Z
gradient port when an S derivative is required.

The actual K3 parent still evaluates

    g0(CW+M* first_block(Ghat_Delta(W,Z)))

on the SAME complete record, and its outer g_t still observes S. A stationary
law for G_Delta does not preserve that arbitrary nonlinear joint observer
merely because its marginal is standard. The next required joint-heat
comparison must keep this source genealogy or price all resulting currents.
The bound(4) is a distinct feedback energy profile; it is not asserted equal
to the smaller actual energy of the complete E_K after all terminal
cancellations. The latter remains the intact marked source where required.
