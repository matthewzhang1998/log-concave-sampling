# A positive projected-gradient split for the terminal-flat fixture

New finite source construction, 2026-10-04. This is a positive escape from using the terminal-flat fixture as an obstruction to every multi-source repair. It is not a general bounded-Hessian decomposition theorem. Its positive energy gain uses the explicit fixture's narrow nonlinear support.

## Result first

For the exact terminal-flat SAME potential, there is a two-original-VALUE full gradient G0, a known coisometry B0, and an actual seven-VALUE remainder R0 such that

    E=B0 G0+R0,
    Lip G0<=C ra,
    ||G0||p<=C_p e_actual,
    ||R0||p<=C_p[a+(a epsilon)^(1/p)] e_actual.

Here e_actual is the centered L2 energy of the ORIGINAL complete E source, not a nominal or independently sampled substitute. The source B0G0 reproduces the original E value, full first, missing two-factor word and physical skew exactly at the terminal-flat selected record.

A fixed-angle rotation applied ONLY to the projected-gradient approximation creates genuine new auxiliary V columns. Because the linear part is unchanged and the nonlinear part is supported on narrow slabs, this also has an actual-E-marked source error of order (a epsilon)^(1/p)e_actual. The original W,V and the entire old nonlinear opposite baseline can remain unchanged and retained in this pathwise comparison. This is not a common reparameterization of both sides.

The remainder has actual first O(a), recorded curl O(a^2), and the positive energy gain above, so it genuinely re-enters the unchanged near-gradient envelope. This provides a fixed-grade no-copy source split and an actual V-column escape on this fixture. It does not supply an arbitrary-rank decomposition for a general potential, nor automatically admit an observer of fine records already integrated by a pair theorem.

## 1. Literal fixture, scales and actual energy

Set d=2, M=I, c=3/5, s=4/5,

    T=[[1/2,1/20],[1/20,1/2]], delta=1/40,
    r=a, epsilon=a^(9/10), q=a epsilon,
    kappa0=ra,
    g(y)=Ty+u(y),
    u(y)=delta q psi(y1/q)e1,
    psi(z)=(z+1/2) chi(10(z+1/2)).

The smooth chi is supported on [-1,1] and equals exp(1-1/(1-t^2)) inside. Thus u is supported where y1 belongs to I_q=[-.6q,-.4q]. As already verified, .4I<=Dg<=.6I, g(0)=0, |u|<=Cq, Lip u<=C, and ||T||=.55<c. Every occurrence uses this SAME primitive.

The actual K3 source on W=(S,U,Z) is

    x=cS+sU,
    Delta=g(S)-g(S+epsilon Z),
    x1=x+a Delta,
    t0=S+a g(x), t1=S+a g(x1),
    E=r[g(t0)-g(t1)].

Its complete source is unchanged. Under the nondegenerate native auxiliary bound (M=I and numerical strong monotonicity),

    c_* r a^2 epsilon sqrt(d) <= e_actual
       :=||E-EE||2 <=C r a^2 epsilon sqrt(d).             (1)

The lower bound is uniform in the other original roots before integrating them. The original pointwise bound |E|<=C r a^2 epsilon |Z| supplies the fixed-p upper envelopes. Consequently one may use (1) to convert the explicit nominal bounds below into one actual-E mark. This step is justified for this fixture; it is not asserted for arbitrary singular M.

## 2. A two-VALUE full-gradient source

Define known d-by-3d maps

    C0=[cI,sI,0],
    C1=[cI,sI,-qT],
    xbar1=C1W=x-qTZ.

They obey C0 C0*=C0 C1*=I. On the original W, define

    GDelta(W)=kappa0[C0* g(C0W)-C1* g(C1W)].             (2)

This is the exact full gradient of kappa0[V(C0W)-V(C1W)]. Potential VALUES are not queried. It uses two original gradient VALUES, and its first/adjoint uses the same two original HVP sites. In particular

    Lip GDelta<=kappa0(.6)(||C0||^2+||C1||^2)<=C kappa0.

Add one unused d-dimensional Gaussian pad and put

    G0=(GDelta,0),
    B0=[T C0, (I-T^2)^(1/2)].

Then B0B0*=I, with a fixed known spectral gap in the fill, and

    E0:=B0 G0=kappa0 T[g(x)-g(xbar1)].                 (3)

There is no inverse q, a, singular value, or small covariance gap. The known T is part of this explicit fixture. It is not a secretly queried unknown Hessian of a general potential.

The complete gradient source has the useful exact decomposition

    GDelta=kappa0{C0*[g(x)-g(xbar1)]
                         +q Zrow* T g(xbar1)},         (4)

where Zrow selects Z. The second term is a real source companion and is retained. Original Lipschitz and Gaussian moment bounds give

    ||G0||p<=C_p kappa0 q sqrt(d)<=C_p e_actual.        (5)

Thus the full source, not just its physical readout, has an actual-E energy certificate. The raw individual gradient terms in (2) need not have that mark; (4) proves it for their SAME-record difference.

## 3. Exact remainder and a positive fixed-moment energy gain

Write

    eta=u(S)-u(S+epsilon Z),
    Delta=-epsilon TZ+eta,
    x1=xbar1+a eta.

The exact finite remainder is

    R0=E-E0=R_feedback+R_terminal,
    R_feedback=ra T[g(xbar1)-g(x1)],
    R_terminal=r[u(t0)-u(t1)].                          (6)

There is no law subtraction in (6). A joint E/E0 evaluation shares g(x) and needs one additional g(xbar1) call beyond E's six original VALUES: seven generic distinct VALUES, plus their actual first/adjoint replay. At the special selected record further aliases may occur, but they are not assumed in this generic count.

Since |eta|<=Cq,

    ||R_feedback||p<=C_p r a^2 q<=C_p a e_actual.       (7)

For the terminal term, global Lipschitz gives

    |R_terminal|<=C r a^2 epsilon |Z|
                   1{t0,1 in I_q or t1,1 in I_q}.     (8)

Condition on (S2,U,Z). The variable S1 remains standard Gaussian. Both maps S1 -> t_i,1 are strictly increasing, with derivatives bounded below by a fixed positive constant for a<=1/16:

    d t0,1/dS1=1+a c H11(x)>=1,
    d t1,1/dS1
      =1+a e1*H(x1)[c e1+a(H(S)-H(S+epsilon Z))e1]
      >=1-2a^2(.6)^2>=1/2.

The second lower bound uses the positive a c H11 term rather than discarding its sign. Hence the conditional probability of either slab event is at most Cq, uniformly in the conditioned U,Z,S2. Applying (8) before integrating the fixed |Z|^p gives

    ||R_terminal||p<=C_p r a^2 epsilon q^(1/p)
                       <=C_p q^(1/p)e_actual.         (9)

Combining (7)-(9),

    ||R0||p<=C_p[a+q^(1/p)]e_actual.                   (10)

For p=2 this is O(A^.95 e_actual), since q=A^1.9 and a=A. For an exterior finite moment list up to p_*, take the actual common gain min(1,1.9/p_*), allowing strict slack for logarithmic factors. An all-p gain independent of p is not claimed.

The actual full first is bounded directly from the two executed programs:

    Lip R0<=C r+C kappa0<=C A.

The original E square lift has curl O(kappa0), while E0 has first O(kappa0); hence the recorded-P_S curl of R0 is O(kappa0). No derivative of (10) is taken. Both sources vanish at the complete fresh origin. For this fixed-data fixture the direct external caller is zero; varying supplied labels require their literal original query paths, not a reset inferred from the energy gain.

Thus R0 is an actual new near-gradient VALUE source in the unchanged A,kappa envelope, with a positive actual energy gain. The identity does not say that applying the same split again gives another gain: the remaining localized nonlinear source has not been recursively decomposed here.

## 4. The lowest terminal-flat word and physical skew are reproduced exactly

At the original selected record S=U=0,Z=e1,

    x1=xbar1=-qT e1,
    H(0)=H(epsilon e1)=T,
    A0=T, A1=T+delta N, N=e1e1*,
    T0=T1=T.

The affine-feedback and original-feedback derivatives also agree at this record. Therefore

    E0=E=ra^2 epsilon T^3 e1,
    D E0=D E,
    D_S E0=-kappa0 c delta TN,
    D_U E0=-kappa0 s delta TN,
    D_Z E0=kappa0 q T(T+delta N)T.                     (11)

In particular R0 and its full first vanish at this selected record. Equation (11) gives the noncommuting physical skew, not only the hidden-column square. With Bword=-kappa0 delta TN, the leading physical orientation is

    Bword Bword* - c^2 Sym(Bword^2),

whose (2,2) entry is (kappa0 delta)^2/400. The small projected-gradient lift (2) has reproduced precisely the term that the common-mode/parity test missed. Thus that test is not an impossibility proof against this source decomposition.

## 5. A genuine new auxiliary V column at a fixed angle

Adjoin an independent standard V, while keeping the original W and its labels unchanged. Fix a numerical angle theta with sin(theta)!=0 and cos(theta)!=0. Let

    U_plus=cos(theta)U+sin(theta)V,
    U_minus=cos(theta)U-sin(theta)V,
    H_plus=E0(S,U_plus,Z), H_minus=E0(S,U_minus,Z).

These are new functions on the SAME master tape (W,V). They are not obtained by simultaneously rotating the old baseline or the original E. Each has an explicit two-VALUE full-gradient lift using

    C_plus,0=[cI,s cos(theta)I,0,s sin(theta)I],
    C_minus,0=[cI,s cos(theta)I,0,-s sin(theta)I],
    C_sign,1=C_sign,0-qT Zrow,
    G_sign=kappa0[C_sign,0* g(C_sign,0 Wpad)
                             -C_sign,1* g(C_sign,1 Wpad)].

The known projections B_sign=[T C_sign,0, (I-T^2)^(1/2)] differ. Both have a numerical fill gap, full source first O(kappa0), and full source energy O(e_actual). The original Gaussian W and the auxiliary V remain distinct owned records.

A generic Lipschitz rotation bound would lose the energy advantage, but this fixture has a sharper exact expansion:

    E0(W)=kappa0 q T^3 Z
        +kappa0 T[u(x)-u(x-qTZ)].                     (12)

The first term in (12) is unchanged by rotating U. The remaining terms have amplitude O(kappa0 q), and each of their first-coordinate slab arguments is a nondegenerate Gaussian with uniformly bounded conditional density given Z. The same slab argument therefore proves, even at a fixed angle,

    ||H_sign-E0||p<=C_p kappa0 q^(1+1/p)
                           <=C_p q^(1/p)e_actual.     (13)

Together with (10), this gives an exact actual-E-marked remainder

    E=H_sign+R_sign,
    ||R_sign||p<=C_p[a+q^(1/p)]e_actual,
    Lip R_sign<=C A, Curl(P_S*R_sign)<=C kappa0.       (14)

Unlike a mere common coordinate change, H_sign has an actual nonzero V column relative to the original owned W. At the terminal-flat record U=V=0, its S derivative is unchanged, while the U/V derivatives split the old hidden row by cos(theta), +/-sin(theta). Their squared sum is the old hidden Gram, so its full leading physical orientation at that record is unchanged. The individual derivatives of R_sign need not be small; (14) is not differentiated.

## 6. Jointly retain the unmarked baseline and price the endpoint comparison

Let N_h(W,V) be the entire old nonlinear midpoint opposite baseline, with all six original queries and the same W,V. The tuples

    (W,V,N_h(W,V),E(W))
    (W,V,N_h(W,V),H_sign(W,V))

have a literal SAME-record coupling in which only their last coordinate differs. Equation (14) supplies its Lp discrepancy. No old root, nonlinear baseline body, Gaussian carrier or query alias is replaced. The baseline still has its actual unmarked O(kappa0) Gaussian row and its actual kappa0/h root/caller bill.

For an exterior observer uniformly L-Lipschitz in that last coordinate, the coupled error is at most L times (14), regardless of its dependence on the unchanged baseline. This is an actual retained-source comparison. It does not silently give the stronger kappa0 times that error for every arbitrary endpoint observer. Such an extra factor must come from the actual consuming row/current.

There is also a useful RETAINED COVARIANCE-FIELD comparison at a supplied common Gaussian covariance clock. Define, for a source F and fixed physical S coisometry P,

    O_J(F;R)=sum_j w_j {J_F,j J_F,j*
                                  -Sym[(J_F,j P*)^2]}.

Use the same literal clock and same roots for E and H_sign, with complete independent fine banks where its products require them. Since H_sign has full first O(kappa0), the residual has curl O(kappa0), and the actual clock sums satisfy sum w_j/v_j<=Lambda, polarization and first-chaos Bessel give

    ||O_J(E)-O_J(H_sign)||_(L2 roots;HS)
       <=Lambda kappa0 ||R_sign||2
       <=Lambda kappa0[a+q^(1/2)]e_actual.             (15)

Indeed O_J(R_sign) is bounded by its curl times its Hilbert mark; each cross term with H_sign has its O(kappa0) operator first on one side and the R_sign Hilbert mark on the other. This is a norm of the actual retained field, not a comparison only after expectation. Fixed-gap conditional Gaussian references differing by this field have the corresponding same-root W2 bound, provided their required positive covariance references have actually been supplied.

In (15), old observers may read the retained coarse roots and baseline values computed from those roots. It does NOT permit appending private fine records that a pair comparison integrated away. Retaining N_h at a finer source record inside that comparison requires the actual additional joint blocks and matched current, which are not supplied by (15). The exact tuple comparison above remains valid for the full original W,V, with its stated energy bill.

## 7. Finite cost and the remaining all-rank boundary

Every call is an original-gradient VALUE program. E0 uses two original VALUES; joint E/E0 uses seven generic VALUES. Each rotated source H_sign uses two at its own actual affine arguments. A joint E, H_plus, H_minus and old baseline call has at most twelve distinct original VALUES before valid exact-key aliases (six for E, four rotated ancestor queries, two additional baseline terminal queries). Source restoration, all known rows/fills, moving callers and any changed source graph are included at their actual versions. One first/adjoint direction uses the corresponding original HVP sites and known row sweeps only.

Known fills have numerical spectral gaps; there is no inverse-q readout in the gradient split or in the fixed-angle source. All finite numerical errors enter additively after the actual output/readout multipliers. To claim an actual-E-marked tolerance, use (1) or retain an absolute floor. Compiler orders, clock cardinalities, share inverses and precision schedules remain fixed-order polynomial-logarithmic gates. A fixed-rank invocation therefore has extra original-query/first-sweep exponent b=0, not an inverse-heat empirical count.

At p=2 the nominal source profile e_actual=Theta(A^3.9) makes the field tolerance in (15) order A^6.85, compared with the original A^5.9 envelope. This is a bounded fixture-specific improvement, with its absolute numerical floors. The highest requested moment can require the smaller gain 1.9/p_*.

The positive result uses a supplied known T and a nonlinear component supported on an interval of width q. A general Hessian-sandwich potential need not have either this affine decomposition with a small supported remainder or a uniform gain replacing (9)/(13). The remaining R_sign is an admitted near-gradient VALUE source, but it is not an original gradient and no recursive decomposition of it is proved here. A stationary skew/Gram channel for H_sign must still retain its own required joint baseline blocks; a marginal gradient mean is not a substitute for that endpoint contract.

Thus the actual-V escape passes the terminal-flat test and produces a real fixed-grade re-entry. It neither proves the all-rank generic compiler nor supports using the earlier fixed-frame/parity results as a no-go against all finite differently projected sources.
