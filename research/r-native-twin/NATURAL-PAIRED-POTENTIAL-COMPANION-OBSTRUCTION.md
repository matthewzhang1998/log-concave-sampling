# Natural paired-potential lift: companion energy and first obstruction

New source-qualified test, 2026-10-04. It concerns the natural whole-composite gradient lift of the ACTUAL scalar K3 source. It does not rule out a different VALUE-only lift, a retained-body completion, or a weaker current comparison.

## 1. One original strongly convex primitive and exact companion

Use the same original gradient at both nodes,

    g(y)=m y+(d/k)sin(ky),       m=1/2, d=1/10,
    H(y)=m+d cos(ky),           .4<=H<=.6=:L.

This is anchored and globally strongly convex with uniformly bounded Hessian, for any frequency k. Let 0<a<=1/16, 0<epsilon<=1, r>0, and k>=1/(a epsilon). Set independent standard S,U,Z, x=cS+sU with c^2+s^2=1, and use M=1. The exact native queries are

 Delta=g(S)-g(S+epsilon Z),  x1=x+a Delta,
 t0=S+a g(x),               t1=S+a g(x1),
 E=r[g(t0)-g(t1)].

For independent analytical coordinates (S,x,Z), consider the paired potential

    Psi=r[V(t0)-V(t1)].

Its x companion is EXACTLY

    C_x=ra[H(x)g(t0)-H(x1)g(t1)]
       =ra{[H(x)-H(x1)]g(t0)+H(x1)[g(t0)-g(t1)]}.     (1)

The known linear change x=cS+sU restores the actual Gaussian law. It only changes these companion rows by fixed known factors and does not erase the bound below.

The actual native mark satisfies, pointwise,

    |E|<=L^3 r a^2 epsilon |Z|,
    e=||E-EE||_2<=L^3 r a^2 epsilon.                   (2)

We do not replace e by the energy of an individual terminal.

## 2. A finite, uniform lower bound for companion energy

Put phi=kx, t=m k a epsilon, and

    eta=a d[sin(kS)-sin(kS+k epsilon Z)],  |eta|<=2ad.

Then exactly k(x1-x)=-tZ+eta. Define

    L0=d m S[cos(phi)-cos(phi-tZ)].

By the Lipschitz bound of cosine,

    |[H(x)-H(x1)]-d[cos(phi)-cos(phi-tZ)]|<=2ad^2.

Also ||g(t0)-mS||_2<=m a L+d/k<=.4a. The second term of (1), divided by ra, is at most L^4 a^2 epsilon in L2. Therefore

    ||C_x/(ra)-L0||_2<=.09a+.1296a^2 epsilon.           (3)

The comparison keeps the exact phase perturbation eta from the SAME terminal/twin pair; it is not dropped as an independently sampled error.

Gaussian integration gives an exact formula for the leading norm. Since x is standard, Corr(S,x)=c and Z is independent,

 E L0^2=d^2 m^2{1-exp(-t^2/2)
       +[1-2exp(-t^2/2)+exp(-2t^2)]/2
                   *(1-4c^2 k^2)exp(-2k^2)}.         (4)

Here t>=m=1/2 and k>=16. The exponentially small final term in braces has absolute value below10^-200, uniformly in c. Thus ||L0||_2>.017. For a<=1/16 the right side of (3) is at most .006132. In particular the explicit uniform lower bound

    ||C_x||_2 >= .01 r a                              (5)

holds throughout this parameter range. At the weak native scales r=a=A, epsilon=A^.9, its ratio to the ACTUAL mark satisfies

    ||C_x||_2/e >= [.01/(.6)^3]/(a epsilon).

This diverges as A tends to zero. The full potential gradient therefore does not have an O(e) one-energy return, let alone an O(kappa e) return, merely because its leading physical VALUE difference is E. The x companion is a genuine additional unmarked body of scale ra.

## 3. Its active first is not controlled by the admitted first-only source class

At fixed a,epsilon, allow k to increase while keeping the SAME uniform original Hessian bound[.4,.6]. Choose S=1, x=pi/(2k), and Z=pi/(m k a epsilon). Then phi=pi/2 and k(x1-x)=-pi+eta with |eta|<=2ad. At these actual records, g(t0) and g(t1) tend to m, while H'(x)=-dk and H'(x1)=dk cos(eta). Differentiating the exact companion (1) in x yields

    (ra)^(-1) D_x C_x
       =H'(x)g(t0)-H'(x1)g(t1)
            +a H(x)^2 H(t0)-a H(x1)^2 H(t1).

The last two terms are uniformly bounded by2aL^3. The first two have magnitude at least c0 k for all sufficiently large k, with a numerical c0>0 independent of k. Hence the companion's complete first can grow like ra k with no uniform bound determined by r,a,epsilon and the original first radius alone.

Equivalently, this derivative differentiates an original HVP and introduces the uncontrolled third derivative of V. The analytical paired potential is a valid scalar function, but its full gradient is not thereby an admitted VALUE-only primitive with the required first/caller contract. Executing the displayed companion itself already uses original Hessian actions; a later first sweep would need more than the permitted original first information.

## 4. What this does and does not settle

The result excludes treating this PARTICULAR natural paired-potential gradient as a free one-mark producer for the K3 remainder. It also explains why listing its scalar potential is not enough to close a repair generator: the companions and their full firsts must be implemented and priced.

It does not exclude retaining a larger unmarked body with a separate valid first/caller implementation, a different source-aware joint-gradient lift, or an analytical current in which those companions cancel BEFORE a norm is taken. The scalar direct true-Gram reserve3e9fb810 does not use this lift and is unaffected. The general matrix same-query word36d25b9e remains a distinct construction problem.
