# Independent literal-source audit of the carrier rescaling cut

2026-10-04. Checked LITERAL-K3-CARRIER-RESCALING-CUT-FAILURE.md against the actual LOW30 formulas, not only its numerical checker.

The calculation is sound within its stated shortcut scope. The exact scalar native source with g(y)=m y is E=bZ, b=ra^2 m^3 epsilon. Its square-lift first is L=b eS eZ*, with L^2=0 and LL*=b^2 eS eS*.

LOW30 lines 208–221 define the finite first-response filter with weights summing to one. For a linear source every response is exactly Lp. Lines 253–270 normalize the clock source to Lz/ell and form the same-tape outer difference; this is exactly L^2p and therefore zero at every finite node and every width. Scaling by alpha preserves the zero. No limit, unknown matrix oracle or derivative-as-VALUE is involved.

LOW30 lines 331–340 then give the literal nine-replay mean

    sum_i wi f(Ui)+eta P−(q2/(2eta))Ccov(f;P)+zeta Z.

Since Ccov=0 on this source, its physical output is exactly C+alpha b sum_i wi Zi, with C=eta P+zeta Z standard and independent of the nine main records. The stated weights have sum squares 1/3. Thus the two claimed scalar Gaussian W2 errors are exact:

    before: sqrt(1+alpha^2 b^2/3)−1,
    after C+(Yalpha−C)/alpha: sqrt(1+b^2/3)−1.

For 0<b<=1, rationalization gives

    error=b^2/[3(sqrt(1+b^2/3)+1)],

which lies between b^2/7 and b^2/6. At r=a=A,epsilon=A^.9 it is indeed Theta(A^7.8), with the fixed m-dependent coefficient. The full residual first and L2 energy are alpha b/sqrt(3) before and b/sqrt(3) after the inverse readout.

The exact source-zero carrier relation is an internal packet wire. The complete packet law can of course be compared after execution; it is then exactly the unchanged law above. Applying a marginal comparison before subtracting its actual carrier does not give a stronger result. Independent full outer packets retain the same within-packet covariance.

The source note now explicitly uses 0<alpha<=1 for attenuation. This is the minor guard clarification requested in the first audit and it has been applied; the underlying original source still obeys the stated seed guard. Numerical floors and all scheduled original calls still need their stated charges.

This does not rule out a genuine covariance repair, an explicit solution of this quadratic fixture, or other no-copy algorithms. It refutes the specific finished-output attenuation / old-carrier rescaling shortcut.
