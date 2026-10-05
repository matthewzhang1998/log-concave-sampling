# Independent audit of the bounded shifted atom

2026-10-05. Verdict: PASS for the stated bounded local source, exact local positive current, and specialized native mean return. This is not a true next-history or m4 completion.

Reviewed document: `../BOUNDED-SHIFTED-CURRENT.md`, SHA256 `4c4c5047d80cc922bdd8baa3edd5b251377b362e07f2e73e6037781f6ffc0787`.

Imported LOW30 original source SHA256: `7767ac389d6b677eb9517b2585f7a24cb9f85ede47d099a653c635b7be750fb8`. The relevant source contracts are `b27:raw:self-reserve`, `b27:raw:split`, and, for its explicit caller/energy conventions, `t30:lem:value-mean`. The original source was read and hashed without modification.

## Findings

1. The graph uses exactly seven original gradient VALUES. Successive original-Lipschitz comparisons prove `|Delta|<=ab A^3|x0|` and `|E|<=abc A^4|x0|` without a derivative of a small remainder. The displayed `LU,LV`, `9A/4` difference-first bound and `21A/8` response-first bound are valid under `A<=1/2`. The retained caller has the same bounds because all four caller rows have norm at most one and the original g is frozen.

2. The complete square lift is essential and correct. With `C3=sigma P`, the term `sigma P*(J_U-J_V)P` is symmetric on the entire private tape. The remaining two terms give full curl at most `2c A(LU+LV)<=13A^2/2`, including every auxiliary off-diagonal block. No Hessian commutation or Hessian-continuity modulus is used.

3. The five-VALUE pair `(V,Delta)` gives the exact positive pushforward `Z_theta=sqrt(v)G-q(V+theta Delta)`. Its test-function current follows by ordinary differentiation in theta and the pointwise mark `-q Delta`. The SAME-tape coupling preserves Delta and all retained labels and has endpoint displacement `-q Delta`. Covariance block contraction retains both cross terms. Cumulant expansions are identities of this actual finite law; the document does not turn these identities into a higher-rank tensor supplier.

4. The direct-current first includes BOTH parts of its complete input. For the Omega/G block derivative, the operator norm is bounded by `sqrt(v+q^2 L_theta^2)`, where `L_theta=(1-theta)LV+theta LU`. The corresponding joint-map square-sum bound with the Delta first is valid. This is stronger bookkeeping than reporting only the source correction's O(A) first while omitting the fresh Gaussian root. The final clarification correctly freezes theta,q,v for these uniform firsts; its optional parameter derivatives are respectively `-q Delta`, `-(V+theta Delta)`, and `G/(2sqrt(v))`, with their actual Gaussian moment profiles.

5. At the native buffer, the actual square source is `P*E/sqrt(u)` with radius `r=O(A/sqrt(u))`, curl relative `O(A)`, and relative energy `delta=O(A^3 m_2(Y)/sqrt(n))`. Substituting `mu=r` into the FIVE raw-split remainder terms, then multiplying by sqrt(u), gives exactly:

       A^6/u, A^6/u, A^6/sqrt(u),
       A^(13/2)/u^(5/4), A^7/u^(3/2).

   The last three divided by `A^6/u` are `sqrt(u)`, `(A^2/u)^(1/4)`, and `A/sqrt(u)`. Therefore the claimed reduction under `u<=1` and `u>=c0 A^2` is valid. The physical complete/caller first is correctly scaled as `Lambda[A+A^(3/2)u^(-1/4)]`. All literal small-radius, gap, numerical and finite-clock guards must still hold; `u=A^2` alone does not prove admission.

6. The zero baseline is correctly specialized at the program level: the genuine-gradient field is identically zero under a valid coisometry, so its allocated Gaussian is exact and has no gradient-prior error. This does not set the coisometry to zero or substitute `rho=0` into an incompatible generic `r<=rho` premise. The remaining self-reserve is a VALUE program, not an executed derivative source.

7. The requested origin clarification is present in the reviewed version. Executing `E(Y,0)` gives the pointwise envelope `abc A^4|D0Y+d0|`; subtracting and restoring it yields the stated anchored profile. Energy is conditional on the actual standard private bank after the caller is fixed, and profiles are squared and averaged under the actual retained-caller law. No uniform small energy at an arbitrary unbounded caller is claimed.

8. The seven local error-readout factors `[abc A^3,bc A^2,bc A^2,c A,c A,1,1]` follow from the actual directed paths. The native cost is the complete occurrence count times seven, plus the stated capture/restoration/numerical work. Private ancestors cannot be reused across changed complete source arguments. All first/adjoint calls occur at recorded original VALUE sites.

## Independent checks

`check_shifted_atom_independently.py` and `independent_shifted_atom_checks.json` are independent of the author's executable source and checker. The run passes 6,413 assertions, including 400 nonlinear shared-bank cases with globally bounded, genuinely noncommuting Hessians. It checks the full lifted derivative, all auxiliary curl blocks, private/caller firsts, pointwise marks, executed origin, fresh-root-inclusive current and joint firsts, and the pointwise current identity. Exact rational checks establish the linear four-force sign. Exact Fraction exponent arithmetic verifies all five native error powers and their reduction ratios.

The maximum directional finite-difference derivative discrepancy is `3.15e-10`; the maximum current-identity discrepancy is `4.50e-11`. These are finite diagnostics supporting the displayed proof, not execution or numerical certification of the imported LOW30 mean compiler, a proof over all nonlinear inputs, or evidence that the atom represents the true OU history.

## Scope

The strongest certified result is the actual seven-VALUE shifted source, its exact five-VALUE positive local marked law/current, and a positive buffered native mean service for that local response. The mean service does not retain the local marked genealogy. The separate exact current does. Neither supplies the missing finite coherent true-history ancestry, the full F3 block covariance law, true third/fourth cumulant packets, an m4 accuracy gain, or a new complete-cost exponent. The document states those boundaries accurately.
