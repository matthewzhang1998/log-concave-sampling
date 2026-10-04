# Independent whole-source OU orientation stability audit

Scoped PASS against `astra-coherent-gram/common-input-orientation-heat/WHOLE-SOURCE-OU-ORIENTATION-STABILITY.md`, SHA-256 `4b2b72a4da96f944b1981de748a4db75b22e22458f20fd6035f7b1f16e5494a6`.

For vector chaos degree k, the source's output-swap convention gives exactly delta(Curl f_k)=k Aop_k f_k. The degree-(k-1) divergence norm is at most sqrt(k), yielding a k^-1/2 curl-row bound. Orthogonality of the curl degrees then gives the row frame kappa sup_k (1-exp(-2kt))/sqrt(k), bounded by kappa min(1,sqrt(2t)). This is a sum of row energies before the supremum; it incurs no dimensional or Hermite-degree count.

OU is self-adjoint and commutes with the degreewise swap. Therefore O_f-O_(P_t f)=Sym E[(f-Ef) tensor Aop(I-P_(2t))f]. The first map has HS norm e, while the second has the row operator bound just proved. This gives the claimed one-energy HS inequality and, using Gaussian Poincare on the first factor, its ordinary A kappa operator counterpart. The weaker E[CC*] operator hypothesis suffices, as printed. Closure from finite Hermite sums is legitimate in Gaussian Sobolev L2.

I reran `check_whole_source_heat.py`: 1,080 mixed-chaos checks passed, maximum divergence-identity discrepancy below 2.3e-16. This confirms the displayed normalization in finite cases and supplements the proof.

The heat acts on the WHOLE same finite input graph, including all twins, anchors and repeated query aliases. The result supplies an analytical target-restoration tolerance only. It does not produce a finite averaged source, preserve the coefficient after independently heating primitive Hessians, or supply its actual retained/caller contract. Its use at a consumer also requires the stated fixed covariance gap or the correct same-endpoint current.
