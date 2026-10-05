# Source provenance and preserved snapshots

2026-10-05. Both requested input manifests have their exact requested hashes:

- dyadic-prefix-law-join-20261005/MANIFEST.json: 130131317183c6f00ea3566c187296912db7850081904322215a3ec9d5eb37f2
- positive-endpoint-mean-join-20261005/MANIFEST.json: af11e79e0b121f38ae9ce3559e9cd6d813d61e1e13a44514622e96e7ea1638f5

An independent verification found one stale child-file entry inside the old endpoint manifest. Its entry for independent-audit/INDEPENDENT-ENDPOINT-JOIN-AUDIT.md records cc4f5bc8d50672a8da492cef001f498f0c9275c837261d59734470a7354f89ef. The actual report hashes to 962eb398c771a25ee4fa1ddeb4a30123fef46a1ea71ab0269031c8c4664c8b10, which agrees with that independent-audit directory's own SHA256SUMS. All other file entries in that old manifest, including the mathematical endpoint source, match. All new dyadic-prefix manifest file entries match.

No source file or manifest was altered to hide or repair this mismatch. INPUT-PINS.json pins the actual old audit report separately, together with the unchanged mathematical endpoint source, the requested old manifest, and the new dyadic-prefix source and its manifest. The new endpoint theorem does not infer audit status solely from the stale old manifest entry: it has its own independent mathematical and exact-snapshot audit.

The proof uses the source-qualified native interfaces and their literal guards. It does not claim that the imported complete native compilers were numerically executed, that a production arbitrary-precision scalar Gauss solver was run, or that unspecified numerical A,D,precision values have an instantiated guard certificate.
