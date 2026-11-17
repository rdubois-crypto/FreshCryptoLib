# Generators

Run from `solidity/tests/WebAuthn_forge`. Random values come from `secrets` (python) and `crypto.randomBytes` (node).

P-256 vectors:

```sh
python3 test/gen/p256.py test/vectors_opt/p256.json
python3 test/gen/p256_s8.py test/vectors_opt/p256_s8.json
python3 test/gen/p256_s8.py test/vectors_opt/p256_s8_2G.json 2G
```

BabyJubJub vectors (circomlibjs signers):

```sh
(cd test/gen && npm install)
node test/gen/gen_vectors.mjs test/vectors_bjj/eddsa.json test/vectors_bjj/keys.json
python3 test/gen/gen_w.py test/vectors_bjj/keys.json test/vectors_bjj/weier.json
```

BabyJubJub libraries (the Weierstrass one reuses the core of `FCL_elliptic.sol`):

```sh
python3 test/gen/build_bjj.py ../../src/FCL_elliptic.sol ../../src/FCL_babyjubjub.sol
python3 test/gen/build_bjj_edwards.py ../../src/FCL_babyjubjub_edwards.sol
```

`bjj.py` holds the curve maps (Edwards, Montgomery, Weierstrass a = -3) and prints the constants;
`rare_cases.py` builds scalars that drive the 2-bit joint window loop through its rare paths.
