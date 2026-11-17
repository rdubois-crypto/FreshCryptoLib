# Versions

## 0.1.1 (2026-11-17)

- README.md, solidity/README.md, solidity/src/README.md: new libraries, gas measured under Osaka and Prague
  (EIP-7883 triples the modexp price in Osaka), how to run the tests without sagemath.
- ecdsa_precomputed_hackmem with a 16kb table: 74 468 -> 61 962 gas under Osaka, 69 028 -> 56 522 under Prague.
- test/gen/.gitignore for node_modules and python caches.

## 0.1.0 (2026-11-17)

### solidity/src/FCL_elliptic.sol (P-256, same API)

All existing functions keep their signature and semantics. Three internal helpers are added:
`ecZZ_mulmuladd_core`, `ecZZ_smallMultipleOfG`, `ecZZ_transpose4x64`.

- `ecZZ_mulmuladd_S_asm`, `ecZZ_mulmuladd`: shared core, joint 2-bit window (2 bits of u and 2 bits of v
  per step, 15 points iG+jQ), Q-dependent entries computed in XYZZ and normalized with a single modexp
  (Montgomery batch inversion), halved doubling with 11 mulmod, mixed addition without `sub` thanks to a
  table stored as [-x, y, -y, x] and an accumulator whose sign alternates every half step.
- `ecZZ_mulmuladd_S8_extcode`, `ecZZ_mulmuladd_S8_hackmem`: scalars transposed once into 64 nibbles,
  one shift and one mask per scalar and per round instead of 8 bit extractions, same doubling and
  addition improvements. Table format unchanged.

Gas (forge 1.7.1, solc 0.8.30, optimizer 1e6 runs, EVM osaka, external call included):

| function | before | after | gain |
|---|---:|---:|---:|
| `FCL_ecdsa.ecdsa_verify` | 212 488 | 157 680 | -25.8 % |
| `ecZZ_mulmuladd_S_asm` | 207 225 | 152 305 | -26.5 % |
| `ecZZ_mulmuladd` | 208 322 | 152 475 | -26.8 % |
| `ecZZ_mulmuladd_S8_extcode` | 76 756 | 63 112 | -17.8 % |
| `ecZZ_mulmuladd_S8_hackmem` | 69 725 | 57 205 | -18.0 % |
| WebAuthn `checkSignature` (existing test) | 227 161 | 169 584 | -25.3 % |

Runtime code of a contract exposing only `ecdsa_verify`: 5 708 bytes before, 10 393 bytes after.

Fixed behaviours:
- `ecZZ_mulmuladd` returned a wrong point for Q = -G (and Q = -c.G, c in {2, 3, 1/2, 3/2, 1/3, 2/3}).
- `ecZZ_mulmuladd_S8_hackmem` returned a wrong x when the leading 8-bit window was zero
  (the loop started at bit 254 whatever the first non zero window).
- `ecZZ_mulmuladd_S8_*` ran out of gas for u = v = 0, they now return 0.

### solidity/src/FCL_babyjubjub.sol (new)

BabyJubJub (circomlib) EdDSA verification through the isomorphism Edwards -> Montgomery -> short
Weierstrass, rescaled to a = -3 (-3/a is a fourth power mod r). The FCL engine above is reused with
lazy reduction (254-bit field). API: `verifyRaw(hm, R8x, R8y, S, Ax, Ay)`, `verifyRawW(hm, R8x, R8y, S, Wx, Wy)`
with `(Wx, Wy) = pubKeyToW8(Ax, Ay)`, plus `ed2W`, `W2ed`, `ecZZ_mulmuladd`.

| function | gas |
|---|---:|
| `verifyRaw` (circomlib key A) | 138 031 |
| `verifyRawW` (precomputed W8) | 131 971 |
| `pubKeyToW8` | 6 698 |

### solidity/src/FCL_babyjubjub_edwards.sol (new)

Same `verifyRaw` on the a = -1 twisted Edwards model (x -> sqrt(-168700).x), extended coordinates,
complete formulas: 120 888 gas against 138 165 for the Weierstrass version in the same run (-12.5 %).

### Tests (solidity/tests/WebAuthn_forge/test)

- `FCL_elliptic_opt.t.sol`: python oracle (309 mulmuladd cases including degenerate keys and every rare
  path of the loop), random ECDSA, Wycheproof valid/invalid, random differential against the original
  library kept in `ref/`, S8 variants on two tables, gas.
- `FCL_babyjubjub.t.sol`: 124 circomlibjs vectors (Poseidon, MiMC7, MiMCSponge, Pedersen signers,
  perturbations, torsion components, small order and degenerate keys), isomorphism, W8 keys,
  206 mulmuladd cases, fuzzing.
- `FCL_babyjubjub_edwards.t.sol`: same vectors, cross fuzzing between both verifiers.
- `gen/`: generators of the vectors and of both BabyJubJub libraries.
