# Architecture decisions

## ADR-1: joint 2-bit window for uG+vQ, table normalized with one inversion (2026-11-17)

Context: the 1-bit Strauss-Shamir loop of `ecZZ_mulmuladd_S_asm` performs 256 doublings and about 192
mixed additions, plus three `if` per bit to select G, Q or G+Q.

Decision: process 2 bits of u and 2 bits of v per step with the 15 points iG+jQ. G, 2G, 3G are constants;
the 12 Q-dependent points are computed in XYZZ (1 doubling, 10 additions) then normalized together
(Montgomery batch inversion, one modexp). Table entry k sits at T + 128k as [-x, y, -y, x], so the
digit is turned into an address with shifts only.

Consequences: about 120 additions instead of 192 and no selection branches. About 7 k gas of
precomputation. Larger joint windows (3 bits, 63 points) cost more in precomputation than they save.
Runtime code grows (5.7 kB to 10.4 kB for a contract with only `ecdsa_verify`).

## ADR-2: halved doubling and sign alternating accumulator (2026-11-17)

Decision: doubling with U = Y instead of U = 2Y, output scaled by 1/2 in XYZZ (3/2 constant in M), and
W.Y replaced by V^2: 11 mulmod. The FCL "-Y" trick is extended to the addition: the accumulator
represents -A, the table stores -x, so P and R come without `sub`; the result sign alternates every
half step, the first half of a 4-bit group adds -T[d] (entry offset 0x40), the second adds T[d] (0x20).

## ADR-3: small multiple keys handled by rewriting, not by complete formulas (2026-11-17)

Context: table construction degenerates when jQ = +/-iG, i.e. Q = +/-c.G with c in
{1, 2, 3, 1/2, 3/2, 1/3, 2/3}.

Decision: compare Qx with the 7 x-coordinates at entry (`ecZZ_smallMultipleOfG`) and rewrite
uG+vQ = (u + c.v)G. Table additions then never meet P = +/-Q and need no branch. In the main loop the
only check is P == 0 (same x); the point at infinity is kept with X = 0 so that it is caught by the same
test. The original `ecZZ_mulmuladd` was wrong on these keys, the original `ecZZ_mulmuladd_S_asm`
only handled Q = -G.

## ADR-4: S8 variants, scalars transposed into nibbles (2026-11-17)

Decision: the 8-dimensional tables (sage `Precalc_8dim` format) are kept. Each scalar is transposed once
(nibble j = bits j, 64+j, 128+j, 192+j) with a 6-step spread of each 64-bit limb, so the table offset of a
round is one shift and one mask per scalar. Doubling and addition as in ADR-2, with the table point
positive. The extcode variant reads into the scratch space 0x00-0x3f to stay under 16 stack slots.

## ADR-5: BabyJubJub through a = -3 short Weierstrass isomorphism (2026-11-17)

Context: same approach as ed25519 in SCL (Wei25519), reuse of the a = -3 FCL engine.

Decision: map Edwards (168700, 168696) to Montgomery (A = 168698, B = 1), then to short Weierstrass,
then rescale (x, y) -> (mu^2 x, mu^3 y) with mu^4 = -3/a_w, which exists since -3/a_w is a fourth
power mod r. The doubling of FCL applies unchanged. Since r < 2^254, values are kept below 4p between
mulmod and most addmod become add (bounds: X < 3p, Y < 2p after an addition).

## ADR-6: raw verification API and projective check of R8 (2026-11-17)

Decision: `verifyRaw` takes the already hashed challenge hm, so one function covers the Poseidon, MiMC7,
MiMCSponge and Pedersen signers of circomlib. It checks S.B8 - hm.(8A) = R8 with v = n - (hm mod n)
(8A lies in the prime subgroup), with the input checks of circomlibjs (A and R8 on curve, S < subOrder).
R8 is compared in projective coordinates, without inversion; R8 = (0,1) is accepted only when the result
is the point at infinity. `pubKeyToW8` precomputes 8A in Weierstrass coordinates (3 Edwards doublings,
one modexp); `verifyRawW` then needs a single modexp.

## ADR-7: Edwards a = -1 version kept as a baseline (2026-11-17)

Context: on BabyJubJub, -168700 is a square and d' = -168696/168700 is not, so the a = -1 model has
complete extended formulas (7 mulmod per doubling, 6 per precomputed mixed addition).

Decision: ship `FCL_babyjubjub_edwards.sol` with the same window and batch normalization. It is 12.5 %
cheaper than the Weierstrass version and has no exceptional case. The Weierstrass version stays the
reference of the isomorphism approach (shared engine, generic curve precompile).
