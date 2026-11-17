// BabyJubJub EdDSA vectors produced by the circomlibjs signers.
// For each signature: hm (already hashed challenge), R8, S, A and the expected verdict.
// The raw oracle follows circomlibjs verifyPoseidon step by step, without the hash.
// usage: node gen_vectors.mjs <eddsa.json> <keys.json>
import { buildEddsa } from "circomlibjs";
import { Scalar } from "ffjavascript";
import crypto from "crypto";
import fs from "fs";

const eddsa = await buildEddsa();
const bj = eddsa.babyJub;
const F = bj.F;
const toB = (e) => F.toObject(e);
const l = bj.subOrder;
const r = Scalar.fromString("21888242871839275222246405745257275088548364400416034343698204186575808495617");

// oracle: raw circomlib equation S*B8 == R8 + (8*hm)*A, same input checks
function rawVerify(hm, R8, S, A) {
    if (!bj.inCurve(R8) || !bj.inCurve(A)) return false;
    if (Scalar.geq(S, l)) return false;
    const Pleft = bj.mulPointEscalar(bj.Base8, S);
    const Pright = bj.addPoint(R8, bj.mulPointEscalar(A, Scalar.mul(hm, 8)));
    return F.eq(Pleft[0], Pright[0]) && F.eq(Pleft[1], Pright[1]);
}

const pt = (x, y) => [F.e(x), F.e(y)];
const rndField = () => Scalar.mod(Scalar.fromRprBE(crypto.randomBytes(32), 0, 32), r);
const rows = [];
function push(tag, hm, R8, S, A, expect) {
    const ok = rawVerify(hm, R8, S, A);
    if (expect !== undefined && ok !== expect) throw new Error("oracle " + tag);
    rows.push({ tag, hm, r8x: toB(R8[0]), r8y: toB(R8[1]), s: S, ax: toB(A[0]), ay: toB(A[1]), ok });
}

const signers = {
    poseidon: (prv, m) => [eddsa.signPoseidon(prv, m), (s, A) => eddsa.poseidon([s.R8[0], s.R8[1], A[0], A[1], m]), (s, A) => eddsa.verifyPoseidon(m, s, A)],
    mimc7: (prv, m) => [eddsa.signMiMC(prv, m), (s, A) => eddsa.mimc7.multiHash([s.R8[0], s.R8[1], A[0], A[1], m]), (s, A) => eddsa.verifyMiMC(m, s, A)],
    mimcsponge: (prv, m) => [eddsa.signMiMCSponge(prv, m), (s, A) => eddsa.mimcSponge.multiHash([s.R8[0], s.R8[1], A[0], A[1], m]), (s, A) => eddsa.verifyMiMCSponge(m, s, A)],
};

const keys = [];
for (const [name, mk] of Object.entries(signers)) {
    const N = name === "poseidon" ? 24 : 6;
    for (let i = 0; i < N; i++) {
        const prv = crypto.randomBytes(32);
        const A = eddsa.prv2pub(prv);
        const m = F.e(rndField());
        const [sig, H, ver] = mk(prv, m);
        const hm = toB(H(sig, A));
        if (!ver(sig, A)) throw new Error("circomlibjs verify " + name);
        push(name, hm, sig.R8, sig.S, A, true);
        keys.push(A);
        if (i % 3 === 0) {
            push(name + ":S+1", hm, sig.R8, Scalar.mod(Scalar.add(sig.S, 1), l), A, false);
            push(name + ":hm+1", Scalar.add(hm, 1), sig.R8, sig.S, A, false);
            push(name + ":S+l", hm, sig.R8, Scalar.add(sig.S, l), A, false);
            push(name + ":R8 neg", hm, [F.neg(sig.R8[0]), sig.R8[1]], sig.S, A, false);
            push(name + ":A neg", hm, sig.R8, sig.S, [F.neg(A[0]), A[1]], false);
        }
    }
}
// Pedersen : hm = pedersen(R8p || Ap || msg), entier 256 bits petit-boutiste
for (let i = 0; i < 4; i++) {
    const prv = crypto.randomBytes(32);
    const A = eddsa.prv2pub(prv);
    const msg = crypto.randomBytes(32);
    const sig = eddsa.signPedersen(prv, msg);
    const buf = new Uint8Array(64 + msg.length);
    buf.set(bj.packPoint(sig.R8), 0);
    buf.set(bj.packPoint(A), 32);
    buf.set(msg, 64);
    const hm = Scalar.fromRprLE(eddsa.pedersenHash.hash(buf), 0, 32);
    if (!eddsa.verifyPedersen(msg, sig, A)) throw new Error("pedersen");
    push("pedersen", hm, sig.R8, sig.S, A, true);
}

// edge cases: small order points and torsion components
const O = pt(0, 1), T2 = pt(0, r - 1n);
// points d'ordre 4 : x^2 = 1/a  (y = 0)
const inva = F.inv(F.e(168700));
const x4 = F.sqrt(inva);
const T4 = [x4, F.e(0)];
const k0 = crypto.randomBytes(32);
const A0 = eddsa.prv2pub(k0);
const m0 = F.e(7);
const s0 = eddsa.signPoseidon(k0, m0);
const h0 = toB(eddsa.poseidon([s0.R8[0], s0.R8[1], A0[0], A0[1], m0]));
push("R8+T2", h0, bj.addPoint(s0.R8, T2), s0.S, A0, false);
push("A+T4 (8A unchanged)", h0, s0.R8, s0.S, bj.addPoint(A0, T4));
push("A+T2", h0, s0.R8, s0.S, bj.addPoint(A0, T2));
push("A=O, S=0, R8=O", 123n, O, 0n, O, true);
push("A=O, S=0, R8=T2", 123n, T2, 0n, O, false);
push("A=T4, S=0, R8=O", 5n, O, 0n, T4, true);
push("A=T2, S=1, R8=B8", 9n, bj.Base8, 1n, T2, true);
push("hm=0, R8=S.B8", 0n, bj.mulPointEscalar(bj.Base8, 77n), 77n, A0, true);
push("hm=l", l, bj.mulPointEscalar(bj.Base8, 77n), 77n, A0, true);
push("S=0, hm=0, R8=O", 0n, O, 0n, A0, true);
push("R8 off curve", h0, [s0.R8[0], F.add(s0.R8[1], F.e(1))], s0.S, A0, false);
push("A off curve", h0, s0.R8, s0.S, [A0[0], F.add(A0[1], F.e(1))], false);
// degenerate keys: 8A = c.B8 for small c (engine fallback path)
for (const c of [1n, 2n, 3n, l - 1n, l - 2n, (l + 1n) / 2n]) {
    const inv8 = (() => { let x = 1n, b = 8n, e = l - 2n; while (e) { if (e & 1n) x = x * b % l; b = b * b % l; e >>= 1n; } return x; })();
    const A = bj.mulPointEscalar(bj.Base8, c * inv8 % l);
    const S = BigInt("0x" + crypto.randomBytes(31).toString("hex")) % l;
    const hm = rndField();
    // R8 = S.B8 - 8 hm A for a valid signature
    const P = bj.mulPointEscalar(bj.Base8, S);
    const Q = bj.mulPointEscalar(A, Scalar.mul(hm, 8));
    const R8 = bj.addPoint(P, [F.neg(Q[0]), Q[1]]);
    push("8A=" + c + ".B8", hm, R8, S, A, true);
    push("8A=" + c + ".B8 wrong", hm + 1n, R8, S, A, false);
}

const out = {};
for (const k of ["hm", "r8x", "r8y", "s", "ax", "ay"]) out[k] = rows.map((x) => "0x" + BigInt(x[k]).toString(16));
out.ok = rows.map((x) => (x.ok ? 1 : 0));
out.tag = rows.map((x) => x.tag);
fs.writeFileSync(process.argv[2], JSON.stringify(out));
fs.writeFileSync(process.argv[3], JSON.stringify(keys.map((A) => [toB(A[0]).toString(), toB(A[1]).toString()])));
console.log(rows.length, "vectors,", rows.filter((x) => x.ok).length, "valid");
