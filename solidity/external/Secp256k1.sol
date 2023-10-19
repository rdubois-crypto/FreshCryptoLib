pragma solidity >=0.8.19 <0.9.0;

import "./ECCMath.sol";

/**
 * @title Secp256k1
 *
 * secp256k1 implementation.
 *
 * The library implements 'Curve' and 'codec/ECCConversion', but since it's a library
 * it does not actually extend the contracts. This is a Solidity thing and will be
 * dealt with later.
 *
 * @author Andreas Olofsson (androlo1980@gmail.com)
 */
library Secp256k1 {
    // TODO separate curve from crypto primitives?

    // Field size
    uint256 constant pp = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F;

    // Base point (generator) G
    uint256 constant Gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798;
    uint256 constant Gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8;

    // Order of G
    uint256 constant nn = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141;

    // Cofactor
    // uint constant hh = 1;

    // Maximum value of s
    uint256 constant lowSmax = 0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5D576E7357A4501DDFE92F46681B20A0;

    // For later
    // uint constant lambda = "0x5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72";
    // uint constant beta = "0x7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee";

    /// @dev See Curve.onCurve
    function onCurve(uint256[2] memory P) internal pure returns (bool) {
        uint256 p = pp;
        if (0 == P[0] || P[0] == p || 0 == P[1] || P[1] == p) {
            return false;
        }
        uint256 LHS = mulmod(P[1], P[1], p);
        uint256 RHS = addmod(mulmod(mulmod(P[0], P[0], p), P[0], p), 7, p);
        return LHS == RHS;
    }

    /// @dev See Curve.isPubKey
    function isPubKey(uint256[2] memory P) internal pure   returns(bool isPK) {
        isPK = onCurve(P);
    }

    /// @dev See Curve.validateSignature
    function validateSignature(bytes32 message, uint256[2] memory rs, uint256[2] memory Q) internal pure returns (bool) {
        uint256 n = nn;
        uint256 p = pp;
        if (rs[0] == 0 || rs[0] >= n || rs[1] == 0 || rs[1] > lowSmax) {
            return false;
        }
        if (!isPubKey(Q)) {
            return false;
        }

        uint256 sInv = ECCMath.invmod(rs[1], n);
        uint256[3] memory u1G = _mul(mulmod(uint256(message), sInv, n), [Gx, Gy]);
        uint256[3] memory u2Q = _mul(mulmod(rs[0], sInv, n), Q);
        uint256[3] memory P = _add(u1G, u2Q);

        if (P[2] == 0) {
            return false;
        }

        uint256 Px = ECCMath.invmod(P[2], p); // need Px/Pz^2
        Px = mulmod(P[0], mulmod(Px, Px, p), p);
        return Px % n == rs[0];
    }

    /// @dev See Curve.compress
    function compress(uint256[2] memory P) internal pure returns (uint8 yBit, uint256 x) {
        x = P[0];
        yBit = P[1] & 1 == 1 ? 1 : 0;
    }

    /// @dev See Curve.decompress
    function decompress(uint8 yBit, uint256 x) internal pure returns (uint256[2] memory P) {
        uint256 p = pp;
        uint256 y2 = addmod(mulmod(x, mulmod(x, x, p), p), 7, p);
        uint256 y_ = ECCMath.expmod(y2, (p + 1) / 4, p);
        uint256 cmp = yBit ^ y_ & 1;
        P[0] = x;
        P[1] = (cmp == 0) ? y_ : p - y_;
    }

    // Point addition, P + Q
    // inData: Px, Py, Pz, Qx, Qy, Qz
    // outData: Rx, Ry, Rz
    function _add(uint256[3] memory P, uint256[3] memory Q) public pure returns (uint256[3] memory R) {
        if (P[2] == 0) {
            return Q;
        }
        if (Q[2] == 0) {
            return P;
        }
        uint256 p = pp;
        uint256[4] memory zs; // Pz^2, Pz^3, Qz^2, Qz^3
        zs[0] = mulmod(P[2], P[2], p);
        zs[1] = mulmod(P[2], zs[0], p);
        zs[2] = mulmod(Q[2], Q[2], p);
        zs[3] = mulmod(Q[2], zs[2], p);
        uint256[4] memory us =
            [mulmod(P[0], zs[2], p), mulmod(P[1], zs[3], p), mulmod(Q[0], zs[0], p), mulmod(Q[1], zs[1], p)]; // Pu, Ps, Qu, Qs
        if (us[0] == us[2]) {
            if (us[1] != us[3]) {
                revert();
            } //not sure
            else {
                return _double(P);
            }
        }
        uint256 h = addmod(us[2], p - us[0], p);
        uint256 r = addmod(us[3], p - us[1], p);
        uint256 h2 = mulmod(h, h, p);
        uint256 h3 = mulmod(h2, h, p);
        uint256 Rx = addmod(mulmod(r, r, p), p - h3, p);
        Rx = addmod(Rx, p - mulmod(2, mulmod(us[0], h2, p), p), p);
        R[0] = Rx;
        R[1] = mulmod(r, addmod(mulmod(us[0], h2, p), p - Rx, p), p);
        R[1] = addmod(R[1], p - mulmod(us[1], h3, p), p);
        R[2] = mulmod(h, mulmod(P[2], Q[2], p), p);
    }

    // Point addition, P + Q. P Jacobian, Q affine.
    // inData: Px, Py, Pz, Qx, Qy
    // outData: Rx, Ry, Rz
    function _addMixed(uint256[3] memory P, uint256[2] memory Q) internal pure returns (uint256[3] memory R) {
        if (P[2] == 0) {
            return [Q[0], Q[1], 1];
        }
        if (Q[1] == 0) {
            return P;
        }
        uint256 p = pp;
        uint256[2] memory zs; // Pz^2, Pz^3, Qz^2, Qz^3
        zs[0] = mulmod(P[2], P[2], p);
        zs[1] = mulmod(P[2], zs[0], p);
        uint256[4] memory us = [P[0], P[1], mulmod(Q[0], zs[0], p), mulmod(Q[1], zs[1], p)]; // Pu, Ps, Qu, Qs
        if (us[0] == us[2]) {
            if (us[1] != us[3]) {
                P[0] = 0;
                P[1] = 0;
                P[2] = 0;
                return P;
            } else {
                _double(P);
                return P;
            }
        }
        uint256 h = addmod(us[2], p - us[0], p);
        uint256 r = addmod(us[3], p - us[1], p);
        uint256 h2 = mulmod(h, h, p);
        uint256 h3 = mulmod(h2, h, p);
        uint256 Rx = addmod(mulmod(r, r, p), p - h3, p);
        Rx = addmod(Rx, p - mulmod(2, mulmod(us[0], h2, p), p), p);
        R[0] = Rx;
        R[1] = mulmod(r, addmod(mulmod(us[0], h2, p), p - Rx, p), p);
        R[1] = addmod(R[1], p - mulmod(us[1], h3, p), p);
        R[2] = mulmod(h, P[2], p);
    }

    // Same as addMixed but params are different and mutates P.
    function _addMixedM(uint256[3] memory P, uint256[2] memory Q) internal pure{
        if (P[1] == 0) {
            P[0] = Q[0];
            P[1] = Q[1];
            P[2] = 1;
            return;
        }
        if (Q[1] == 0) {
            return;
        }
        uint256 p = pp;
        uint256[2] memory zs; // Pz^2, Pz^3, Qz^2, Qz^3
        zs[0] = mulmod(P[2], P[2], p);
        zs[1] = mulmod(P[2], zs[0], p);
        uint256[4] memory us = [P[0], P[1], mulmod(Q[0], zs[0], p), mulmod(Q[1], zs[1], p)]; // Pu, Ps, Qu, Qs
        if (us[0] == us[2]) {
            if (us[1] != us[3]) {
                P[0] = 0;
                P[1] = 0;
                P[2] = 0;
                return;
            } else {
                _doubleM(P);
                return;
            }
        }
        uint256 h = addmod(us[2], p - us[0], p);
        uint256 r = addmod(us[3], p - us[1], p);
        uint256 h2 = mulmod(h, h, p);
        uint256 h3 = mulmod(h2, h, p);
        uint256 Rx = addmod(mulmod(r, r, p), p - h3, p);
        Rx = addmod(Rx, p - mulmod(2, mulmod(us[0], h2, p), p), p);
        P[0] = Rx;
        P[1] = mulmod(r, addmod(mulmod(us[0], h2, p), p - Rx, p), p);
        P[1] = addmod(P[1], p - mulmod(us[1], h3, p), p);
        P[2] = mulmod(h, P[2], p);
    }

    // Point doubling, 2*P
    // Params: Px, Py, Pz
    // Not concerned about the 1 extra mulmod.
    function _double(uint256[3] memory P) internal pure returns (uint256[3] memory Q) {
        uint256 p = pp;
        if (P[2] == 0) {
            return P;
        }
        uint256 Px = P[0];
        uint256 Py = P[1];
        uint256 Py2 = mulmod(Py, Py, p);
        uint256 s = mulmod(4, mulmod(Px, Py2, p), p);
        uint256 m = mulmod(3, mulmod(Px, Px, p), p);
        uint256 Qx = addmod(mulmod(m, m, p), p - addmod(s, s, p), p);
        Q[0] = Qx;
        Q[1] = addmod(mulmod(m, addmod(s, p - Qx, p), p), p - mulmod(8, mulmod(Py2, Py2, p), p), p);
        Q[2] = mulmod(2, mulmod(Py, P[2], p), p);
    }

    // Same as double but mutates P and is internal only.
    function _doubleM(uint256[3] memory P) internal pure {
        uint256 p = pp;
        if (P[2] == 0) {
            return;
        }
        uint256 Px = P[0];
        uint256 Py = P[1];
        uint256 Py2 = mulmod(Py, Py, p);
        uint256 s = mulmod(4, mulmod(Px, Py2, p), p);
        uint256 m = mulmod(3, mulmod(Px, Px, p), p);
        uint256 PxTemp = addmod(mulmod(m, m, p), p - addmod(s, s, p), p);
        P[0] = PxTemp;
        P[1] = addmod(mulmod(m, addmod(s, p - PxTemp, p), p), p - mulmod(8, mulmod(Py2, Py2, p), p), p);
        P[2] = mulmod(2, mulmod(Py, P[2], p), p);
    }

    // Multiplication dP. P affine, wNAF: w=5
    // Params: d, Px, Py
    // Output: Jacobian Q
    function _mul(uint256 d, uint256[2] memory P) public pure returns (uint256[3] memory Q) {
        uint256 p = pp;
        if (
            d == 0 // TODO
        ) {
            return [P[0] * 0, P[0] * 0, P[0] * 0];
        }
        uint256 dwPtr; // points to array of NAF coefficients.
        uint256 i;

        // wNAF
        assembly {
            let dm := 0
            dwPtr := mload(0x40)
            mstore(0x40, add(dwPtr, 512)) // Should lower this.
            for {} gt(d, 0) {} {
                if iszero(and(d, 1)) {
                    d := div(d, 2)
                    i := add(i, 1)
                    continue
                }

                dm := mod(d, 32)
                mstore8(add(dwPtr, i), dm) // Don't store as signed - convert when reading.
                d := add(sub(d, dm), mul(gt(dm, 16), 32))
            }
        }

        // Pre calculation
        uint256[3][8] memory PREC; // P, 3P, 5P, 7P, 9P, 11P, 13P, 15P
        PREC[0] = [P[0], P[1], 1];
        uint256[3] memory X = _double(PREC[0]);
        PREC[1] = _addMixed(X, P);
        PREC[2] = _add(X, PREC[1]);
        PREC[3] = _add(X, PREC[2]);
        PREC[4] = _add(X, PREC[3]);
        PREC[5] = _add(X, PREC[4]);
        PREC[6] = _add(X, PREC[5]);
        PREC[7] = _add(X, PREC[6]);

        uint256[16] memory INV;
        INV[0] = PREC[1][2]; // a1
        INV[1] = mulmod(PREC[2][2], INV[0], p); // a2
        INV[2] = mulmod(PREC[3][2], INV[1], p); // a3
        INV[3] = mulmod(PREC[4][2], INV[2], p); // a4
        INV[4] = mulmod(PREC[5][2], INV[3], p); // a5
        INV[5] = mulmod(PREC[6][2], INV[4], p); // a6
        INV[6] = mulmod(PREC[7][2], INV[5], p); // a7

        INV[7] = ECCMath.invmod(INV[6], p); // a7inv
        INV[8] = INV[7]; // aNinv (a7inv)

        INV[15] = mulmod(INV[5], INV[8], p); // z7inv
        for (uint256 k = 6; k >= 2; k--) {
            // z6inv to z2inv
            INV[8] = mulmod(PREC[k + 1][2], INV[8], p);
            INV[8 + k] = mulmod(INV[k - 2], INV[8], p);
        }
        INV[9] = mulmod(PREC[2][2], INV[8], p); // z1Inv
        for (uint256 k = 0; k < 7; k++) {
            ECCMath.toZ1(PREC[k + 1], INV[k + 9], mulmod(INV[k + 9], INV[k + 9], p), p);
        }

        // Mult loop
        while (i > 0) {
            uint256 dj;
            uint256 pIdx;
            i--;
            assembly {
                dj := byte(0, mload(add(dwPtr, i)))
            }
            _doubleM(Q);
            if (dj > 16) {
                pIdx = (31 - dj) / 2; // These are the "negative ones", so invert y.
                _addMixedM(Q, [PREC[pIdx][0], p - PREC[pIdx][1]]);
            } else if (dj > 0) {
                pIdx = (dj - 1) / 2;
                _addMixedM(Q, [PREC[pIdx][0], PREC[pIdx][1]]);
            }
        }
    }
}
