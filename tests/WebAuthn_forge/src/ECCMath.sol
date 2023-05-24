pragma solidity ^0.8.13;

/**
 * @title ECCMath
 *
 * Functions for working with integers, curve-points, etc.
 *
 * @author Andreas Olofsson (androlo1980@gmail.com)
 */
library ECCMath {
    /// @dev Modular inverse of a (mod p) using euclid.
    /// 'a' and 'p' must be co-prime.
    /// @param a The number.
    /// @param p The mmodulus.
    /// @return x such that ax = 1 (mod p)
    function invmod(uint a, uint p) internal returns (uint) {
       unchecked{
        if (a == 0 || a == p || p == 0)
            revert();
        if (a > p)
            a = a % p;
        int t1;
        int t2 = 1;
        uint r1 = p;
        uint r2 = a;
        uint q;
        while (r2 != 0) {
            q = r1 / r2;
            (t1, t2, r1, r2) = (t2, t1 - int(q) * t2, r2, r1 - q * r2);
        }
        if (t1 < 0)
            return (p - uint(-t1));
        return uint(t1);
        }
    }

    /// @dev Modular exponentiation, b^e % m
    /// Basically the same as can be found here:
    /// https://github.com/ethereum/serpent/blob/develop/examples/ecc/modexp.se
    /// @param b The base.
    /// @param e The exponent.
    /// @param m The modulus.
    /// @return r such that x = b**e (mod m)
    function expmod(uint b, uint e, uint m) internal  returns (uint r) {
        if (b == 0)
            return 0;
        if (e == 0)
            return 1;
        if (m == 0)
            revert();
        r = 1;
        uint bit = 2 ** 255;
        assembly {
            for{} gt(bit,0) {}{ 
                r := mulmod(mulmod(r, r, m), exp(b, iszero(iszero(and(e, bit)))), m)
                r := mulmod(mulmod(r, r, m), exp(b, iszero(iszero(and(e, div(bit, 2))))), m)
                r := mulmod(mulmod(r, r, m), exp(b, iszero(iszero(and(e, div(bit, 4))))), m)
                r := mulmod(mulmod(r, r, m), exp(b, iszero(iszero(and(e, div(bit, 8))))), m)
                bit := div(bit, 16)
               
        }
        }
    }
  
    ///  @dev Converts a point (Px, Py, Pz) expressed in Jacobian coordinates to affine coordinates
    /// Mutates P.
    /// @param P The point.
    /// @param zInv The modular inverse of 'Pz'.
    /// @param z2Inv The square of zInv
    /// @param prime The prime modulus.
    function toZ1(uint[3] memory P, uint zInv, uint z2Inv, uint prime) internal  {
        P[0] = mulmod(P[0], z2Inv, prime);
        P[1] = mulmod(P[1], mulmod(zInv, z2Inv, prime), prime);
        P[2] = 1;
    }

    /// @dev See _toZ1(uint[3], uint, uint).
    /// Warning: Computes a modular inverse.
    /// @param PJ The point.
    /// @param prime The prime modulus.
    function toZ1(uint[3] memory PJ, uint prime) internal  {
        uint zInv = invmod(PJ[2], prime);
        uint zInv2 = mulmod(zInv, zInv, prime);
        PJ[0] = mulmod(PJ[0], zInv2, prime);
        PJ[1] = mulmod(PJ[1], mulmod(zInv, zInv2, prime), prime);
        PJ[2] = 1;
    }

}
