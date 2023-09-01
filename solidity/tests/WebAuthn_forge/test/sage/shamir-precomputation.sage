#//**********************************************************************************************/
#//  ___           _       ___               _         _    _ _
#// | __| _ ___ __| |_    / __|_ _ _  _ _ __| |_ ___  | |  (_) |__
#// | _| '_/ -_|_-< ' \  | (__| '_| || | '_ \  _/ _ \ | |__| | '_ \
#// |_||_| \___/__/_||_|  \___|_|  \_, | .__/\__\___/ |____|_|_.__/
#//                                |__/|_|
#// Copyright (C) 2022 - Renaud Dubois - This file is part of FCL (Fresh CryptoLib) project     */
#// License: This software is licensed under MIT License 	                                    */
#// See LICENSE file at the root folder of the project.				                            */
#// FILE: FCL_ecdsa_precompute.sage						                                        */
#// 											                                                */
#// 											                                                */
#// DESCRIPTION: precompute a 8 dimensional table for Shamir's trick from a public key
#//
#//**********************************************************************************************/

import os
import sys

# load he library needed to perform the precomputation
load("../../sage/FCL_common/FCL_elliptic.sage");

#//curve secp256r1, aka p256
#//curve prime field modulus
sec256p_p = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF;
#//short weierstrass first coefficient
sec256p_a =0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC;
#//short weierstrass second coefficient
sec256p_b =0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B;
#//generating point affine coordinates
sec256p_gx =0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296;
sec256p_gy =0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5;
#//curve order (number of points)
sec256p_n =0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551;

#//Init the curve
secp256r1, G = Init_Curve(sec256p_p, sec256p_a, sec256p_b, sec256p_gx, sec256p_gy, sec256p_n);

def precompute_pubkey(Q, Curve):
  # Initialize the precomputed table and the powers of 64*q
  Pow64_PQ=[ Q for i in range(0,16)];
  Prec=[ Curve(0) for i in range(0,256)];
  Pow64_PQ[0]=Curve([sec256p_gx, sec256p_gy]);
  Pow64_PQ[4]=Q;

   # Compute the powers of 64*q
  for j in [1..3]:
    Pow64_PQ[j]=2^64*Pow64_PQ[j-1];
    Pow64_PQ[j+4]=2^64*Pow64_PQ[j+3];

  # Compute the precomputed table
  Prec[0]=Curve(0);
  for i in range(1,256):
    Prec[i]=Curve(0);
    for j in [0..7]:
      if( (i&(1<<j))!=0):
        (Prec[i])=(Pow64_PQ[j]+ Prec[i]);

  # Return the precomputed table
  return Prec;

def print_setlength(X, n):
    # Convert X to a hexadecimal string and remove the '0x' prefix
    hex_str = hex(X)[2:]

    # Add leading zeros to the string to make it length n
    return hex_str.rjust(n, '0')

def get_concatenate_point(Prec):
  coords = ""

  # Concatenate the x and y coordinates of each point in Prec
  for i in [0..255]:
    px=print_setlength( Prec[i][0], 64);
    py=print_setlength( Prec[i][1], 64);
    coords=coords+px+py;

  # Join the coordinates into a single string
  return coords

if __name__ == '__main__':
  # Load the C0 and C1 environment variables
  C0 = int(os.environ['C0'])
  C1 = int(os.environ['C1'])

  # Compute the precomputed table
  q = secp256r1([C0, C1])
  Prec = precompute_pubkey(q, secp256r1)

  # Get the concatenated points and write them to stdout
  concatenated_points = get_concatenate_point(Prec)
  sys.stdout.write(concatenated_points)
