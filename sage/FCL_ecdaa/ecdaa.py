##*************************************************************************************/
##/* Copyright (C) 2022 - Renaud Dubois - This file is part of cairo_musig2 project	 */
##/* License: This software is licensed under a dual BSD and GPL v2 license. 	 */
##/* See LICENSE file at the root folder of the project.				 */
##/* FILE: ecdaa.py							             	  */
##/* 											  */
##/* 											  */
##/* DESCRIPTION: ECDAA algorithm*/
##/* This is a high level simulation for validation purpose				  */
##/* 
##https://fidoalliance.org/specs/fido-v2.0-id-20180227/fido-ecdaa-algorithm-v2.0-id-20180227.html#ecdaa-join-algorithm           				  */
## The only difference lies in the use of cofactor cleaning for non BN curves
##/* note that some constant aggregating values could be precomputed			  */
##**************************************************************************************/
from sage.all_cmdline import *   # import sage library

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.functional import cyclotomic_polynomial
from sage.rings.finite_rings.finite_field_constructor import FiniteField, GF
from sage.schemes.elliptic_curves.constructor import EllipticCurve

# this is much much faster with this statement:
# proof.arithmetic(False)
from sage.structure.proof.all import arithmetic

from sage.crypto.util import bin_to_ascii

from hashlib import *
from sha3 import keccak_256    

#uncomment the curve to be used
#BLS12_381 : future of Ethereum precompiled contracts

#from common.arithmetic.curves.bls12_381 import long_pairing as _e #dunnow why short name fail
#from common.arithmetic.curves.bls12_381 import *

#ALTBN128 aka BN254 : current Ethereum precompiled contracts
from FCL_pairings.arithmetic.curves.atlbn128 import long_pairing as _e
from FCL_pairings.arithmetic.curves.atlbn128 import *



##################################################################
################################# 3.1 Object Encodings
##################################################################


def to_bytes(n, length, endianess='big'):
    h = '%x' % n
    s = ('0'*(len(h) % 2) + h).zfill(length*2).decode('hex')
    return s if endianess == 'big' else s[::-1]

#Convert an int to the input of sha function
#inVal: integer reprensenting the number to convert
#size: size in bytes of the output 
def BigNumberToB(inVal, size):
 bin_a=(bin(inVal)[2:]);
 bin_a=bin_a.zfill(size*8);

 return bin_to_ascii(bin_a).encode();
 
 

#3.1.2 Encoding EcPoint values as byte strings
def ECPointToB():
 return();
 
#3.1.3 Encoding EcPoint2 values as byte strings
def ECPoint2ToB():
 return();
 
    
#returns a nonce in ZP (beware that p is the order of the curve, noted r elsewhere)
def get_ZPnonce():
 return randint(0,r-1) 


##################################################################
################################# 3.2 Global ECDAA System Parameters
##################################################################
# 3.2.1 Groups G1, G2 and GT of sufficiently large prime order
G1=E1;
G2=E2;
# 3.2.2 Two generators P1 and P2 such that G1=<P1> and G2=<P2>
P1=Gen1;
P2=Gen2;
order_S1=len(bin(r)[2:]);#bitsize of the curve order
order_S8=(order_S1//8)+( (order_S1%8)!=0); #bytesize of curve order, (left padding to 0)
modulus_S1=len(bin(p)[2:]);
modulus_S8=(modulus_S1//8)+( (modulus_S1%8)!=0);



Fq=GF(p);
Fp=GF(r); #yes i hate those notations

# 3.2.3 A bilinear pairing e:G1xG2->GT, we use "ate pairing" over BLS12_381 or BN254 (see imports)
# note: implemented as _e(P,Q)
# 3.2.3 A bilinear pairing e:G1xG2->GT, we use "ate pairing" over BLS12_381 or BN254 (see imports)


##################################################################
################################# HASHING
#/*https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-13#name-bls12-381-g1 specifies use of SHA256 for H*/
##################################################################

#################################
##Hash on Zp
# 3.2.4 hash function hash with H:{0,1}* -> Zp, to reduce bias, SHA512 is recommended, thus EVM compatibility requires sha256
#################################
# the core hash function used is defined here
def H_Zp_init():
 ctx = keccak_256(); 
 return ctx;


#update state with an input of given bytesize 
def H_Zp_update(ctx, input,bytesize):
 ctx.update(BigNumberToB(input,bytesize ));
 return ctx;

 
#update state with an element of Zp (order)
def H_Zp_update_Zp(ctx, input):
 ctx.update(BigNumberToB(input,order_S8 ));
 return ctx;


#update state with an element of Fq (characteristic)
def H_Zp_update_Fq(ctx, input):
 ctx.update(BigNumberToB(input,modulus_S8 ));
 return ctx;


#update state with an element of G1
def H_Zp_update_G1(ctx, input):
 H_Zp_update_Fq(ctx, input[0] );
 H_Zp_update_Fq(ctx, input[1] );
 
 return ctx;
	
#update state with an element of G2
def H_Zp_update_G2(ctx, input):
 H_Zp_update_Fq(ctx, input[0].polynomial().list()[1] );
 H_Zp_update_Fq(ctx, input[0].polynomial().list()[0] );
 H_Zp_update_Fq(ctx, input[1].polynomial().list()[1] );
 H_Zp_update_Fq(ctx, input[1].polynomial().list()[0] );
 
 return ctx;

def H_Zp_final(ctx):
  return int('0x'+ctx.hexdigest(),16)%r;
 
def H_Zp(int_a):
 ctx=H_Zp_init();
 ctx=H_Zp_update_Zp(ctx, int_a);
 return H_Zp_final(ctx)

def H_Zp_long(int_a, size_S8):
  ctx=H_Zp_init();
  ctx=H_Zp_update(ctx, int_a, size_S8);	
  return int('0x'+ctx.hexdigest(),16)%r;
  
#################################
#Hash_onG1
#################################
#Hash on curve
def Hash_onG1(m):
 flag=false;
 i=0;
  
 while (flag==false):
  #H((i,4)|m)
  x=H_Zp(i+(m*2^32));
  z=Fq(x)**3+Fq(b);
  if(is_square(z)==true):
    flag=true
  else:
    i=i+1
    if(i==232):
      return false
 y=int(sqrt(z))
 minus_y=p-y
 Gy=min(y, minus_y)
 Pcleaned=c*E1([x, Gy]);
 
 return Pcleaned
 

#################################
#Hash_pre
#################################
def HG1_pre(m):

 flag=false;
 i=0;
  
 while (flag==false):
  #H((i,4)|m)
  sc=i+(m*2^32)
  x=H_Zp(sc);
  z=Fq(x)**3+Fq(b);
  if(is_square(z)==true):
    flag=true
  else:
    i=i+1
    if(i==232):
      return false
 y=int(sqrt(z))
 minus_y=p-y
 Gy=min(y, minus_y)
 Pcleaned=c*E1([x, Gy]); #for BLS curves, cofactor cleaning is required
  
 return sc, int(Pcleaned[1])
 
#################################
#Deriv B
#################################
def Deriv_B(sc,yc):
 m=sc//2^32;
 x=H_Zp(sc);
 z=Fq(x)**3+Fq(b);
 if(is_square(z)==false):
  return E1([0,1,0]) ## error, the value is not square
 y=int(sqrt(z))
 minus_y=p-y;
 Gy=min(y, minus_y)

 Pcleaned=c*E1([x, Gy]);
 if(yc!=Pcleaned[1]):
  return [0,1,0] ## error, the value is not expected value
 B=Pcleaned;
 
 return B

##################################################################
################################# SETUP
##################################################################
#3.3 Issuer Specific parameters
#Description: this function generates the secret and public parameters
# of the Issuer. All privacy security depends on the sk generated.
#isk_x, isk_y, r_x, r_y, Ux, Uy, c,

#Generate Issuer secret parameter
def SetUp_Priv():
#1,2. Randomly generated ECDAA Issuer private key isk=(x,y) with 
#[x,y=RAND(p)]
 isk_x=get_ZPnonce();
 isk_y=get_ZPnonce();
#[rx=RAND(p)]
 r_x=get_ZPnonce();
#[ry=RAND(p)]
 r_y=get_ZPnonce();
 return isk_x,isk_y, r_x, r_y;

#Derivate Issuer public parameters
def SetUp_DerivPub(isk_x, isk_y, r_x, r_y):
#Issuer public Key
 X=isk_x*P2;
 Y=isk_y*P2;
#3. Ux=rx.P2 
 U_x=r_x*P2;
#4. Uy=rx.P2 
 U_y=r_y*P2;
 ctx_c=H_Zp_init();
#5. c=H(Ux|Uy|P2|X|Y)
 ctx_c=H_Zp_update_G2(ctx_c, U_x);#Ux
 ctx_c=H_Zp_update_G2(ctx_c, U_y);#Uy
 ctx_c=H_Zp_update_G2(ctx_c, P2);#P2 
 ctx_c=H_Zp_update_G2(ctx_c, X);#X
 ctx_c=H_Zp_update_G2(ctx_c, Y);#Y
 i_c=H_Zp_final(ctx_c);
#6. sx=r_x+c.x mod p
 s_x=(r_x+i_c*isk_x)%r 
#7. sy=r_y+c.y mod p
 s_y=(r_y+i_c*isk_y)%r 
    
 return(X,Y,i_c,s_x, s_y)

#H(sx.P2+c*X|syP2-cY|P2|X|Y)== c ?
def CheckSetup(X,Y,i_c,s_x, s_y):
 ctx_c=H_Zp_init();
 chunk1=s_x*P2-i_c*X;#sx.P2+c*X
 
 ctx_c=H_Zp_update_G2(ctx_c, chunk1)
 chunk2=s_y*P2-i_c*Y;#sy.P2-c*Y
 ctx_c=H_Zp_update_G2(ctx_c, chunk2)
 ctx_c=H_Zp_update_G2(ctx_c, P2);#P2 
 ctx_c=H_Zp_update_G2(ctx_c, X);#X
 ctx_c=H_Zp_update_G2(ctx_c, Y);#Y
 i_c2=H_Zp_final(ctx_c);

 return(i_c==i_c2);



##################################################################
################################# 3.4 JOIN
##################################################################
# The join algorithm is split in 4 steps:
# Step 2-3: Issuer generates credential B
# Step 4-8: Authenticator generates secret key and Proof of Knowledge
# Step 9-12: Issuers issues credential A,B,C,D
# Step 13+: ASM/Authenticator checks validity of received credentials


def Issuer_Join_Generate_B():
#2.The ECDAA Issuer chooses a nonce BigInteger m=RAND(p)
 m=get_ZPnonce();
#3.The ECDAA Issuer computes the B value of the credential as B
 sc, yc = HG1_pre(m);
 return sc, yc;


def Authenticator_Join_GenPriv( sc, yc):
 m=sc//2^32;
#4. Authenticator chooses and stores ECDAA private key sk=Rand(p) 
 sk=get_ZPnonce();
#5.The authenticator re-computes B = (H(sc),yc)
 B=Deriv_B(sc,yc);
#6. The authenticator computes its public key ECPoint Q=B.sk
 
 Q=sk*B;
#7.The authenticator proves knowledge of sk as follows
#7.1.BigInteger r1​​ =RAND(p)
 r1=get_ZPnonce();
#7.2 U1=B.r1
 U1=B*r1;
#7.3.​​BigInteger c​2​​ =H(U​1​​ ∣P​​​ ∣Q∣m)
 ctx=H_Zp_init();
 
 H_Zp_update_G1(ctx, U1);
 H_Zp_update_G1(ctx, P1);
 H_Zp_update_G1(ctx, Q);
 H_Zp_update_Zp(ctx, m);
 i_c2=H_Zp_final(ctx);
 #7.4.BigInteger n=RAND(p)
 n=get_ZPnonce();
#7.5. BigInteger ​​ c1=H(n∣c​2​​ )
 ctx=H_Zp_init();
 H_Zp_update_Zp(ctx, n);
 H_Zp_update_Zp(ctx, i_c2);
 i_c1=H_Zp_final(ctx);
#7.6. BigInteger s​1​​ =r​1​​ +c​1​​ ⋅sk
 s1=r1+i_c1*sk;
#8 Authenticator sends Q, c1,s1,n to Issuer
 return sk,Q, i_c1, s1, n; 

#10. issuer verifies that Q in G1 and H(n|H(s1B-c1Q|P1|Q|m): check proof of possesion of private key 
def Check_Proof(Q, B, m, c1, s1, n):
 U1=s1*B-c1*Q;
 
 ctx=H_Zp_init();
 H_Zp_update_G1(ctx, U1);
 H_Zp_update_G1(ctx, P1);
 H_Zp_update_G1(ctx, Q);
 H_Zp_update_Zp(ctx, m);
 i_c2=H_Zp_final(ctx);
 
 ctx=H_Zp_init();
 H_Zp_update_Zp(ctx, n);
 H_Zp_update_Zp(ctx, i_c2);
 i_c1=H_Zp_final(ctx);
 
 return (i_c1==c1)

#input (sk_x, sk_y): issuer secret key
def Issuer_Gen_Credentials(sk_x, sk_y, m, B, Q, c1, s1, n):
#11 The ECDAA Issuer creates credential (A,B,C,D) as follows
#11.1 A=B^(1/y)
 A=int(Fp(sk_y)**(-1))*B;
#11.2 B, la réponse B B=B
#11.3 C=(A*Q)^x;
 C=sk_x*(A+Q);
 return A,B,C,Q;

# 12, 13, 14, Authenticator checks parameters validity
def Check_Credentials(A,B,C,D, X, Y):
 flag=(_e(A,Y)==_e(B,P2)); #checks that e(A,Y)==e(B,P2) ?
 flag=flag and ( _e(C,P2)==_e(A+D, X));# checks that e(C,P2)==e(A+D, X) ?
 return flag;
 

##################################################################
################################# 3.5 SIGN
################################################################## 
#Input: 
#sk=user secret key
#A,B,C,D: credentials
#Data: some additional data of bytesize Data_s8
#h_KRD : hash of KRD (message) of bytesize Data_s8
def ECDAA_Sign(sk, A, B, C, D, Data, Data_s8, h_KRD, h_s8):
 #2. BigNumber l=rand(p)
 l=get_ZPnonce();
 #3. ECPoint R=A^l
 R=l*A;
 #4. EcPoint S=B^l
 S=l*B;
 #5. EcPoint T=C^l
 T=l*C;
 #6. EcPoint W=D^l
 W=l*D;
 #7. Big integer r = rand(p)
 rr=get_ZPnonce();
 #8. Ecpoint U=S^r
 U=rr*S;

 #9. BigInteger c2=H(U|S|W|AdditionalData|h_KRD) 
 ctx=H_Zp_init();
 H_Zp_update_G1(ctx, U);
 H_Zp_update_G1(ctx, S);
 H_Zp_update_G1(ctx, W);
 H_Zp_update(ctx, Data, Data_s8);
 H_Zp_update(ctx, h_KRD, h_s8);
 i_c2=H_Zp_final(ctx);
 
 #10 BigInteger n=rand(p)
 n= get_ZPnonce();
 #11 c=H(n|c2)
 ctx=H_Zp_init();
 H_Zp_update_Zp(ctx, n);
 H_Zp_update_Zp(ctx, i_c2);
 i_c=H_Zp_final(ctx);
 #12 BigInteger s=r+c.sk mod order
 s=(rr+i_c*sk)%r;
 
 return i_c,s,R,S,T,W,n


##################################################################
################################# 3.6 Verify
################################################################## 
#Input:
#X,Y: issuer public key
# Data, Data_s8, h_KRD, h_s8: additional data and message
#c,s,R,S,T,W,n : signature
def ECDAA_Verify(X,Y, Data, Data_s8, h_KRD, h_s8,  i_c,s,R,S,T,W,n):
 if R==E1([0,1,0]): 
  return false #Check R is not infinity point (1_G1)
 if S==E1([0,1,0]): 
  return false #Check S is not infinity point (1_G1)
 
 #3. H(n|H(S^s.W^-c|S|W|Data|H(KRD))) == c ? fail if not equal
 U=s*S-i_c*W;
 U2=s*S+(r-i_c)*W;

 #from here same as signing process hashing from U, S, W, Data, HKRD
 #9. BigInteger c2=H(U|S|W|AdditionalData|h_KRD) 
 ctx=H_Zp_init();
 H_Zp_update_G1(ctx, U);
 H_Zp_update_G1(ctx, S);
 H_Zp_update_G1(ctx, W);
 H_Zp_update(ctx, Data, Data_s8);
 H_Zp_update(ctx, h_KRD, h_s8);
 i_c2=H_Zp_final(ctx);
 
 #11 c=H(n|c2)
 ctx=H_Zp_init();
 H_Zp_update_Zp(ctx, n);
 H_Zp_update_Zp(ctx, i_c2);
 i_cprime=H_Zp_final(ctx);
 
 flag=true;
 if(i_cprime!=i_c):
  flag=false;
 #4 e(R,Y)!=e(S,P2), fail if not equal
 if(_e(R,Y)!=_e(S,P2)):
   flag=false;
 
 #5 e(T,P2)==e(R.W, X) ? fail if not equal  
 if(_e(T,P2)!=_e(R+W,X)):
   flag=false;

 #6 for all rogues, if W=S^sk, fail
 # It is assumed that rogue checking is done prior to the call
 
 return flag;

# Revocation check
#Description: This function parse a list of compromised key to reject
#signatures issued by revoked users
#Input: 
# Rogs: Rogue list of compromised sk
# W,S : W and S parts of signature
def ECDAA_KillRogue(W,S, Rogs):
 #6 for all rogues, if W=S^sk, fail
 for i in Rogs:
   if(W==Rogs[i]*S):
     return false;
   
   return true;
















 
 
 
