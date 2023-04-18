// SPDX-License-Identifier: UNLICENSED
pragma solidity ^0.8.13;

import "forge-std/Test.sol";
import "src/FCL_elliptic.sol";

contract Wrap_ecdsa{

function wrap_ecdsa_core( bytes32 message,
        uint[2] calldata rs,
        uint[2] calldata Q) public returns (bool)
 {
  return FCL_Elliptic_ZZ.ecdsa_verify(message,rs,Q);
 }
 
 constructor(){
 
 }
 
}

contract CounterTest is Test {

   //curve prime field modulus
    uint constant p = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF;
    //short weierstrass first coefficient
    uint constant a =
        0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC;
    //short weierstrass second coefficient    
    uint constant b =
        0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B;
    //generating point affine coordinates    
    uint constant gx =
        0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296;
    uint constant gy =
        0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5;
    //curve order (number of points)
    uint constant n =
        0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551;    
    /* -2 mod p constant, used to speed up inversion and doubling (avoid negation)*/
    uint constant minus_2 = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFD;
    /* -2 mod n constant, used to speed up inversion*/
    uint constant minus_2modn = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC63254F; 
       
    uint constant minus_1=      0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF;
 function setUp() public {
    
    uint256 i_u256_a=0x703b99d600000000000000000000000000000000000000000000000000000000;
    
    uint256 res=FCL_Elliptic_ZZ.FCL_nModInv(i_u256_a);
  //  console.log("inv:",mulmod(i_u256_a, res,FCL_Elliptic_ZZ.n ));
  
 
//     wx=vm.parseJsonUint(deployData, ".key[0]");
//    console.log("Read:",wx);

   
 }
 
 function uintToBytes(uint v) public returns (bytes32 ret) {
    if (v == 0) {
        ret = '0';
    }
    else {
        while (v > 0) {
            ret = bytes32(uint(ret) / (2 ** 8));
            ret |= bytes32(((v % 10) + 48) * 2 ** (8 * 31));
            v /= 10;
        }
    }
    return ret;
}

 function test_Fuzz_InVmodn(uint256 i_u256_a) public
 { 
  vm.assume(i_u256_a<FCL_Elliptic_ZZ.n );
  vm.assume(i_u256_a!=0 );
  
  uint256 res=FCL_Elliptic_ZZ.FCL_nModInv(i_u256_a);   	
  
  assertEq(mulmod(res, i_u256_a,FCL_Elliptic_ZZ.n),1);
 
 }

 function test_Fuzz_InVmodp(uint256 i_u256_a) public
 { 
  vm.assume(i_u256_a<FCL_Elliptic_ZZ.p );
  vm.assume(i_u256_a!=0 );
  
  uint256 res=FCL_Elliptic_ZZ.FCL_pModInv(i_u256_a);   	
  
  assertEq(mulmod(res, i_u256_a,FCL_Elliptic_ZZ.p),1);
 
 }
 //ecAff_isOnCurve
  
  function ReadField_th(string memory json, string memory field, uint256 num) public returns (uint256)
  {
   string memory snum=string(abi.encodePacked(uintToBytes(num)));
   console.log("snum:",snum);
 
   string memory access=string.concat(field, snum);	
   console.log("access:",access);
   
   return vm.parseJsonUint(json, access);
   return 0;
  }
  
  
    /**
     * @dev Computation of uG+vQ using Strauss-Shamir's trick, G basepoint, Q public key
     */
     function ecZZ_mulmuladd_S_asm(
        uint Q0, uint Q1,// Point G and Q stored in one memory for stack optimization
        uint scalar_u,
        uint scalar_v
    ) internal returns (uint X) {
     uint zz;
     uint zzz;
     uint Y;
     uint index=255;
     uint[6] memory T;
     uint H0;
     uint H1;   
     
     unchecked {
     
     if(scalar_u==0 && scalar_v==0) return 0;
     
     (H0,H1 )=FCL_Elliptic_ZZ.ecAff_add(gx,gy,Q0, Q1);//will not work if Q=P, obvious forbidden private key
   
   /*
     while( ( ((scalar_u>>index)&1)+2*((scalar_v>>index)&1) ) ==0){
      index=index-1; 
     }
     */
         
      assembly{
      
      for{  let T4:=add( shl(1, and(shr(index, scalar_v),1)), and(shr(index, scalar_u),1) )
      } eq(T4,0) {index := sub(index, 1)}
      {
          T4:=add( shl(1, and(shr(index, scalar_v),1)), and(shr(index, scalar_u),1) )
     
      }
      
       zz:=add( shl(1, and(shr(index, scalar_v),1)), and(shr(index, scalar_u),1) )
           
      if eq(zz,1) {
      	X:=gx
      	Y:=gy
      	}
      if eq(zz,2) {
       X:=Q0
      	Y:=Q1
      }
      if eq(zz,3) {
      	 X:=H0
      	 Y:= H1
      }
     
      index:=sub(index,1)
      zz:=1
      zzz:=1
      
      
      for {   } gt( minus_1, index) { index := sub(index, 1) } 
      {
      // inlined EcZZ_Dbl
      let T1:=mulmod(2, Y, p) //U = 2*Y1, y free
      let T2:=mulmod(T1,T1,p)  // V=U^2
      let T3:=mulmod(X, T2,p)// S = X1*V
      T1:=mulmod(T1, T2,p) // W=UV
      let T4:=mulmod(3, mulmod(addmod(X,sub(p,zz),p), addmod(X,zz,p),p) ,p) //M=3*(X1-ZZ1)*(X1+ZZ1)
      zzz:=mulmod(T1,zzz,p)//zzz3=W*zzz1
      zz:=mulmod(T2, zz, p) //zz3=V*ZZ1, V free
     
      X:=addmod(mulmod(T4,T4,p), mulmod(minus_2, T3,p),p) //X3=M^2-2S
      //T2:=mulmod(T4,addmod(T3, sub(p, X),p),p)//M(S-X3)
      T2:=mulmod(T4,addmod(X, sub(p, T3),p),p)//-M(S-X3)=M(X3-S)
      
      //Y:= addmod(T2, sub(p, mulmod(T1, Y ,p)),p  )//Y3= M(S-X3)-W*Y1
      Y:= addmod(mulmod(T1, Y ,p), T2,p  )//-Y3= W*Y1-M(S-X3), we replace Y by -Y to avoid a sub in ecAdd
      
     { 
      //value of dibit	
      T4:=add( shl(1, and(shr(index, scalar_v),1)), and(shr(index, scalar_u),1) )
      
      if iszero(T4){
       Y:=sub(p,Y)//restore the -Y inversion 
       continue
      }// if T4!=0
        
      if eq(T4,1) {
      	T1:=gx
      	T2:=gy
      	
      	}
      if eq(T4,2) {
        T1:=Q0
      	T2:=Q1
      }
      if eq(T4,3) {
      	 T1:=H0
      	 T2:= H1
      	 }
      	 	 
       // inlined EcZZ_AddN
      //T3:=sub(p, Y)
      //T3:=Y
      let y2:=addmod(mulmod(T2, zzz,p),Y,p)  
      T2:=addmod(mulmod(T1, zz,p),sub(p,X),p)  
      
      T4:=mulmod(T2, T2, p)//PP
      let TT1:=mulmod(T4,T2,p)//PPP, this one could be spared, but adding this register spare gas
      zz:=mulmod(zz,T4,p) 
      zzz:= mulmod(zzz,TT1,p) //zz3=V*ZZ1
      let TT2:=mulmod(X, T4, p)
      T4:=addmod(addmod(mulmod(y2,y2, p), sub(p,TT1),p ), mulmod(minus_2, TT2,p) ,p )
      Y:=addmod(mulmod(addmod(TT2, sub(p,T4),p), y2, p), mulmod(Y, TT1,p),p)
     
      X:=T4
       }
          
           }//end loop
        mstore(add(T, 0x60),zz)
      //(X,Y)=ecZZ_SetAff(X,Y,zz, zzz);
      //T[0] = inverseModp_Hard(T[0], p); //1/zzz, inline modular inversion using precompile:
     // Define length of base, exponent and modulus. 0x20 == 32 bytes
      mstore(T, 0x20)
      mstore(add(T, 0x20), 0x20)
      mstore(add(T, 0x40), 0x20)
      // Define variables base, exponent and modulus
      //mstore(add(pointer, 0x60), u)
      mstore(add(T, 0x80), minus_2)
      mstore(add(T, 0xa0), p)
               
      // Call the precompiled contract 0x05 = ModExp
      if iszero(call(not(0), 0x05, 0, T, 0xc0, T, 0x20)) {
            revert(0, 0)
      }
       
      //Y:=mulmod(Y,zzz,p)//Y/zzz
      //zz :=mulmod(zz, mload(T),p) //1/z
      //zz:= mulmod(zz,zz,p) //1/zz
      X:=mulmod(X,mload(T),p)//X/zz
      } //end assembly
     }//end unchecked
     
      return X;
    }
    
       
      //8 dimensions Shamir's trick, using precomputations stored in Shamir8,  stored as Bytecode of an external
      //contract at given address dataPointer
      //(thx to Lakhdar https://github.com/Kelvyne for EVM storage explanations and tricks)
      // the external tool to generate tables from public key is in the /sage directory
    function ecZZ_mulmuladd_S8_extcode(uint scalar_u, uint scalar_v, address dataPointer) 
    internal  returns(uint X/*, uint Y*/)
    {
    
      unchecked{ 
      uint zz; // third and  coordinates of the point
     
      uint[6] memory T;
      zz=256;//start index
      
      
      while(T[0]==0)
      {
      zz=zz-1;
      //tbd case of msb octobit is null
      T[0]=64*(128*((scalar_v>>zz)&1)+64*((scalar_v>>(zz-64))&1)+
           32*((scalar_v>>(zz-128))&1)+16*((scalar_v>>(zz-192))&1)+
               8*((scalar_u>>zz)&1)+4*((scalar_u>>(zz-64))&1)+2*((scalar_u>>(zz-128))&1)+((scalar_u>>(zz-192))&1));
      }
     assembly{
   
      extcodecopy(dataPointer, T, mload(T), 64)
      
      X:= mload(T)
      let Y:= mload(add(T,32))
      let zzz:=1
      zz:=1
     
      //loop over 1/4 of scalars thx to Shamir's trick over 8 points
      for { let index := 254 } gt(index, 191) { index := add(index, 191) } 
      { 
   	{
      let TT1:=mulmod(2, Y, p) //U = 2*Y1, y free
      let T2:=mulmod(TT1,TT1,p)  // V=U^2
      let T3:=mulmod(X, T2,p)// S = X1*V
      let T1:=mulmod(TT1, T2,p) // W=UV
      let T4:=mulmod(3, mulmod(addmod(X,sub(p,zz),p), addmod(X,zz,p),p) ,p) //M=3*(X1-ZZ1)*(X1+ZZ1)
      zzz:=mulmod(T1,zzz,p)//zzz3=W*zzz1
      zz:=mulmod(T2, zz, p) //zz3=V*ZZ1, V free
     
      X:=addmod(mulmod(T4,T4,p), mulmod(minus_2, T3,p),p) //X3=M^2-2S
      //T2:=mulmod(T4,addmod(T3, sub(p, X),p),p)//M(S-X3)
      let T5:=mulmod(T4,addmod(X, sub(p, T3),p),p)//-M(S-X3)=M(X3-S)
      
      //Y:= addmod(T2, sub(p, mulmod(T1, Y ,p)),p  )//Y3= M(S-X3)-W*Y1
      Y:= addmod(mulmod(T1, Y ,p), T5,p  )//-Y3= W*Y1-M(S-X3), we replace Y by -Y to avoid a sub in ecAdd
       
      /* compute element to access in precomputed table */
      }
      {
      let T4:= add( shl(13, and(shr(index, scalar_v),1)), shl(9, and(shr(index, scalar_u),1)) )
      let index2:=sub(index, 64)
      let T3:=add(T4, add( shl(12, and(shr(index2, scalar_v),1)), shl(8, and(shr(index2, scalar_u),1)) ))
      let index3:=sub(index2, 64)
      let T2:=add(T3,add( shl(11, and(shr(index3, scalar_v),1)), shl(7, and(shr(index3, scalar_u),1)) ))
      index:=sub(index3, 64)
      let T1:=add(T2,add( shl(10, and(shr(index, scalar_v),1)), shl(6, and(shr(index, scalar_u),1)) ))
      
      //index:=add(index,192), restore index, interleaved with loop
      
      //tbd: check validity of formulae with (0,1) to remove conditional jump
         if iszero(T1){
           Y:=sub(p, Y)
    
         continue
         }
       extcodecopy(dataPointer, T,T1, 64)
     }
     
     {
     
         /* Access to precomputed table using extcodecopy hack */
          
      // inlined EcZZ_AddN
      
      
      let y2:=addmod(mulmod(mload(add(T,32)), zzz,p),Y,p)  
      let T2:=addmod(mulmod(mload(T), zz,p),sub(p,X),p)  
      let T4:=mulmod(T2, T2, p)
      let T1:=mulmod(T4,T2,p)//
       zz:=mulmod(zz,T4,p) //zzz3=V*ZZ1
      zzz:= mulmod(zzz,T1,p) // W=UV/
      let zz1:=mulmod(X, T4, p)
      X:=addmod(addmod(mulmod(y2,y2, p), sub(p,T1),p ), mulmod(minus_2, zz1,p) ,p )
      Y:=addmod(mulmod(addmod(zz1, sub(p,X),p), y2, p), mulmod(Y, T1,p),p)
      
    
      }
      
     }//end loop
      mstore(add(T, 0x60),zz)
        
      //(X,Y)=ecZZ_SetAff(X,Y,zz, zzz);
      //T[0] = inverseModp_Hard(T[0], p); //1/zzz, inline modular inversion using precompile:
     // Define length of base, exponent and modulus. 0x20 == 32 bytes
      mstore(T, 0x20)
      mstore(add(T, 0x20), 0x20)
      mstore(add(T, 0x40), 0x20)
      // Define variables base, exponent and modulus
      //mstore(add(pointer, 0x60), u)
      mstore(add(T, 0x80), minus_2)
      mstore(add(T, 0xa0), p)
               
      // Call the precompiled contract 0x05 = ModExp
      if iszero(call(not(0), 0x05, 0, T, 0xc0, T, 0x20)) {
            revert(0, 0)
      }
      
      zz:=mload(T)
      X:=mulmod(X,zz,p)//X/zz
       }       
      }//end unchecked
    }
 function ecdsa_verify(
        bytes32 message,
        uint[2] memory rs,
        uint[2] memory Q
    ) internal  returns (bool) {
        uint256 n=FCL_Elliptic_ZZ.n;
        if (rs[0] == 0 || rs[0] >= n || rs[1] == 0) {
            return false;
        }
        
        
        if (!FCL_Elliptic_ZZ.ecAff_isOnCurve(Q[0], Q[1])) {
            return false;
        }
  	
        uint sInv = FCL_Elliptic_ZZ.FCL_nModInv(rs[1]);
        
        uint scalar_u=mulmod(uint(message), sInv, n);
        uint scalar_v= mulmod(rs[0], sInv, n);
        uint x1;
	
       x1=ecZZ_mulmuladd_S_asm(Q[0], Q[1],scalar_u, scalar_v);
       	
       	console.log("x1=",x1);
        
        return x1 == rs[0];
        
       }
     
  function test_Invariant_ecZZ_mulmuladd_S_asm(
       
    ) public returns (bool)
 {
//    string memory deployData = vm.readFile("test/vectors_wychproof/vec_sec256r1_valid.json");
    string memory deployData = vm.readFile("test/vectors_sage/vec_sage_ecdsa_sec256r1.json");
   
   
    uint256 wx=vm.parseJsonUint(deployData, ".NumberOfTests");
    console.log("NumberOfTests:",wx);
     uint256[2] memory key;
     key[0]=vm.parseJsonUint(deployData, ".keyx");
    console.log("key_x:", key[0]);
       key[1]=vm.parseJsonUint(deployData, ".keyy");
     console.log("key_y:", key[1]);
     bool res=FCL_Elliptic_ZZ.ecAff_isOnCurve(key[0],key[1]);
     console.log("Is key on curve:",res);
     
     //ReadField_th(deployData,".sigx_",0x1);
     
   uint256[2] memory rs;
   string memory snum="7";
  for(uint256 i=1;i<2;i++)
  {
   //snum=string(abi.encodePacked(uintToBytes(i)));
   string memory title=string(vm.parseJson(deployData, string.concat(".test_", snum)));
     console.log("test:",title);

   
    //title=string(vm.parseJson(deployData, string.concat(".test_", string(abi.encodePacked(uintToBytes(1))))));
   
  rs[0]=vm.parseJsonUint(deployData, string.concat(".sigx_", snum));
       console.log("sigx:",rs[0]);
  rs[1]=vm.parseJsonUint(deployData, string.concat(".sigy_", snum));
       console.log("sigy:",rs[1]);
  
    
   uint256 msg=vm.parseJsonUint(deployData, string.concat(".msg_", snum));
       console.log("msg:",msg);
       
       
     vm.prank(vm.addr(5));
     Wrap_ecdsa wrap=new Wrap_ecdsa();
     
      uint sInv = FCL_Elliptic_ZZ.FCL_nModInv(rs[1]);
        
        uint scalar_u=mulmod(uint(msg), sInv, FCL_Elliptic_ZZ.n);
        uint scalar_v= mulmod(rs[0], sInv, FCL_Elliptic_ZZ.n);
    
    console.log("u1:%x",scalar_u);
    console.log("u2:%x",scalar_v);
    console.log("-u2:%x",FCL_Elliptic_ZZ.n-scalar_v);
      
    res=ecdsa_verify(bytes32(msg), rs,key );
     console.log("Sig verif:",res);
     
     rs[1]=FCL_Elliptic_ZZ.n-rs[1];
     res=ecdsa_verify(bytes32(msg), rs,key );
     console.log("Sig verif:",res);

   }
  }
}
