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

contract wrap_ecdsa_precal{
 address precomputations;

function wrap_ecdsa_core( bytes32 message,
        uint[2] calldata rs,
        uint[2] calldata Q) public returns (bool)
 {
  return FCL_Elliptic_ZZ.ecdsa_precomputed_verify(message,rs, precomputations);
 }
 
 constructor(address bytecode){
   precomputations=bytecode;
 }
}

contract EcdsaTest is Test {

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
  
  function do_precalc() public returns (bool)
  {
    bytes memory bytecode ;
    
    string[] memory inputs = new string[](1);
    inputs[0] = "echo toto  ";
/*    inputs[1] = "114874632398302156264159990279427641021947882640101801130664833947273521181002";
    inputs[2] = ";C1=";
    inputs[3] = "32136952818958550240756825111900051564117520891182470183735244184006536587423'";
    inputs[4] = ";load(\"test/FCL_ecdsa_precompute.sage\")'>>fcl_ecdsa_precbytecode.js" ;*/
	
    vm.ffi(inputs);    
    console.log("precalc done");
  }
  
  function read_precalc() public returns (bool)
  {
   
  }
  
  
  function test_Invariant_edge() public returns (bool)
  {
   //choose Q=2P, then verify duplication is ok
   uint256[4] memory Q;
   (Q[0], Q[1],Q[2], Q[3])=FCL_Elliptic_ZZ.ecZZ_Dbl(gx, gy,1,1);
   uint256[4] memory _4P;
   (_4P[0], _4P[1],_4P[2], _4P[3])=FCL_Elliptic_ZZ.ecZZ_Dbl(Q[0], Q[1],Q[2], Q[3]);
   uint256 _4P_res1;
   
   (_4P_res1,)=FCL_Elliptic_ZZ.ecZZ_SetAff(_4P[0], _4P[1],_4P[2], _4P[3]);
   
   
   uint256  _4P_res2=FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(gx,gy, 4, 0);
   assertEq(_4P_res1, _4P_res2);
   
   uint256[2] memory nQ;
   (nQ[0],nQ[1])=FCL_Elliptic_ZZ.ecZZ_SetAff(Q[0], Q[1],Q[2], Q[3]);
   uint256 _4P_res3=FCL_Elliptic_ZZ.ecZZ_mulmuladd_S_asm(nQ[0],nQ[1], 2, 1);
  
   assertEq(_4P_res1, _4P_res3);
   
  }
  
  //testing Wychproof vectors
  function test_Invariant_ecZZ_mulmuladd_S_asm(
       
    ) public returns (bool)
 {
    string memory deployData = vm.readFile("test/vectors_wychproof/vec_sec256r1_valid.json");
    //string memory deployData = vm.readFile("test/vectors_sage/vec_sage_ecdsa_sec256r1.json");
   
   
    uint256 wx=vm.parseJsonUint(deployData, ".NumberOfTests");
    console.log("NumberOfTests:",wx);
    uint256[2] memory key;
    key[0]=vm.parseJsonUint(deployData, ".keyx");
    console.log("key_x:", key[0]);
    key[1]=vm.parseJsonUint(deployData, ".keyy");
    console.log("key_y:", key[1]);
    bool res=FCL_Elliptic_ZZ.ecAff_isOnCurve(key[0],key[1]);
    assertEq(res,true);
     
     console.log("Is key on curve:",res);
     
     //ReadField_th(deployData,".sigx_",0x1);
     
     
     
    uint256[2] memory rs;
    string memory snum="1";
    for(uint256 i=1;i<10;i++)
    {
  
   
   //snum=string(abi.encodePacked(uintToBytes(i)));
  
  
      string memory title=string(vm.parseJson(deployData, string.concat(".test_", snum)));
      console.log("\n test:",snum, title);
       console.log("\n ||");
    
  
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
    
      
       res=wrap.wrap_ecdsa_core(bytes32(msg), rs,key );
       console.log("Sig verif:",res);
     
     /*
     console.log("-u2:%x",FCL_Elliptic_ZZ.n-scalar_v); 
     rs[1]=FCL_Elliptic_ZZ.n-rs[1];
     res=ecdsa_verify(bytes32(msg), rs,key );
     console.log("Sig verif:",res);*/
       assembly{
        mstore(add(snum,32),add(0x0100000000000000000000000000000000000000000000000000000000000000, mload(add(snum,32))))
      }
   }
  }
  
  
  
}
