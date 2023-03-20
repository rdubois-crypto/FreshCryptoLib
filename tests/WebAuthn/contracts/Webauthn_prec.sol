// SPDX-License-Identifier: Apache-2.0
pragma solidity ^0.8.0;

import {Base64URL} from "./Base64URL.sol";
import {Ec_ZZ} from "./Elliptic_ZZ.sol";

import {FCL_Elliptic_ZZ} from "./FCL_elliptic.sol";

import "hardhat/console.sol";

import  "solmate/src/utils/SSTORE2.sol";



error InvalidAuthenticatorData();
error InvalidClientData();
error InvalidSignature();

contract replaceable{
 string constant x="replace me please";
}

contract BytecodeTable {
 //precomputations
   
  
  function Reader(address selfe) public
  {
    uint[2] memory px;
    
    
    console.log("address=",uint256(uint160(selfe)));
    
     
    
    for(uint i=0;i<256;i++){  
    uint offset=64*i;
    assembly{
       extcodecopy(selfe, px, offset, 64)
      }
      //console.log("Test on curve of point ",i, px[0], px[1]);
      //console.log(Ec_ZZ.ecAff_isOnCurve(px[0], px[1]));
  }    
  }
}


contract Webauthn_prec3{

 address   dataPointer;
 uint256 public counter;
  
  constructor( address Precompute_contract)
  { 
	dataPointer=Precompute_contract;
  }
  
  

    function checkSignature_prec(
        bytes memory authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes memory clientData,
        bytes32 clientChallenge,
        uint clientChallengeDataOffset,
        uint[2] memory rs       
    ) public  returns (bool) {
        // Let the caller check if User Presence (0x01) or User Verification (0x04) are set
        if (
            (authenticatorData[32] & authenticatorDataFlagMask) !=
            authenticatorDataFlagMask
        ) {
            revert InvalidAuthenticatorData();
        }
        // Verify that clientData commits to the expected client challenge
        string memory challengeEncoded = Base64URL.encode32(
            abi.encodePacked(clientChallenge)
        );
        bytes memory challengeExtracted = new bytes(
            bytes(challengeEncoded).length
        );
        copyBytes(
            clientData,
            clientChallengeDataOffset,
            challengeExtracted.length,
            challengeExtracted,
            0
        );
        if (
            keccak256(abi.encodePacked(bytes(challengeEncoded))) !=
            keccak256(abi.encodePacked(challengeExtracted))
        ) {
            revert InvalidClientData();
        }      
        // Verify the signature over sha256(authenticatorData || sha256(clientData))
        bytes memory verifyData = new bytes(authenticatorData.length + 32);
        copyBytes(
            authenticatorData,
            0,
            authenticatorData.length,
            verifyData,
            0
        );
        copyBytes(
            abi.encodePacked(sha256(clientData)),
            0,
            32,
            verifyData,
            authenticatorData.length
        );
        bytes32 message = sha256(verifyData);
	bool result=FCL_Elliptic_ZZ.ecdsa_precomputed_verify(message, rs,  dataPointer);
	console.log("result= %s", result);

        return result;
    }
    
 function validate_prec(
        bytes memory authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes memory clientData,
        bytes32 clientChallenge,
        uint clientChallengeDataOffset,
        uint[2] memory rs
       
    ) public {
        if (
            !checkSignature_prec(
                authenticatorData,
                authenticatorDataFlagMask,
                clientData,
                clientChallenge,
                clientChallengeDataOffset,
                rs              
            )
        ) {
            revert InvalidSignature();
        }
        counter++;
    }
    
    
    function copyBytes(
        bytes memory _from,
        uint _fromOffset,
        uint _length,
        bytes memory _to,
        uint _toOffset
    ) internal pure returns (bytes memory _copiedBytes) {
        uint minLength = _length + _toOffset;
        require(_to.length >= minLength); // Buffer too small. Should be a better way?
        uint i = 32 + _fromOffset; // NOTE: the offset 32 is added to skip the `size` field of both bytes variables
        uint j = 32 + _toOffset;
        while (i < (32 + _fromOffset + _length)) {
            assembly {
                let tmp := mload(add(_from, i))
                mstore(add(_to, j), tmp)
            }
            i += 32;
            j += 32;
        }
        return _to;
    }



  
}


contract Webauthn_prec2 {

  address   dataPointer;
  uint256 public counter;
 
      
  constructor( bytes memory Shamir8_ss2)
  {
     uint taille;
     assembly{
     	taille:=mload(Shamir8_ss2)
     	}
      
     console.log("Ecriture de la table input  dans le contrat, taille :", taille);
   
     dataPointer=SSTORE2.write(Shamir8_ss2);
     bytes memory ec_point=new bytes(64);
    
     uint256[2] memory xy;
     
     uint offset;
     uint endoffset;
     uint px;
     uint py;
     address cpy_datap=dataPointer;
     
     /*
     console.log("Relecture de la table :", taille);

     
     for(uint i=0;i<256;i++)
     {
      offset=64*i;
      endoffset=offset+64;
      
      ec_point=SSTORE2.read(dataPointer, 64*i, 64*i+64);
      assembly{
        px:=mload(add(ec_point,32))
        py:=mload(add(ec_point,64))
        
      }
      
      console.log("Read1 :", px,py);
      console.log(Ec_ZZ.ecAff_isOnCurve(px, py));
     } */
  }

  

    function checkSignature_prec(
        bytes memory authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes memory clientData,
        bytes32 clientChallenge,
        uint clientChallengeDataOffset,
        uint[2] memory rs       
    ) public  returns (bool) {
        // Let the caller check if User Presence (0x01) or User Verification (0x04) are set
        if (
            (authenticatorData[32] & authenticatorDataFlagMask) !=
            authenticatorDataFlagMask
        ) {
            revert InvalidAuthenticatorData();
        }
        // Verify that clientData commits to the expected client challenge
        string memory challengeEncoded = Base64URL.encode32(
            abi.encodePacked(clientChallenge)
        );
        bytes memory challengeExtracted = new bytes(
            bytes(challengeEncoded).length
        );
        copyBytes(
            clientData,
            clientChallengeDataOffset,
            challengeExtracted.length,
            challengeExtracted,
            0
        );
        if (
            keccak256(abi.encodePacked(bytes(challengeEncoded))) !=
            keccak256(abi.encodePacked(challengeExtracted))
        ) {
            revert InvalidClientData();
        }      
        // Verify the signature over sha256(authenticatorData || sha256(clientData))
        bytes memory verifyData = new bytes(authenticatorData.length + 32);
        copyBytes(
            authenticatorData,
            0,
            authenticatorData.length,
            verifyData,
            0
        );
        copyBytes(
            abi.encodePacked(sha256(clientData)),
            0,
            32,
            verifyData,
            authenticatorData.length
        );
        bytes32 message = sha256(verifyData);
	bool result=Ec_ZZ.validateSignature_Precomputed_ss2(message, rs,  dataPointer);
	console.log("result= %s", result);

        return result;
    }
    
 function validate_prec(
        bytes memory authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes memory clientData,
        bytes32 clientChallenge,
        uint clientChallengeDataOffset,
        uint[2] memory rs
       
    ) public {
        if (
            !checkSignature_prec(
                authenticatorData,
                authenticatorDataFlagMask,
                clientData,
                clientChallenge,
                clientChallengeDataOffset,
                rs              
            )
        ) {
            revert InvalidSignature();
        }
        counter++;
    }
    
    
    function copyBytes(
        bytes memory _from,
        uint _fromOffset,
        uint _length,
        bytes memory _to,
        uint _toOffset
    ) internal pure returns (bytes memory _copiedBytes) {
        uint minLength = _length + _toOffset;
        require(_to.length >= minLength); // Buffer too small. Should be a better way?
        uint i = 32 + _fromOffset; // NOTE: the offset 32 is added to skip the `size` field of both bytes variables
        uint j = 32 + _toOffset;
        while (i < (32 + _fromOffset + _length)) {
            assembly {
                let tmp := mload(add(_from, i))
                mstore(add(_to, j), tmp)
            }
            i += 32;
            j += 32;
        }
        return _to;
    }



}

contract Webauthn_prec {
    uint256 public counter;
     
    
    //precomputations
    uint[2][256] public Shamir8;
    
    constructor(uint[2] memory Q){
    
    console.log("Precompute for Public key:", Q[0], Q[1]);
      	
    Shamir8=Ec_ZZ.Precalc_Shamir8_part1(Q);
    console.log("Precompute part 1 done");
    }
    
    
    
    function deploy_part2(uint[2] memory Q) public  returns ( uint[2][256] memory res) 
    {
     unchecked{  
     
     uint[2][128] memory Prec2=Ec_ZZ.Precalc_Shamir8_part2(Q);
     
       for(uint i=128;i<256;i++)
     {       
        Shamir8[i][0]=Prec2[i-128][0];
        Shamir8[i][1]=Prec2[i-128][1];
     }
     
        console.log("Precompute part 2 done:");
  	
    }
    
 
     /* 256 point, 2 coordinates, 32 bytes each =16384 bytes*/
     //bytes memory Shamir8_ss2= new bytes(1);
   
     for(uint i=0; i<256; i++)
     {
         console.log("\"", Shamir8[i][0],"\",");
         console.log("\"", Shamir8[i][1],"\",");
	//res[i][0]=Shamir8[i][0];
	//res[i][1]=Shamir8[i][1];
     }
     
     //cette ligne pète
     //dataPointer=SSTORE2.write(Shamir8_ss2);
     return res=Shamir8;
    }

    function checkSignature_prec(
        bytes memory authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes memory clientData,
        bytes32 clientChallenge,
        uint clientChallengeDataOffset,
        uint[2] memory rs,
        uint[2][256] memory Shamir8
    ) public  returns (bool) {
        // Let the caller check if User Presence (0x01) or User Verification (0x04) are set
        if (
            (authenticatorData[32] & authenticatorDataFlagMask) !=
            authenticatorDataFlagMask
        ) {
            revert InvalidAuthenticatorData();
        }
        // Verify that clientData commits to the expected client challenge
        string memory challengeEncoded = Base64URL.encode32(
            abi.encodePacked(clientChallenge)
        );
        bytes memory challengeExtracted = new bytes(
            bytes(challengeEncoded).length
        );
        copyBytes(
            clientData,
            clientChallengeDataOffset,
            challengeExtracted.length,
            challengeExtracted,
            0
        );
        if (
            keccak256(abi.encodePacked(bytes(challengeEncoded))) !=
            keccak256(abi.encodePacked(challengeExtracted))
        ) {
            revert InvalidClientData();
        }      
        // Verify the signature over sha256(authenticatorData || sha256(clientData))
        bytes memory verifyData = new bytes(authenticatorData.length + 32);
        copyBytes(
            authenticatorData,
            0,
            authenticatorData.length,
            verifyData,
            0
        );
        copyBytes(
            abi.encodePacked(sha256(clientData)),
            0,
            32,
            verifyData,
            authenticatorData.length
        );
        bytes32 message = sha256(verifyData);
	bool result=Ec_ZZ.validateSignature_Precomputed(message, rs,  Shamir8);
	console.log("result= %s", result);

        return result;
    }
    
 function validate_prec(
        bytes memory authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes memory clientData,
        bytes32 clientChallenge,
        uint clientChallengeDataOffset,
        uint[2] memory rs
       
    ) public {
        if (
            !checkSignature_prec(
                authenticatorData,
                authenticatorDataFlagMask,
                clientData,
                clientChallenge,
                clientChallengeDataOffset,
                rs,
                Shamir8
            )
        ) {
            revert InvalidSignature();
        }
        counter++;
    }
    
    
    function copyBytes(
        bytes memory _from,
        uint _fromOffset,
        uint _length,
        bytes memory _to,
        uint _toOffset
    ) internal pure returns (bytes memory _copiedBytes) {
        uint minLength = _length + _toOffset;
        require(_to.length >= minLength); // Buffer too small. Should be a better way?
        uint i = 32 + _fromOffset; // NOTE: the offset 32 is added to skip the `size` field of both bytes variables
        uint j = 32 + _toOffset;
        while (i < (32 + _fromOffset + _length)) {
            assembly {
                let tmp := mload(add(_from, i))
                mstore(add(_to, j), tmp)
            }
            i += 32;
            j += 32;
        }
        return _to;
    }
}
