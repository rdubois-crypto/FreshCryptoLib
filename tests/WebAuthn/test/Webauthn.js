const { expect } = require("chai");
const crypto = require("crypto");

const {precompute} = require("./precomputed");

function derToRS(der) {
  var offset = 3;
  var dataOffset;

  if (der[offset] == 0x21) {
    dataOffset = offset + 2;
  }
  else {
    dataOffset = offset + 1;
  }
  const r = der.slice(dataOffset, dataOffset + 32);
  offset = offset + der[offset] + 1 + 1
  if (der[offset] == 0x21) {
    dataOffset = offset + 2;
  }
  else {
    dataOffset = offset + 1;
  }
  const s = der.slice(dataOffset, dataOffset + 32);
  return [ r, s ]
}


describe("Webauthn", function() {

  it("Check message", async function() {
    
console.log("\n***************************************** \n Validating SSTORE 2 writing \n*****************************************" );

const precomputations_buff = Buffer.from(precompute, "hex");
  
  
 console.log("prec buff:",precomputations_buff);
  
  
console.log("length:",precomputations_buff.length);
 
const wo_core = await ethers.getContractFactory("Webauthn_prec2");

const wo = await wo_core.deploy(precomputations_buff);


console.log("\n***************************************** \n Direct table adressing \n*****************************************" );


const wo_table = await ethers.getContractFactory("BytecodeTable");
const deployed = await wo_table.deploy();

const replaceable = await ethers.getContractFactory("BytecodeTable");
const deployed_2 = await replaceable.deploy();

await deployed_2.deployed();

 console.log(
    `replaceable contract deployed to  ${deployed_2.address}`
  );


await network.provider.send("hardhat_setCode", [deployed_2.address, "0x"+precompute]);

const result_show = await deployed.Reader(deployed_2.address);



console.log("\n***************************************** \n Validating ECDSA Core verification \n*****************************************" );
  /* I Validation of Core ecdsa verification (no webauthn encoding) without precomputations */
  
const Ecdsa_core = await ethers.getContractFactory("Webauthn");

  let elliptic = require('elliptic');
  let sha3 = require('js-sha3');
 // let ec = new elliptic.ec('p256');
  let ec = new elliptic.ec('p256');
// let keyPair = ec.genKeyPair(); // Generate random keys
 // let keyPair = ec.keyFromPrivate( "97ddae0f3a25b92268175400149d65d6887b9cefaf28ea2c078e05cdc15a3c01");

 let keyPair = ec.keyFromPrivate(
 ethers.BigNumber.from(ethers.utils.randomBytes(32)).toHexString());


//console.log("random sk=",ethers.BigNumber.from(ethers.utils.randomBytes(32)).toHexString());

let privKey = keyPair.getPrivate("hex");
let pubKey = keyPair.getPublic();
console.log(`Private key: ${privKey}`);
console.log("Public key :", pubKey.encode("hex").substr(2));
console.log("Public key (compressed):",
    pubKey.encodeCompressed("hex"));

 const publicKey_ecdsa =Buffer.from(pubKey.encode("hex").substr(2), "hex");
 
 
console.log("publicKey_ecdsa:", publicKey_ecdsa);
 
 let msg = 'Message for signing';
let msgHash = sha3.keccak256(msg);
let hash=Buffer.from(msgHash, "hex");

let signature_js =  ec.sign(msgHash, privKey, "hex", {canonical: true});

console.log(`Msg: ${msg}`);
console.log(`Msg hash: ${msgHash}`);

console.log("hash:" , hash);


console.log("Signature:", signature_js.toDER('hex').toString(16));

console.log("Signature:",Buffer.from(signature_js.toDER('hex'), "hex"));

const ecdsaParsed = derToRS(Buffer.from(signature_js.toDER('hex'), "hex"));

console.log("Signature parsed:", ecdsaParsed);
    
    

    const ecdsa = await Ecdsa_core.deploy();
    await ecdsa.deployed();
    
    
    const result_ecdsa = await ecdsa.ecdsa_verif(hash,
        [ ethers.BigNumber.from("0x" + ecdsaParsed[0].toString('hex')), ethers.BigNumber.from("0x" + ecdsaParsed[1].toString('hex'))],
        [ ethers.BigNumber.from("0x" + publicKey_ecdsa.slice(0, 32).toString('hex')), ethers.BigNumber.from("0x" + publicKey_ecdsa.slice(32).toString('hex'))]
    );
    
    await result_ecdsa.wait();
    
    
console.log("\n***************************************** \n Validating WebAuthn with XYZZ coordinates \n*****************************************" );
  /* II Validation of Webauthn verification without precomputations */
    const Webauthn = await ethers.getContractFactory("Webauthn");
   
    const publicKey = Buffer.from("fdf8bce27f54e06f3aee3b6a542db1ab1f2418d7370a78b150d06965f942b14a470cdee69ab50e610c39b840681bf816b030f4a0a5d5af02ce27dcce6bede89f", "hex");
    const signature = Buffer.from("30440220655c9a457615aac594d92fb6d842f0e910e5ee6677cddbcddaea624f3203f0e702207b71a302b06c91a52b9c4ba5a7fb85226738b02c144e8ee177d034022a79c946", "hex");
    
console.log("Signature:", signature);
    const authenticatorData = Buffer.from("f8e4b678e1c62f7355266eaa4dc1148573440937063a46d848da1e25babbd20b010000004d", "hex");
    const clientData = Buffer.from("7b2274797065223a22776562617574686e2e676574222c226368616c6c656e6765223a224e546f2d3161424547526e78786a6d6b61544865687972444e5833697a6c7169316f776d4f643955474a30222c226f726967696e223a2268747470733a2f2f66726573682e6c65646765722e636f6d222c2263726f73734f726967696e223a66616c73657d", "hex");
    const clientChallenge = Buffer.from("353a3ed5a0441919f1c639a46931de872ac3357de2ce5aa2d68c2639df54189d", "hex");

    // Get on the "challenge" key value following "
    const challengeOffset = clientData.indexOf("226368616c6c656e6765223a", 0, "hex") + 12 + 1;    
    const signatureParsed = derToRS(signature);


console.log("Signature parsed:", signatureParsed);
	/*
    const result = await webauthn.checkSignature(authenticatorData, 0x01, clientData, clientChallenge, challengeOffset,
        [ ethers.BigNumber.from("0x" + signatureParsed[0].toString('hex')), ethers.BigNumber.from("0x" + signatureParsed[1].toString('hex'))],
        [ ethers.BigNumber.from("0x" + publicKey.slice(0, 32).toString('hex')), ethers.BigNumber.from("0x" + publicKey.slice(32).toString('hex'))]
    );
    expect(result);
*/

//uncomment for no precomputation validation



    const webauthn = await Webauthn.deploy();
    await webauthn.deployed();

    const result = await webauthn.validate(authenticatorData, 0x01, clientData, clientChallenge, challengeOffset,
        [ ethers.BigNumber.from("0x" + signatureParsed[0].toString('hex')), ethers.BigNumber.from("0x" + signatureParsed[1].toString('hex'))],
        [ ethers.BigNumber.from("0x" + publicKey.slice(0, 32).toString('hex')), ethers.BigNumber.from("0x" + publicKey.slice(32).toString('hex'))]
    );
    await result.wait();
   

//uncomment for with precomputation validation
    
    
/*
console.log("\n***************************************** \n Validating WebAuthn with XYZZ coordinates and Precomputations in memory \n*****************************************" );    
  
      const Webauthn_prec = await ethers.getContractFactory("Webauthn_prec");
   
    const webauthn_prec = await Webauthn_prec.deploy([ ethers.BigNumber.from("0x" + publicKey.slice(0, 32).toString('hex')), ethers.BigNumber.from("0x" + publicKey.slice(32).toString('hex'))]);
    
    
      const webauthn_prec2 = await webauthn_prec.deploy_part2([ ethers.BigNumber.from("0x" + publicKey.slice(0, 32).toString('hex')), ethers.BigNumber.from("0x" + publicKey.slice(32).toString('hex'))]);
      
      
 
      const result2 = await webauthn_prec.validate_prec(authenticatorData, 0x01, clientData, clientChallenge, challengeOffset,
        [ ethers.BigNumber.from("0x" + signatureParsed[0].toString('hex')), ethers.BigNumber.from("0x" + signatureParsed[1].toString('hex'))]
        
    );
    
 await result2.wait();
    
*/

console.log("\n***************************************** \n Validating WebAuthn with XYZZ coordinates and Precomputations in precomputed contract with sstore2 \n*****************************************" );    
  /* III Validation of Webauthn verification with precomputations */
  
  
  
 console.log("prec buff:",precomputations_buff);
  
  
console.log("length:",precomputations_buff.length);
 
const wo_core2 = await ethers.getContractFactory("Webauthn_prec2");

const wo2 = await wo_core2.deploy(precomputations_buff);


      
      
 
      const result3 = await wo2.validate_prec(authenticatorData, 0x01, clientData, clientChallenge, challengeOffset,
        [ ethers.BigNumber.from("0x" + signatureParsed[0].toString('hex')), ethers.BigNumber.from("0x" + signatureParsed[1].toString('hex'))]
        
    );
    
 await result3.wait();
        

console.log("\n***************************************** \n Validating WebAuthn with XYZZ coordinates and Precomputations in precomputed contract with extcodecopy \n*****************************************" );    
  /* IV Validation of Webauthn verification with precomputations */
  



  
const wo_core3 = await ethers.getContractFactory("Webauthn_prec3");
 
 
const wo3 = await wo_core3.deploy(deployed_2.address);
  
const codebyte = await hre.network.provider.send("eth_getCode", [
        wo3.address,
    ]);
    
console.log("\n Bytecode size:\n",codebyte.length); 

    
      const result4 = await wo3.validate_prec(authenticatorData, 0x01, clientData, clientChallenge, challengeOffset,
        [ ethers.BigNumber.from("0x" + signatureParsed[0].toString('hex')), ethers.BigNumber.from("0x" + signatureParsed[1].toString('hex'))]
        
    );
    
 await result4.wait();      
 
    const result5=  wo3.ecdsa_verif_prec(clientChallenge, 
        [ ethers.BigNumber.from("0x" + signatureParsed[0].toString('hex')), ethers.BigNumber.from("0x" + signatureParsed[1].toString('hex'))]
        
    );
   
  })
  
//console.log("\n Bytecode:\n",codebyte); 
  
});
