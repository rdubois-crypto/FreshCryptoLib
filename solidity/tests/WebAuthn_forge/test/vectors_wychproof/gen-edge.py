# A little script to reformat valid wychproof input to ecdsa_core input

import hashlib
import json
from re import sub
import binascii
from pyasn1.codec.der.decoder import decode as der_decoder
from pyasn1.type.univ import Sequence
from pyasn1.type.univ import Integer
from pyasn1.type.namedtype import NamedTypes
from pyasn1.type.namedtype import NamedType
import string



TEST_CASE = """
@view
func test_{title}{{syscall_ptr : felt*, range_check_ptr, pedersen_ptr : HashBuiltin*}}() {
    {expect_revert}
    let public_key_pt = EcPoint(
        BigInt3({x0},{x1},{x2}),
        BigInt3({y0},{y1},{y2}));
    let r = BigInt3({r0},{r1},{r2});
    let s = BigInt3({s0},{s1},{s2});
    let msg_hash = BigInt3({m0},{m1},{m2});
    verify_ecdsa(public_key_pt=public_key_pt, msg_hash=msg_hash, r=r, s=s);
    return ();
}
"""



def snake_case(s):
    return '_'.join(
        sub('([A-Z][a-z]+)', r' \1',
            sub('([A-Z]+)', r' \1',
                s.replace('-', ' '))).split()).lower()


class DERSig(Sequence):
    componentType = NamedTypes(
        NamedType('r', Integer()),
        NamedType('s', Integer())
    )


input = open('ecdsa_secp256r1_sha256_test.json')
f = open("ecdsa_secp256r1_sha256_test.cairo", "w")

data = json.load(input)


test="\n";

numvec=0;
for tg in data['testGroups']:
    print("{\n \"keyx\":",int(tg['key']['wx'], 16), ",\"keyy\":", int(tg['key']['wy'], 16), ",")

    
    for j, tc in enumerate(tg['tests']):
        try:
        
            title = str(
                j) + "_" + snake_case(tc['comment'].translate(str.maketrans('', '', string.punctuation)))
                
           
            
            msg = hashlib.sha256(binascii.unhexlify(tc["msg"])).hexdigest()
            
            
            
            sig, rest = der_decoder(
                binascii.unhexlify(tc["sig"]), asn1Spec=DERSig())
                
                
                
           # if len(rest) != 0:
            #    raise Exception('Bad encoding')

            expect_revert = ""
            if tc["result"] == "invalid":
                numvec=numvec+1;
                charvec=str(numvec);
                print(" \"test_"+charvec+"\":\"", title, "\", \"msg_"+charvec+"\": \"0x"+msg+"\"",", \"sigx_"+charvec+"\":", int(sig['r']), ", \"sigy_"+charvec+"\":", int(sig['s']),",");
                
        except Exception as err:
            title

    break

print(" \n \"NumberOfTests\":",numvec,"\n}"); 

f.write(test)
f.close()
input.close()
