from fido2.hid import CtapHidDevice
from fido2.client import Fido2Client, UserInteraction 
from fido2.server import Fido2Server
from fido2.webauthn import PublicKeyCredentialParameters, PublicKeyCredentialType
from fido2.cose import ES256
import sys
import ctypes
import getpass
import binascii

class CliInteraction(UserInteraction):
    def prompt_up(self):
        print("\nTouch your authenticator device now...\n")

    def request_pin(self, permissions, rd_id):
        return getpass.getpass("Enter PIN: ")

    def request_uv(self, permissions, rd_id):
        print("User Verification required.")
        return True

CHALLENGE = binascii.unhexlify("353a3ed5a0441919f1c639a46931de872ac3357de2ce5aa2d68c2639df54189d")

device = next(CtapHidDevice.list_devices(), None)
client = Fido2Client(device, "https://fresh.ledger.com", user_interaction=CliInteraction())
server = Fido2Server({"id": "fresh.ledger.com", "name" : "Webauthn"})
server.allowed_algorithms = [ PublicKeyCredentialParameters(PublicKeyCredentialType.PUBLIC_KEY, ES256.ALGORITHM) ]
user = {"id":b"test", "name" : "test"}
create_options, state = server.register_begin(user, user_verification="discouraged", authenticator_attachment="cross-platform")
result = client.make_credential(create_options["publicKey"])
auth_data = server.register_complete(state, result.client_data, result.attestation_object)
credentials = [auth_data.credential_data]
publicKey = credentials[0].public_key[-2] + credentials[0].public_key[-3]
print("public key " + binascii.hexlify(publicKey).decode('utf-8'))

request_options, state = server.authenticate_begin(credentials, user_verification="discouraged", challenge=CHALLENGE)
result = client.get_assertion(request_options["publicKey"])
result = result.get_response(0)

server.authenticate_complete(
    state,
    credentials,
    result.credential_id,
    result.client_data,
    result.authenticator_data,
    result.signature,
)

print("authenticator data " + binascii.hexlify(result.authenticator_data).decode('utf-8'))
print("client data " + binascii.hexlify(result.client_data).decode('utf-8'))
print("client data " + str(bytes(result.client_data).decode('utf-8')))
print("client data hash " + binascii.hexlify(result.client_data.hash).decode('utf-8'))
print("signature " + binascii.hexlify(result.signature).decode('utf-8'))


