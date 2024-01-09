// SPDX-License-Identifier: MIT

pragma solidity ^0.8.0;

import {Test, console2} from "forge-std/Test.sol";
import {FCL_WebAuthn} from "@solidity/FCL_Webauthn.sol";
import {FCL_Elliptic_ZZ} from "@solidity/FCL_elliptic.sol";
import {Base64Url} from "@solidity/utils/Base64Url.sol";

/// @dev This test contract is used to show the need for Base64Url encoding in the WebAuthn_Format function
/**
 * @title FCL_WebAuthn_Base64Url_Test
 * @author evmBrahmin
 * @notice Mock data is used in these tests. Here I will explain how this mock data was derived to ensure its accuracy.
 * This is the response object from the WebAuthn Authentication API call:
 * {
 *      response:
 *         authenticatorData: "SZYN5YgOjGh0NBcPZHZgW4_krrmihjLHmVzzuoMdl2MFAAAAAA"
 *         clientDataJSON: "eyJ0eXBlIjoid2ViYXV0aG4uZ2V0IiwiY2hhbGxlbmdlIjoibkNMX1h5SHd1QnNSUG1QMzIyMnBULTN2RWJJUm0wQ0l1SlprLTVvOHRsZyIsIm9yaWdpbiI6Imh0dHA6Ly9sb2NhbGhvc3Q6MzAwMCIsImNyb3NzT3JpZ2luIjpmYWxzZX0"
 *         signature: "MEUCIECRE8S97mXV1Dwqqp3uF_CW3c6XvQMQrkrgjnx1lVnLAiEA00ucboY5T_qXn5MJdpYyzvid-8MROOS9-Q3QRPvqsl4"
 *         userHandle:"cUtqZTRBNGk0TTdDTFBGTVE4UFVOam5PU3RsRUlMdDRyOWpMdG00amtDRT0"
 *     }
 *  The authenticatorDataMock is the hex representation of the bytes retrieved by Base64URL decoding the authenticatorData from the response object above.
 *  The clientDataMock provided below is retrieved by Base64URL decoding the clientDataJSON from the response object above, and then converting it to the hex representation of those bytes.
 *  The challengeMock is the hex representation of the bytes retrieved by Base64URL decoding the clientDataJSON from the response object above.
 *  The rsMock is the r and s values of the signature. Note that the signature provided above is the Base64URL encoding of the DER encoding of the signature.
 *  The QMock is the x and y values of the public key of the key pair used to sign the data. This is generated during the registration process
 */
contract FCL_Webauthn_Base64Url is Test {
    Helper public helper;
    WebAuthn_base64URL public webauthn;

    function setUp() public {
        helper = new Helper();
        webauthn = new WebAuthn_base64URL();
    }

    function test_base64URL_format() external view {
        // mock data, see notice above for details on how this was derived
        bytes memory authenticatorDataMock =
            hex"49960de5880e8c687434170f6476605b8fe4aeb9a28632c7995cf3ba831d97630500000000";
        bytes1 authenticatorDataFlagMaskMock = hex"01";
        bytes memory clientDataMock =
            hex"7b2274797065223a22776562617574686e2e676574222c226368616c6c656e6765223a226e434c5f5879487775427352506d503332323270542d3376456249526d304349754a5a6b2d356f38746c67222c226f726967696e223a22687474703a2f2f6c6f63616c686f73743a33303030222c2263726f73734f726967696e223a66616c73657d";
        bytes32 clientChallengeMock = hex"9c22ff5f21f0b81b113e63f7db6da94fedef11b2119b4088b89664fb9a3cb658";
        uint256 clientChallengeDataOffsetMock = 36;
        uint256[2] memory rsMock = [
            29204351571054144655406732941989447033933540609767730374087271035220690033099,
            95571604233087243530638576272724546495790691135210520108106485630494705365598
        ];
        /*uint256[2] memory QMock = [
            81682839938742543555082486110423905347664508641828409199518693630277409128887,
            17376078742205537968081477783760696603470208142043731074843150602836373389861
        ];*/

        // if this function call doesn't revert with an error, we believe it to have passed
        // the next test ensures that data is properly formated and successfully verified
        bytes32 message = webauthn.format(
            authenticatorDataMock,
            authenticatorDataFlagMaskMock,
            clientDataMock,
            clientChallengeMock,
            clientChallengeDataOffsetMock,
            rsMock
        );
        console2.logString("Message: ");
        console2.logBytes32(message);
    }

    function test_webauthn_Base64URL_checkSignature() external {
        // mock data, see notice above for details on how this was derived
        bytes memory authenticatorDataMock =
            hex"49960de5880e8c687434170f6476605b8fe4aeb9a28632c7995cf3ba831d97630500000000";
        bytes1 authenticatorDataFlagMaskMock = hex"01";
        bytes memory clientDataMock =
            hex"7b2274797065223a22776562617574686e2e676574222c226368616c6c656e6765223a226e434c5f5879487775427352506d503332323270542d3376456249526d304349754a5a6b2d356f38746c67222c226f726967696e223a22687474703a2f2f6c6f63616c686f73743a33303030222c2263726f73734f726967696e223a66616c73657d";
        bytes32 clientChallengeMock = hex"9c22ff5f21f0b81b113e63f7db6da94fedef11b2119b4088b89664fb9a3cb658";
        uint256 clientChallengeDataOffsetMock = 36;
        uint256[2] memory rsMock = [
            29204351571054144655406732941989447033933540609767730374087271035220690033099,
            95571604233087243530638576272724546495790691135210520108106485630494705365598
        ];
        uint256[2] memory QMock = [
            81682839938742543555082486110423905347664508641828409199518693630277409128887,
            17376078742205537968081477783760696603470208142043731074843150602836373389861
        ];

        bool result = webauthn.verify(
            authenticatorDataMock,
            authenticatorDataFlagMaskMock,
            clientDataMock,
            clientChallengeMock,
            clientChallengeDataOffsetMock,
            rsMock,
            QMock
        );

        assertTrue(result, "authentication failed");
    }

    function test_webauthn_format_details() public view {
        // mock data, see notice above for details on how this was derived
        bytes memory authenticatorDataMock =
            hex"49960de5880e8c687434170f6476605b8fe4aeb9a28632c7995cf3ba831d97630500000000";
        bytes1 authenticatorDataFlagMaskMock = hex"01";
        bytes memory clientDataMock =
            hex"7b2274797065223a22776562617574686e2e676574222c226368616c6c656e6765223a226e434c5f5879487775427352506d503332323270542d3376456249526d304349754a5a6b2d356f38746c67222c226f726967696e223a22687474703a2f2f6c6f63616c686f73743a33303030222c2263726f73734f726967696e223a66616c73657d";
        bytes32 clientChallengeMock = hex"9c22ff5f21f0b81b113e63f7db6da94fedef11b2119b4088b89664fb9a3cb658";
        uint256 clientChallengeDataOffsetMock = 36;
        uint256[2] memory rsMock = [
            29204351571054144655406732941989447033933540609767730374087271035220690033099,
            95571604233087243530638576272724546495790691135210520108106485630494705365598
        ];
        uint256[2] memory QMock = [
            81682839938742543555082486110423905347664508641828409199518693630277409128887,
            17376078742205537968081477783760696603470208142043731074843150602836373389861
        ];

        // call helper function for logs of the new encoding process
        bytes32 message = helper.WebAuthn_format(
            authenticatorDataMock, // authenticator data
            authenticatorDataFlagMaskMock, // authenticator data flag mask
            clientDataMock, // client data
            clientChallengeMock, // client challenge
            clientChallengeDataOffsetMock, // client challenge data offset
            rsMock // signature r and s
        );
    }
}

// Contract implementing the FCL_Webauthn library functions to test
contract WebAuthn_base64URL {
    function format(
        bytes calldata authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes calldata clientData,
        bytes32 clientChallenge,
        uint256 clientChallengeDataOffset,
        uint256[2] calldata rs
    ) external pure returns (bytes32 message) {
        message = FCL_WebAuthn.WebAuthn_format(
            authenticatorData, // authenticator data
            authenticatorDataFlagMask, // authenticator data flag mask
            clientData, // client data
            clientChallenge, // client challenge
            clientChallengeDataOffset, // client challenge data offset
            rs // signature r and s
        );
    }

    function verify(
        bytes calldata authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes calldata clientData,
        bytes32 clientChallenge,
        uint256 clientChallengeDataOffset,
        uint256[2] calldata rs,
        uint256[2] calldata Q
    ) external view returns (bool result) {
        result = FCL_WebAuthn.checkSignature(
            authenticatorData, // authenticator data
            authenticatorDataFlagMask, // authenticator data flag mask
            clientData, // client data
            clientChallenge, // client challenge
            clientChallengeDataOffset, // client challenge data offset
            rs, // signature r and s
            Q // public key
        );
    }
}

// A contract with the logic directly embedded to provide logs to see how the encoding is working
contract Helper is Test {
    error InvalidAuthenticatorData();
    error InvalidClientData();
    error InvalidSignature();

    uint256 constant n = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551;

    string internal constant ENCODING_TABLE = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";

    // removed pure due to console2.logs
    function WebAuthn_format(
        bytes calldata authenticatorData,
        bytes1 authenticatorDataFlagMask,
        bytes calldata clientData,
        bytes32 clientChallenge,
        uint256 clientChallengeDataOffset,
        uint256[2] calldata // rs
    ) public view returns (bytes32 result) {
        // Let the caller check if User Presence (0x01) or User Verification (0x04) are set
        {
            if ((authenticatorData[32] & authenticatorDataFlagMask) != authenticatorDataFlagMask) {
                revert InvalidAuthenticatorData();
            }
            // Verify that clientData commits to the expected client challenge
            // Use the Base64Url encoding which omits padding characters to match Webauthn Specification
            string memory challengeEncoded = Base64Url.encode(abi.encodePacked(clientChallenge));
            console2.logString("challengeEncoded: ");
            console2.logString(challengeEncoded);
            bytes memory challengeExtracted = new bytes(bytes(challengeEncoded).length);

            assembly {
                calldatacopy(
                    add(challengeExtracted, 32),
                    add(clientData.offset, clientChallengeDataOffset),
                    mload(challengeExtracted)
                )
            }
            // logs to compare the bytes retrieved against the bytes expected
            console2.logString("bytes(challengeEncoded): ");
            console2.logBytes(bytes(challengeEncoded));

            console2.logString("challengeExtracted: ");
            console2.logBytes(challengeExtracted);

            bytes32 moreData; //=keccak256(abi.encodePacked(challengeExtracted));
            assembly {
                moreData := keccak256(add(challengeExtracted, 32), mload(challengeExtracted))
            }
            if (keccak256(abi.encodePacked(bytes(challengeEncoded))) != moreData) {
                revert InvalidClientData();
            }
        } //avoid stack full

        // Verify the signature over sha256(authenticatorData || sha256(clientData))
        bytes memory verifyData = new bytes(authenticatorData.length + 32);

        assembly {
            calldatacopy(add(verifyData, 32), authenticatorData.offset, authenticatorData.length)
        }

        bytes32 more = sha256(clientData);
        assembly {
            mstore(add(verifyData, add(authenticatorData.length, 32)), more)
        }
        console2.logBytes32(sha256(verifyData));

        return sha256(verifyData);
    }
}
