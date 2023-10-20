"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.revertedWithCustomErrorWithArgs = exports.supportRevertedWithCustomError = exports.REVERTED_WITH_CUSTOM_ERROR_CALLED = void 0;
const chai_1 = require("chai");
const utils_1 = require("../../utils");
const utils_2 = require("./utils");
exports.REVERTED_WITH_CUSTOM_ERROR_CALLED = "customErrorAssertionCalled";
function supportRevertedWithCustomError(Assertion, utils) {
    Assertion.addMethod("revertedWithCustomError", function (contract, expectedCustomErrorName) {
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        // check the case where users forget to pass the contract as the first
        // argument
        if (typeof contract === "string" || contract?.interface === undefined) {
            throw new TypeError("The first argument of .revertedWithCustomError must be the contract that defines the custom error");
        }
        // validate custom error name
        if (typeof expectedCustomErrorName !== "string") {
            throw new TypeError("Expected the custom error name to be a string");
        }
        const iface = contract.interface;
        const expectedCustomError = findCustomErrorByName(iface, expectedCustomErrorName);
        // check that interface contains the given custom error
        if (expectedCustomError === undefined) {
            throw new Error(`The given contract doesn't have a custom error named '${expectedCustomErrorName}'`);
        }
        const onSuccess = () => {
            const assert = (0, utils_1.buildAssert)(negated, onSuccess);
            assert(false, `Expected transaction to be reverted with custom error '${expectedCustomErrorName}', but it didn't revert`);
        };
        const onError = (error) => {
            const assert = (0, utils_1.buildAssert)(negated, onError);
            const returnData = (0, utils_2.getReturnDataFromError)(error);
            const decodedReturnData = (0, utils_2.decodeReturnData)(returnData);
            if (decodedReturnData.kind === "Empty") {
                assert(false, `Expected transaction to be reverted with custom error '${expectedCustomErrorName}', but it reverted without a reason`);
            }
            else if (decodedReturnData.kind === "Error") {
                assert(false, `Expected transaction to be reverted with custom error '${expectedCustomErrorName}', but it reverted with reason '${decodedReturnData.reason}'`);
            }
            else if (decodedReturnData.kind === "Panic") {
                assert(false, `Expected transaction to be reverted with custom error '${expectedCustomErrorName}', but it reverted with panic code ${decodedReturnData.code.toHexString()} (${decodedReturnData.description})`);
            }
            else if (decodedReturnData.kind === "Custom") {
                if (decodedReturnData.id === expectedCustomError.id) {
                    // add flag with the data needed for .withArgs
                    const customErrorAssertionData = {
                        contractInterface: iface,
                        customError: expectedCustomError,
                        returnData,
                    };
                    this.customErrorData = customErrorAssertionData;
                    assert(true, undefined, `Expected transaction NOT to be reverted with custom error '${expectedCustomErrorName}', but it was`);
                }
                else {
                    // try to decode the actual custom error
                    // this will only work when the error comes from the given contract
                    const actualCustomError = findCustomErrorById(iface, decodedReturnData.id);
                    if (actualCustomError === undefined) {
                        assert(false, `Expected transaction to be reverted with custom error '${expectedCustomErrorName}', but it reverted with a different custom error`);
                    }
                    else {
                        assert(false, `Expected transaction to be reverted with custom error '${expectedCustomErrorName}', but it reverted with custom error '${actualCustomError.name}'`);
                    }
                }
            }
            else {
                const _exhaustiveCheck = decodedReturnData;
            }
        };
        const derivedPromise = Promise.resolve(this._obj).then(onSuccess, onError);
        // needed for .withArgs
        utils.flag(this, exports.REVERTED_WITH_CUSTOM_ERROR_CALLED, true);
        this.promise = derivedPromise;
        this.then = derivedPromise.then.bind(derivedPromise);
        this.catch = derivedPromise.catch.bind(derivedPromise);
        return this;
    });
}
exports.supportRevertedWithCustomError = supportRevertedWithCustomError;
async function revertedWithCustomErrorWithArgs(context, Assertion, utils, expectedArgs, ssfi) {
    const negated = false; // .withArgs cannot be negated
    const assert = (0, utils_1.buildAssert)(negated, ssfi);
    const customErrorAssertionData = context.customErrorData;
    if (customErrorAssertionData === undefined) {
        throw new Error("[.withArgs] should never happen, please submit an issue to the Hardhat repository");
    }
    const { contractInterface, customError, returnData } = customErrorAssertionData;
    const errorFragment = contractInterface.errors[customError.signature];
    // We transform ether's Array-like object into an actual array as it's safer
    const actualArgs = Array.from(contractInterface.decodeErrorResult(errorFragment, returnData));
    new Assertion(actualArgs).to.have.same.length(expectedArgs.length, `expected ${expectedArgs.length} args but got ${actualArgs.length}`);
    for (const [i, actualArg] of actualArgs.entries()) {
        const expectedArg = expectedArgs[i];
        if (typeof expectedArg === "function") {
            const errorPrefix = `The predicate for custom error argument with index ${i}`;
            try {
                assert(expectedArg(actualArg), `${errorPrefix} returned false`
                // no need for a negated message, since we disallow mixing .not. with
                // .withArgs
                );
            }
            catch (e) {
                if (e instanceof chai_1.AssertionError) {
                    assert(false, `${errorPrefix} threw an AssertionError: ${e.message}`
                    // no need for a negated message, since we disallow mixing .not. with
                    // .withArgs
                    );
                }
                throw e;
            }
        }
        else if (Array.isArray(expectedArg)) {
            new Assertion(actualArg).to.deep.equal(expectedArg);
        }
        else {
            new Assertion(actualArg).to.equal(expectedArg);
        }
    }
}
exports.revertedWithCustomErrorWithArgs = revertedWithCustomErrorWithArgs;
function findCustomErrorByName(iface, name) {
    const ethers = require("ethers");
    const customErrorEntry = Object.entries(iface.errors).find(([, fragment]) => fragment.name === name);
    if (customErrorEntry === undefined) {
        return undefined;
    }
    const [customErrorSignature] = customErrorEntry;
    const customErrorId = ethers.utils.id(customErrorSignature).slice(0, 10);
    return {
        id: customErrorId,
        name,
        signature: customErrorSignature,
    };
}
function findCustomErrorById(iface, id) {
    const ethers = require("ethers");
    const customErrorEntry = Object.entries(iface.errors).find(([signature]) => ethers.utils.id(signature).slice(0, 10) === id);
    if (customErrorEntry === undefined) {
        return undefined;
    }
    return {
        id,
        name: customErrorEntry[1].name,
        signature: customErrorEntry[0],
    };
}
//# sourceMappingURL=revertedWithCustomError.js.map