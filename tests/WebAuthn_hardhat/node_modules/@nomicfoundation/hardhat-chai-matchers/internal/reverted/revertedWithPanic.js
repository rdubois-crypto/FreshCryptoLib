"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.supportRevertedWithPanic = void 0;
const common_1 = require("hardhat/common");
const utils_1 = require("../../utils");
const panic_1 = require("./panic");
const utils_2 = require("./utils");
function supportRevertedWithPanic(Assertion) {
    Assertion.addMethod("revertedWithPanic", function (expectedCodeArg) {
        const ethers = require("ethers");
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        let expectedCode;
        try {
            if (expectedCodeArg !== undefined) {
                const normalizedCode = (0, common_1.normalizeToBigInt)(expectedCodeArg);
                expectedCode = ethers.BigNumber.from(normalizedCode);
            }
        }
        catch {
            throw new TypeError(`Expected the given panic code to be a number-like value, but got '${expectedCodeArg}'`);
        }
        const code = expectedCode;
        let description;
        let formattedPanicCode;
        if (code === undefined) {
            formattedPanicCode = "some panic code";
        }
        else {
            const codeBN = ethers.BigNumber.from(code);
            description = (0, panic_1.panicErrorCodeToReason)(codeBN) ?? "unknown panic code";
            formattedPanicCode = `panic code ${codeBN.toHexString()} (${description})`;
        }
        const onSuccess = () => {
            const assert = (0, utils_1.buildAssert)(negated, onSuccess);
            assert(false, `Expected transaction to be reverted with ${formattedPanicCode}, but it didn't revert`);
        };
        const onError = (error) => {
            const assert = (0, utils_1.buildAssert)(negated, onError);
            const returnData = (0, utils_2.getReturnDataFromError)(error);
            const decodedReturnData = (0, utils_2.decodeReturnData)(returnData);
            if (decodedReturnData.kind === "Empty") {
                assert(false, `Expected transaction to be reverted with ${formattedPanicCode}, but it reverted without a reason`);
            }
            else if (decodedReturnData.kind === "Error") {
                assert(false, `Expected transaction to be reverted with ${formattedPanicCode}, but it reverted with reason '${decodedReturnData.reason}'`);
            }
            else if (decodedReturnData.kind === "Panic") {
                if (code !== undefined) {
                    assert(decodedReturnData.code.eq(code), `Expected transaction to be reverted with ${formattedPanicCode}, but it reverted with panic code ${decodedReturnData.code.toHexString()} (${decodedReturnData.description})`, `Expected transaction NOT to be reverted with ${formattedPanicCode}, but it was`);
                }
                else {
                    assert(true, undefined, `Expected transaction NOT to be reverted with ${formattedPanicCode}, but it reverted with panic code ${decodedReturnData.code.toHexString()} (${decodedReturnData.description})`);
                }
            }
            else if (decodedReturnData.kind === "Custom") {
                assert(false, `Expected transaction to be reverted with ${formattedPanicCode}, but it reverted with a custom error`);
            }
            else {
                const _exhaustiveCheck = decodedReturnData;
            }
        };
        const derivedPromise = Promise.resolve(this._obj).then(onSuccess, onError);
        this.then = derivedPromise.then.bind(derivedPromise);
        this.catch = derivedPromise.catch.bind(derivedPromise);
        return this;
    });
}
exports.supportRevertedWithPanic = supportRevertedWithPanic;
//# sourceMappingURL=revertedWithPanic.js.map