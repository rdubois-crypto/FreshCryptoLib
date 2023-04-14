"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.supportRevertedWithoutReason = void 0;
const utils_1 = require("../../utils");
const utils_2 = require("./utils");
function supportRevertedWithoutReason(Assertion) {
    Assertion.addMethod("revertedWithoutReason", function () {
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        const onSuccess = () => {
            const assert = (0, utils_1.buildAssert)(negated, onSuccess);
            assert(false, `Expected transaction to be reverted without a reason, but it didn't revert`);
        };
        const onError = (error) => {
            const assert = (0, utils_1.buildAssert)(negated, onError);
            const returnData = (0, utils_2.getReturnDataFromError)(error);
            const decodedReturnData = (0, utils_2.decodeReturnData)(returnData);
            if (decodedReturnData.kind === "Error") {
                assert(false, `Expected transaction to be reverted without a reason, but it reverted with reason '${decodedReturnData.reason}'`);
            }
            else if (decodedReturnData.kind === "Empty") {
                assert(true, undefined, "Expected transaction NOT to be reverted without a reason, but it was");
            }
            else if (decodedReturnData.kind === "Panic") {
                assert(false, `Expected transaction to be reverted without a reason, but it reverted with panic code ${decodedReturnData.code.toHexString()} (${decodedReturnData.description})`);
            }
            else if (decodedReturnData.kind === "Custom") {
                assert(false, `Expected transaction to be reverted without a reason, but it reverted with a custom error`);
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
exports.supportRevertedWithoutReason = supportRevertedWithoutReason;
//# sourceMappingURL=revertedWithoutReason.js.map