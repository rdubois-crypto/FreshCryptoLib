"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.emitWithArgs = exports.supportEmit = exports.EMIT_CALLED = void 0;
const chai_1 = require("chai");
const util_1 = __importDefault(require("util"));
const ordinal_1 = __importDefault(require("ordinal"));
const utils_1 = require("../utils");
exports.EMIT_CALLED = "emitAssertionCalled";
async function waitForPendingTransaction(tx, provider) {
    let hash;
    if (tx instanceof Promise) {
        ({ hash } = await tx);
    }
    else if (typeof tx === "string") {
        hash = tx;
    }
    else {
        ({ hash } = tx);
    }
    if (hash === undefined) {
        throw new Error(`${JSON.stringify(tx)} is not a valid transaction`);
    }
    return provider.waitForTransaction(hash);
}
function supportEmit(Assertion, chaiUtils) {
    Assertion.addMethod("emit", function (contract, eventName) {
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        const tx = this._obj;
        const promise = this.then === undefined ? Promise.resolve() : this;
        const onSuccess = (receipt) => {
            const assert = (0, utils_1.buildAssert)(negated, onSuccess);
            let eventFragment;
            try {
                eventFragment = contract.interface.getEvent(eventName);
            }
            catch (e) {
                // ignore error
            }
            if (eventFragment === undefined) {
                throw new chai_1.AssertionError(`Event "${eventName}" doesn't exist in the contract`);
            }
            const topic = contract.interface.getEventTopic(eventFragment);
            this.logs = receipt.logs
                .filter((log) => log.topics.includes(topic))
                .filter((log) => log.address.toLowerCase() === contract.address.toLowerCase());
            assert(this.logs.length > 0, `Expected event "${eventName}" to be emitted, but it wasn't`, `Expected event "${eventName}" NOT to be emitted, but it was`);
            chaiUtils.flag(this, "eventName", eventName);
            chaiUtils.flag(this, "contract", contract);
        };
        const derivedPromise = promise
            .then(() => waitForPendingTransaction(tx, contract.provider))
            .then(onSuccess);
        chaiUtils.flag(this, exports.EMIT_CALLED, true);
        this.then = derivedPromise.then.bind(derivedPromise);
        this.catch = derivedPromise.catch.bind(derivedPromise);
        this.promise = derivedPromise;
        return this;
    });
}
exports.supportEmit = supportEmit;
async function emitWithArgs(context, Assertion, chaiUtils, expectedArgs, ssfi) {
    const negated = false; // .withArgs cannot be negated
    const assert = (0, utils_1.buildAssert)(negated, ssfi);
    tryAssertArgsArraysEqual(context, Assertion, chaiUtils, expectedArgs, context.logs, assert, ssfi);
}
exports.emitWithArgs = emitWithArgs;
function assertArgsArraysEqual(context, Assertion, chaiUtils, expectedArgs, log, assert, ssfi) {
    const { utils } = require("ethers");
    const actualArgs = chaiUtils.flag(context, "contract").interface.parseLog(log).args;
    const eventName = chaiUtils.flag(context, "eventName");
    assert(actualArgs.length === expectedArgs.length, `Expected "${eventName}" event to have ${expectedArgs.length} argument(s), but it has ${actualArgs.length}`);
    for (let index = 0; index < expectedArgs.length; index++) {
        if (typeof expectedArgs[index] === "function") {
            const errorPrefix = `The predicate for the ${(0, ordinal_1.default)(index + 1)} event argument`;
            try {
                assert(expectedArgs[index](actualArgs[index]), `${errorPrefix} returned false`
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
        else if (expectedArgs[index] instanceof Uint8Array) {
            new Assertion(actualArgs[index], undefined, ssfi, true).equal(utils.hexlify(expectedArgs[index]));
        }
        else if (expectedArgs[index]?.length !== undefined &&
            typeof expectedArgs[index] !== "string") {
            for (let j = 0; j < expectedArgs[index].length; j++) {
                new Assertion(actualArgs[index][j], undefined, ssfi, true).equal(expectedArgs[index][j]);
            }
        }
        else {
            if (actualArgs[index].hash !== undefined &&
                actualArgs[index]._isIndexed === true) {
                new Assertion(actualArgs[index].hash, undefined, ssfi, true).to.not.equal(expectedArgs[index], "The actual value was an indexed and hashed value of the event argument. The expected value provided to the assertion should be the actual event argument (the pre-image of the hash). You provided the hash itself. Please supply the the actual event argument (the pre-image of the hash) instead.");
                const expectedArgBytes = utils.isHexString(expectedArgs[index])
                    ? utils.arrayify(expectedArgs[index])
                    : utils.toUtf8Bytes(expectedArgs[index]);
                const expectedHash = utils.keccak256(expectedArgBytes);
                new Assertion(actualArgs[index].hash, undefined, ssfi, true).to.equal(expectedHash, `The actual value was an indexed and hashed value of the event argument. The expected value provided to the assertion was hashed to produce ${expectedHash}. The actual hash and the expected hash did not match`);
            }
            else {
                new Assertion(actualArgs[index], undefined, ssfi, true).equal(expectedArgs[index]);
            }
        }
    }
}
const tryAssertArgsArraysEqual = (context, Assertion, chaiUtils, expectedArgs, logs, assert, ssfi) => {
    if (logs.length === 1)
        return assertArgsArraysEqual(context, Assertion, chaiUtils, expectedArgs, logs[0], assert, ssfi);
    for (const index in logs) {
        if (index === undefined) {
            break;
        }
        else {
            try {
                assertArgsArraysEqual(context, Assertion, chaiUtils, expectedArgs, logs[index], assert, ssfi);
                return;
            }
            catch { }
        }
    }
    const eventName = chaiUtils.flag(context, "eventName");
    assert(false, `The specified arguments (${util_1.default.inspect(expectedArgs)}) were not included in any of the ${context.logs.length} emitted "${eventName}" events`);
};
//# sourceMappingURL=emit.js.map