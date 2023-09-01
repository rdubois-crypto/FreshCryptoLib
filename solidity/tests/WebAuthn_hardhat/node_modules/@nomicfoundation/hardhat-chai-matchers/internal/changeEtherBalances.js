"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.getBalanceChanges = exports.supportChangeEtherBalances = void 0;
const ordinal_1 = __importDefault(require("ordinal"));
const utils_1 = require("../utils");
const account_1 = require("./misc/account");
const balance_1 = require("./misc/balance");
function supportChangeEtherBalances(Assertion) {
    Assertion.addMethod("changeEtherBalances", function (accounts, balanceChanges, options) {
        const { BigNumber } = require("ethers");
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        let subject = this._obj;
        if (typeof subject === "function") {
            subject = subject();
        }
        const checkBalanceChanges = ([actualChanges, accountAddresses]) => {
            const assert = (0, utils_1.buildAssert)(negated, checkBalanceChanges);
            assert(actualChanges.every((change, ind) => change.eq(BigNumber.from(balanceChanges[ind]))), () => {
                const lines = [];
                actualChanges.forEach((change, i) => {
                    if (!change.eq(BigNumber.from(balanceChanges[i]))) {
                        lines.push(`Expected the ether balance of ${accountAddresses[i]} (the ${(0, ordinal_1.default)(i + 1)} address in the list) to change by ${balanceChanges[i].toString()} wei, but it changed by ${change.toString()} wei`);
                    }
                });
                return lines.join("\n");
            }, () => {
                const lines = [];
                actualChanges.forEach((change, i) => {
                    if (change.eq(BigNumber.from(balanceChanges[i]))) {
                        lines.push(`Expected the ether balance of ${accountAddresses[i]} (the ${(0, ordinal_1.default)(i + 1)} address in the list) NOT to change by ${balanceChanges[i].toString()} wei, but it did`);
                    }
                });
                return lines.join("\n");
            });
        };
        const derivedPromise = Promise.all([
            getBalanceChanges(subject, accounts, options),
            (0, balance_1.getAddresses)(accounts),
        ]).then(checkBalanceChanges);
        this.then = derivedPromise.then.bind(derivedPromise);
        this.catch = derivedPromise.catch.bind(derivedPromise);
        this.promise = derivedPromise;
        return this;
    });
}
exports.supportChangeEtherBalances = supportChangeEtherBalances;
async function getBalanceChanges(transaction, accounts, options) {
    const txResponse = await transaction;
    const txReceipt = await txResponse.wait();
    const txBlockNumber = txReceipt.blockNumber;
    const balancesAfter = await (0, balance_1.getBalances)(accounts, txBlockNumber);
    const balancesBefore = await (0, balance_1.getBalances)(accounts, txBlockNumber - 1);
    const txFees = await getTxFees(accounts, txResponse, options);
    return balancesAfter.map((balance, ind) => balance.add(txFees[ind]).sub(balancesBefore[ind]));
}
exports.getBalanceChanges = getBalanceChanges;
async function getTxFees(accounts, txResponse, options) {
    return Promise.all(accounts.map(async (account) => {
        if (options?.includeFee !== true &&
            (await (0, account_1.getAddressOf)(account)) === txResponse.from) {
            const txReceipt = await txResponse.wait();
            const gasPrice = txReceipt.effectiveGasPrice ?? txResponse.gasPrice;
            const gasUsed = txReceipt.gasUsed;
            const txFee = gasPrice.mul(gasUsed);
            return txFee;
        }
        return 0;
    }));
}
//# sourceMappingURL=changeEtherBalances.js.map