"use strict";
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    Object.defineProperty(o, k2, { enumerable: true, get: function() { return m[k]; } });
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (k !== "default" && Object.prototype.hasOwnProperty.call(mod, k)) __createBinding(result, mod, k);
    __setModuleDefault(result, mod);
    return result;
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.getBalanceChange = exports.supportChangeEtherBalance = void 0;
const utils_1 = require("../utils");
const utils_2 = require("./calledOnContract/utils");
const account_1 = require("./misc/account");
function supportChangeEtherBalance(Assertion) {
    Assertion.addMethod("changeEtherBalance", function (account, balanceChange, options) {
        const { BigNumber } = require("ethers");
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        const subject = this._obj;
        const checkBalanceChange = ([actualChange, address]) => {
            const assert = (0, utils_1.buildAssert)(negated, checkBalanceChange);
            assert(actualChange.eq(BigNumber.from(balanceChange)), `Expected the ether balance of "${address}" to change by ${balanceChange.toString()} wei, but it changed by ${actualChange.toString()} wei`, `Expected the ether balance of "${address}" NOT to change by ${balanceChange.toString()} wei, but it did`);
        };
        const derivedPromise = Promise.all([
            getBalanceChange(subject, account, options),
            (0, account_1.getAddressOf)(account),
        ]).then(checkBalanceChange);
        this.then = derivedPromise.then.bind(derivedPromise);
        this.catch = derivedPromise.catch.bind(derivedPromise);
        this.promise = derivedPromise;
        return this;
    });
}
exports.supportChangeEtherBalance = supportChangeEtherBalance;
async function getBalanceChange(transaction, account, options) {
    const { BigNumber } = await Promise.resolve().then(() => __importStar(require("ethers")));
    const hre = await Promise.resolve().then(() => __importStar(require("hardhat")));
    const provider = hre.network.provider;
    let txResponse;
    if (typeof transaction === "function") {
        txResponse = await transaction();
    }
    else {
        txResponse = await transaction;
    }
    const txReceipt = await txResponse.wait();
    const txBlockNumber = txReceipt.blockNumber;
    const block = await provider.send("eth_getBlockByHash", [
        txReceipt.blockHash,
        false,
    ]);
    (0, utils_2.ensure)(block.transactions.length === 1, Error, "Multiple transactions found in block");
    const address = await (0, account_1.getAddressOf)(account);
    const balanceAfter = await provider.send("eth_getBalance", [
        address,
        `0x${txBlockNumber.toString(16)}`,
    ]);
    const balanceBefore = await provider.send("eth_getBalance", [
        address,
        `0x${(txBlockNumber - 1).toString(16)}`,
    ]);
    if (options?.includeFee !== true && address === txResponse.from) {
        const gasPrice = txReceipt.effectiveGasPrice ?? txResponse.gasPrice;
        const gasUsed = txReceipt.gasUsed;
        const txFee = gasPrice.mul(gasUsed);
        return BigNumber.from(balanceAfter).add(txFee).sub(balanceBefore);
    }
    else {
        return BigNumber.from(balanceAfter).sub(balanceBefore);
    }
}
exports.getBalanceChange = getBalanceChange;
//# sourceMappingURL=changeEtherBalance.js.map