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
exports.clearTokenDescriptionsCache = exports.getBalanceChange = exports.supportChangeTokenBalance = void 0;
const utils_1 = require("../utils");
const utils_2 = require("./calledOnContract/utils");
const account_1 = require("./misc/account");
function supportChangeTokenBalance(Assertion) {
    Assertion.addMethod("changeTokenBalance", function (token, account, balanceChange) {
        const ethers = require("ethers");
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        let subject = this._obj;
        if (typeof subject === "function") {
            subject = subject();
        }
        checkToken(token, "changeTokenBalance");
        const checkBalanceChange = ([actualChange, address, tokenDescription]) => {
            const assert = (0, utils_1.buildAssert)(negated, checkBalanceChange);
            assert(actualChange.eq(ethers.BigNumber.from(balanceChange)), `Expected the balance of ${tokenDescription} tokens for "${address}" to change by ${balanceChange.toString()}, but it changed by ${actualChange.toString()}`, `Expected the balance of ${tokenDescription} tokens for "${address}" NOT to change by ${balanceChange.toString()}, but it did`);
        };
        const derivedPromise = Promise.all([
            getBalanceChange(subject, token, account),
            (0, account_1.getAddressOf)(account),
            getTokenDescription(token),
        ]).then(checkBalanceChange);
        this.then = derivedPromise.then.bind(derivedPromise);
        this.catch = derivedPromise.catch.bind(derivedPromise);
        return this;
    });
    Assertion.addMethod("changeTokenBalances", function (token, accounts, balanceChanges) {
        const ethers = require("ethers");
        // capture negated flag before async code executes; see buildAssert's jsdoc
        const negated = this.__flags.negate;
        let subject = this._obj;
        if (typeof subject === "function") {
            subject = subject();
        }
        checkToken(token, "changeTokenBalances");
        if (accounts.length !== balanceChanges.length) {
            throw new Error(`The number of accounts (${accounts.length}) is different than the number of expected balance changes (${balanceChanges.length})`);
        }
        const balanceChangesPromise = Promise.all(accounts.map((account) => getBalanceChange(subject, token, account)));
        const addressesPromise = Promise.all(accounts.map(account_1.getAddressOf));
        const checkBalanceChanges = ([actualChanges, addresses, tokenDescription,]) => {
            const assert = (0, utils_1.buildAssert)(negated, checkBalanceChanges);
            assert(actualChanges.every((change, ind) => change.eq(ethers.BigNumber.from(balanceChanges[ind]))), `Expected the balances of ${tokenDescription} tokens for ${addresses} to change by ${balanceChanges}, respectively, but they changed by ${actualChanges}`, `Expected the balances of ${tokenDescription} tokens for ${addresses} NOT to change by ${balanceChanges}, respectively, but they did`);
        };
        const derivedPromise = Promise.all([
            balanceChangesPromise,
            addressesPromise,
            getTokenDescription(token),
        ]).then(checkBalanceChanges);
        this.then = derivedPromise.then.bind(derivedPromise);
        this.catch = derivedPromise.catch.bind(derivedPromise);
        return this;
    });
}
exports.supportChangeTokenBalance = supportChangeTokenBalance;
function checkToken(token, method) {
    if (typeof token !== "object" || token === null || !("functions" in token)) {
        throw new Error(`The first argument of ${method} must be the contract instance of the token`);
    }
    else if (token.functions.balanceOf === undefined) {
        throw new Error("The given contract instance is not an ERC20 token");
    }
}
async function getBalanceChange(transaction, token, account) {
    const ethers = require("ethers");
    const hre = await Promise.resolve().then(() => __importStar(require("hardhat")));
    const provider = hre.network.provider;
    const txResponse = await transaction;
    const txReceipt = await txResponse.wait();
    const txBlockNumber = txReceipt.blockNumber;
    const block = await provider.send("eth_getBlockByHash", [
        txReceipt.blockHash,
        false,
    ]);
    (0, utils_2.ensure)(block.transactions.length === 1, Error, "Multiple transactions found in block");
    const address = await (0, account_1.getAddressOf)(account);
    const balanceAfter = await token.balanceOf(address, {
        blockTag: txBlockNumber,
    });
    const balanceBefore = await token.balanceOf(address, {
        blockTag: txBlockNumber - 1,
    });
    return ethers.BigNumber.from(balanceAfter).sub(balanceBefore);
}
exports.getBalanceChange = getBalanceChange;
let tokenDescriptionsCache = {};
/**
 * Get a description for the given token. Use the symbol of the token if
 * possible; if it doesn't exist, the name is used; if the name doesn't
 * exist, the address of the token is used.
 */
async function getTokenDescription(token) {
    if (tokenDescriptionsCache[token.address] === undefined) {
        let tokenDescription = `<token at ${token.address}>`;
        try {
            tokenDescription = await token.symbol();
        }
        catch (e) {
            try {
                tokenDescription = await token.name();
            }
            catch (e2) { }
        }
        tokenDescriptionsCache[token.address] = tokenDescription;
    }
    return tokenDescriptionsCache[token.address];
}
// only used by tests
function clearTokenDescriptionsCache() {
    tokenDescriptionsCache = {};
}
exports.clearTokenDescriptionsCache = clearTokenDescriptionsCache;
//# sourceMappingURL=changeTokenBalance.js.map