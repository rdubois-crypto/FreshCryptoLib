"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
};
Object.defineProperty(exports, "__esModule", { value: true });
exports.getAddressOf = exports.isAccount = void 0;
const assert_1 = __importDefault(require("assert"));
function isAccount(account) {
    const ethers = require("ethers");
    return account instanceof ethers.Contract || account instanceof ethers.Wallet;
}
exports.isAccount = isAccount;
async function getAddressOf(account) {
    if (typeof account === "string") {
        (0, assert_1.default)(/^0x[0-9a-fA-F]{40}$/.test(account), `Invalid address ${account}`);
        return account;
    }
    else if (isAccount(account)) {
        return account.address;
    }
    else {
        return account.getAddress();
    }
}
exports.getAddressOf = getAddressOf;
//# sourceMappingURL=account.js.map