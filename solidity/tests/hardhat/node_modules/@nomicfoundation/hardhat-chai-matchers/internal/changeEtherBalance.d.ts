/// <reference types="chai" />
import type { providers } from "ethers";
import { Account } from "./misc/account";
import { BalanceChangeOptions } from "./misc/balance";
export declare function supportChangeEtherBalance(Assertion: Chai.AssertionStatic): void;
export declare function getBalanceChange(transaction: providers.TransactionResponse | Promise<providers.TransactionResponse> | (() => Promise<providers.TransactionResponse> | providers.TransactionResponse), account: Account | string, options?: BalanceChangeOptions): Promise<import("ethers").BigNumber>;
//# sourceMappingURL=changeEtherBalance.d.ts.map