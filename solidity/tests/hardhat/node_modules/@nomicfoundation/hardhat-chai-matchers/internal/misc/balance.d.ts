import { Account } from "./account";
export interface BalanceChangeOptions {
    includeFee?: boolean;
}
export declare function getAddresses(accounts: Array<Account | string>): Promise<string[]>;
export declare function getBalances(accounts: Array<Account | string>, blockNumber?: number): Promise<import("ethers").BigNumber[]>;
//# sourceMappingURL=balance.d.ts.map