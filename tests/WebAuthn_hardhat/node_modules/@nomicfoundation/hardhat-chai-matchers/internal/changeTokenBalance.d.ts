/// <reference types="chai" />
import type EthersT from "ethers";
import { Account } from "./misc/account";
declare type TransactionResponse = EthersT.providers.TransactionResponse;
interface Token extends EthersT.Contract {
    balanceOf(address: string, overrides?: any): Promise<EthersT.BigNumber>;
}
export declare function supportChangeTokenBalance(Assertion: Chai.AssertionStatic): void;
export declare function getBalanceChange(transaction: TransactionResponse | Promise<TransactionResponse>, token: Token, account: Account | string): Promise<EthersT.ethers.BigNumber>;
export declare function clearTokenDescriptionsCache(): void;
export {};
//# sourceMappingURL=changeTokenBalance.d.ts.map