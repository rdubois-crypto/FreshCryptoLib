import type { BigNumber } from "ethers";
/**
 * Try to obtain the return data of a transaction from the given value.
 *
 * If the value is an error but it doesn't have data, we assume it's not related
 * to a reverted transaction and we re-throw it.
 */
export declare function getReturnDataFromError(error: any): string;
declare type DecodedReturnData = {
    kind: "Error";
    reason: string;
} | {
    kind: "Empty";
} | {
    kind: "Panic";
    code: BigNumber;
    description: string;
} | {
    kind: "Custom";
    id: string;
    data: string;
};
export declare function decodeReturnData(returnData: string): DecodedReturnData;
export {};
//# sourceMappingURL=utils.d.ts.map