interface EthersBigNumberLike {
    toHexString(): string;
}
interface BNLike {
    toNumber(): number;
    toString(base?: number): string;
}
export declare type NumberLike = number | bigint | string | EthersBigNumberLike | BNLike;
export declare type BlockTag = "latest" | "earliest" | "pending";
export {};
//# sourceMappingURL=types.d.ts.map