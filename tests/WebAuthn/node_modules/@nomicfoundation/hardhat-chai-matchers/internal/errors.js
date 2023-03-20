"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.HardhatChaiMatchersDecodingError = void 0;
const common_1 = require("hardhat/common");
class HardhatChaiMatchersDecodingError extends common_1.CustomError {
    constructor(encodedData, type, parent) {
        const message = `There was an error decoding '${encodedData}' as a ${type}`;
        super(message, parent);
    }
}
exports.HardhatChaiMatchersDecodingError = HardhatChaiMatchersDecodingError;
//# sourceMappingURL=errors.js.map