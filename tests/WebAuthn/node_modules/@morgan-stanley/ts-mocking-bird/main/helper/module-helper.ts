import { v4 } from 'uuid';

export type WrappedModule = {
    ___moduleId: string;
};

export function isWrappedModule(value: any): value is WrappedModule {
    return value != null && typeof (value as WrappedModule).___moduleId === 'string';
}

const MockMap = new Map<string, any[] | undefined>();

/**
 * Takes an existing function and wraps all functions and constructors.
 * This allows us to replace the function or constructor after it has been imported into the system under test
 * This will not replace statics on classes or replace constant values in the module
 * @param originalModule
 */
export function proxyModule<T>(originalModule: T): WrappedModule & T {
    const wrappedModule = { ___moduleId: v4() } as WrappedModule & T;

    Object.keys(originalModule).forEach((key) => {
        const originalValue = (originalModule as any)[key];
        (wrappedModule as any)[key] =
            typeof originalValue === 'function' ? wrapFunction(wrappedModule, key, originalValue) : originalValue;
    });

    return wrappedModule;
}

/**
 * Replaces functions or constructors in a previously wrapped module.
 * If a module is not wrapped a warning is logged and no replacements occur
 * @param moduleProxy
 * @param mock
 */
export function registerMock<T>(moduleProxy: T, mock: Partial<T>) {
    if (!isWrappedModule(moduleProxy)) {
        console.warn(`Not registering mock as supplied module is not a wrapped module`);
        return;
    }
    const mocks = lookupMocksForModule(moduleProxy);

    mocks.push(mock);
}

/**
 * Replaces previously replaced module members
 * @param moduleProxy
 */
export function reset<T>(moduleProxy: T) {
    if (!isWrappedModule(moduleProxy)) {
        return;
    }

    MockMap.set(moduleProxy.___moduleId, undefined);
}

function lookupMocksForModule(originalModule: WrappedModule): any[] {
    let mocks = MockMap.get(originalModule.___moduleId);

    if (mocks === undefined) {
        mocks = [];
        MockMap.set(originalModule.___moduleId, mocks);
    }

    return mocks;
}

// eslint-disable-next-line @typescript-eslint/ban-types
function wrapFunction(this: any, mockModule: WrappedModule, key: string, func: Function) {
    function wrapped(this: any, ...args: any[]) {
        const returnFunc = getMocked(mockModule, key) || func;
        return returnFunc.apply(this, args);
    }

    wrapped.prototype = func.prototype;

    return wrapped;
}

function getMocked(mockModule: WrappedModule, key: string): ((...args: any[]) => any) | undefined {
    const mocks = lookupMocksForModule(mockModule);

    for (let index = 0; index < mocks.length; index++) {
        const mock = mocks[index];
        if (typeof mock[key] === 'function') {
            return mock[key];
        }
    }

    return undefined;
}
