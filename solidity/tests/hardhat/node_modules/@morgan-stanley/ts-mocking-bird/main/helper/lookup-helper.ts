import { ConstructorFunction, FunctionCallLookup, IMocked, LookupType } from '../mock';

export function getLookup<T, C extends ConstructorFunction<T>, U extends LookupType>(
    mock: IMocked<T, C>,
    lookupType: U,
): FunctionCallLookup<T, C, U> {
    switch (lookupType) {
        case 'function':
            return mock.functionCallLookup as FunctionCallLookup<T, C, U>;
        case 'getter':
            return mock.getterCallLookup as FunctionCallLookup<T, C, U>;
        case 'setter':
            return mock.setterCallLookup as FunctionCallLookup<T, C, U>;
        case 'staticFunction':
            return mock.staticFunctionCallLookup as FunctionCallLookup<T, C, U>;
        case 'staticGetter':
            return mock.staticGetterCallLookup as FunctionCallLookup<T, C, U>;
        case 'staticSetter':
            return mock.staticSetterCallLookup as FunctionCallLookup<T, C, U>;
    }
    throw new Error(`Unknown lookup type: ${lookupType}`);
}
