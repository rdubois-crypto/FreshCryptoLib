import { ConstructorFunction, FunctionsOnly, IMocked, OperatorFunction } from './contracts';
import { addMatchers } from './matchers';
import {
    defineProperty,
    defineStaticProperty,
    setupFunction,
    setupProperty,
    setupStaticFunction,
    setupStaticProperty,
} from './operators';
import { createFunctionParameterVerifier, createFunctionVerifier } from './verifiers';

export class Mock {
    public static create<T, C extends ConstructorFunction<T> = never>(): IMocked<T, C> {
        addMatchers();
        const mocked: IMocked<T, C> = {
            functionCallLookup: {},
            setterCallLookup: {},
            getterCallLookup: {},

            staticFunctionCallLookup: {},
            staticSetterCallLookup: {},
            staticGetterCallLookup: {},

            functionReplacementLookup: {},

            mock: {} as T,

            // eslint-disable-next-line @typescript-eslint/no-empty-function
            mockConstructor: ((..._args: any[]) => {}) as any,

            setup: (...operators: OperatorFunction<T, C>[]) => {
                let operatorMocked = mocked;
                operators.forEach((operator) => (operatorMocked = operator(mocked)));
                return operatorMocked;
            },

            setupFunction: <K extends keyof FunctionsOnly<T>>(functionName: K, mockFunction?: any) => {
                setupFunction<T, C, K>(functionName, mockFunction)(mocked);
                return mocked.withFunction(functionName);
            },
            setupProperty: <K extends keyof T>(propertyName: K, value?: T[K]) => {
                setupProperty<T, C, K>(propertyName, value)(mocked);
                return mocked.withGetter(propertyName);
            },
            defineProperty: <K extends keyof T>(
                propertyName: K,
                getter?: () => T[K],
                setter?: (value: T[K]) => void,
            ) => {
                defineProperty<T, C, K>(propertyName, getter, setter)(mocked);
                return mocked.withGetter(propertyName);
            },

            setupStaticFunction: <K extends keyof FunctionsOnly<C>>(functionName: K, mockFunction?: any) => {
                setupStaticFunction<T, C, K>(functionName, mockFunction)(mocked);
                return mocked.withStaticFunction(functionName);
            },
            setupStaticProperty: <K extends keyof C>(propertyName: K, value?: C[K]) => {
                setupStaticProperty<T, C, K>(propertyName, value)(mocked);
                return mocked.withStaticGetter(propertyName);
            },
            defineStaticProperty: <K extends keyof C>(
                propertyName: K,
                getter?: () => C[K],
                setter?: (value: C[K]) => void,
            ) => {
                defineStaticProperty<T, C, K>(propertyName, getter, setter)(mocked);
                return mocked.withStaticGetter(propertyName);
            },

            withFunction: <U extends keyof FunctionsOnly<T>>(functionName: U) =>
                createFunctionParameterVerifier(mocked, 'function', functionName),
            withSetter: <U extends keyof T>(functionName: U) =>
                createFunctionParameterVerifier(mocked, 'setter', functionName),
            withGetter: <U extends keyof T>(functionName: U) => createFunctionVerifier(mocked, 'getter', functionName),

            withStaticFunction: <U extends keyof FunctionsOnly<C>>(functionName: U) =>
                createFunctionParameterVerifier(mocked, 'staticFunction', functionName),
            withStaticSetter: <U extends keyof C>(functionName: U) =>
                createFunctionParameterVerifier(mocked, 'staticSetter', functionName),
            withStaticGetter: <U extends keyof C>(functionName: U) =>
                createFunctionVerifier(mocked, 'staticGetter', functionName),
        };

        mocked.mockConstructor.prototype = mocked.mock;

        return mocked;
    }
}
