export type OperatorFunction<T, C extends ConstructorFunction<T>> = (value: IMocked<T, C>) => IMocked<T, C>;

export type MatchFunction<T> = (passedValue: T) => boolean;

// eslint-disable-next-line @typescript-eslint/ban-types
export type FunctionsOnly<T> = Pick<T, { [K in keyof T]: Required<T>[K] extends Function ? K : never }[keyof T]>;

// eslint-disable-next-line @typescript-eslint/ban-types
export type PropertiesOnly<T> = Pick<T, { [K in keyof T]: Required<T>[K] extends Function ? never : K }[keyof T]>;

/**
 * Allows custom logic to verify that a function parameter has the expected value.
 */
export interface IParameterMatcher<T = any> {
    /**
     * function that takes the actual parameter value and returns true if it was the expected value
     */
    readonly isExpectedValue: MatchFunction<T>;
    /**
     * A string representation of the expected value to be used in failure messages
     */
    readonly expectedDisplayValue: string;
    /**
     * Used to format a value into a string to display value of actual parameters passed in failure messages
     */
    readonly parameterToString?: (value: T) => string;
}

export type ParameterMatcher<T> = IParameterMatcher<T> | MatchFunction<T> | T;

type FunctionOrConstructor = (new (...params: any[]) => any) | ((...params: any[]) => any);

export type FunctionParameterMatchers<T extends any[]> = {
    [P in keyof T]: T[P] extends FunctionOrConstructor ? IParameterMatcher<T[P]> : ParameterMatcher<T[P]>;
};

export type FunctionCallLookup<T, C extends ConstructorFunction<T>, U extends LookupType> = {
    [P in FunctionName<T, C, U>]?: LookupParams<T, C, U, P>[];
};

export type LookupParams<
    T,
    C extends ConstructorFunction<T>,
    U extends LookupType,
    K extends FunctionName<T, C, U>,
> = U extends FunctionTypes
    ? FunctionParams<VerifierTarget<T, C, U>[K]>
    : U extends SetterTypes
    ? [VerifierTarget<T, C, U>[K]]
    : [];
export type FunctionParams<T> = T extends (...args: infer P) => any ? P : never;

export type ConstructorFunction<T> = (abstract new (...args: any[]) => T) | (new (...args: any[]) => T);

export type StaticLookupTypes = 'staticFunction' | 'staticGetter' | 'staticSetter';
export type InstanceLookupTypes = 'function' | 'getter' | 'setter';
export type SetterTypes = 'staticSetter' | 'setter';
export type GetterTypes = 'staticGetter' | 'getter';
export type FunctionTypes = 'staticFunction' | 'function';
export type LookupType = StaticLookupTypes | InstanceLookupTypes;

export type VerifierTarget<T, C extends ConstructorFunction<T>, U extends LookupType> = U extends StaticLookupTypes
    ? U extends FunctionTypes
        ? FunctionsOnly<C>
        : C
    : U extends FunctionTypes
    ? FunctionsOnly<T>
    : T;
export type FunctionName<T, C extends ConstructorFunction<T>, U extends LookupType> = keyof VerifierTarget<T, C, U>;

export interface IFunctionWithParametersVerification<
    P extends Array<any>,
    T,
    U extends LookupType,
    C extends ConstructorFunction<T> = never,
> extends IFunctionVerifier<T, U, C> {
    /**
     * Checks the parameters in a non-strict equality way.
     * defaults to the toEqual() matcher
     * Equivalent to expected == actual
     *
     * Expected parameter values can be passed that uses the default matcher (toEqual)
     * Other matchers, a comparison function or a custom IParameterMatcher can also be passed
     *
     * If your function accepts functions as parameters an IParameterMatcher must be used
     *
     * Note: if a matcher function ((value: T) => true) no information will be provided about what
     * value the parameter was expected to be in a test failure message.
     * To provide this expected value string implement IParameterMatcher instead
     *
     * @param args list of parameters to compare against
     */
    withParameters(...args: FunctionParameterMatchers<P>): IStrictFunctionVerification<T, U, C>;
    /**
     * Checks the parameters in a strict euqlity way.
     * defaults to the toBe() matcher
     * Equivalent to expected == actual
     *
     * Expected parameter values can be passed that uses the default matcher (toBe)
     * Other matchers, a comparison function or a custom IParameterMatcher can also be passed
     *
     * If your function accepts functions as parameters an IParameterMatcher must be used
     *
     * Note: if a matcher function ((value: T) => true) no information will be provided about what
     * value the parameter was expected to be in a test failure message.
     * To provide this expected value string implement IParameterMatcher instead
     *
     * @param args list of parameters to compare against
     */
    withParametersEqualTo(...args: FunctionParameterMatchers<P>): IStrictFunctionVerification<T, U, C>;
}

export interface IStrictFunctionVerification<T, U extends LookupType, C extends ConstructorFunction<T> = never>
    extends IFunctionVerifier<T, U, C> {
    /**
     * verify that the function has been called ONLY with the specified parameters and never without
     */
    strict(): IFunctionVerifier<T, U, C>;
}

export interface IFunctionVerifier<T, U extends LookupType, C extends ConstructorFunction<T> = never> {
    type: U;
    functionName: FunctionName<T, C, U>;
    parameterMatchers?: (MatchFunction<any> | IParameterMatcher<any>)[];
    strictCallCount: boolean;
    getMock(): IMocked<T, C>;
}

export interface IMocked<T, C extends ConstructorFunction<T> = never> {
    /**
     * The mocked object. This should be passed to your SUT.
     */
    mock: T;
    /**
     * The mocked constructor (if statics have been mocked).
     * This can be used to return the mocked instance with new myMock.constructor() and for accessing statics.
     */
    mockConstructor: C;

    functionCallLookup: FunctionCallLookup<T, C, 'function'>;
    setterCallLookup: FunctionCallLookup<T, C, 'setter'>;
    getterCallLookup: FunctionCallLookup<T, C, 'getter'>;

    staticFunctionCallLookup: FunctionCallLookup<T, C, 'staticFunction'>;
    staticSetterCallLookup: FunctionCallLookup<T, C, 'staticSetter'>;
    staticGetterCallLookup: FunctionCallLookup<T, C, 'staticGetter'>;

    functionReplacementLookup: Partial<Record<LookupType, Record<string, ((...args: any[]) => any) | undefined>>>;

    /**
     * Used to setup the mock with multiple operators.
     *
     * Mock.create<IMyService>().setup(
     *     setupFunction("functionName"),
     *     setupFunction("otherFunction"),
     *     setupProperty("propertyName"),
     * );
     *
     * @param operators
     */
    setup(...operators: OperatorFunction<T, C>[]): IMocked<T, C>;

    /**
     * Sets up a single function and returns a function verifier to verify calls made and parameters passed.
     *
     * @param functionName
     * @param mockFunction
     */
    setupFunction<K extends keyof FunctionsOnly<T>>(
        functionName: K,
        mockFunction?: T[K],
    ): IFunctionWithParametersVerification<FunctionParams<T[K]>, T, 'function', C>;
    /**
     * Sets up a single property and returns a function verifier to verify value get or set operations.
     *
     * @param propertyname
     * @param value
     */
    setupProperty<K extends keyof T>(propertyname: K, value?: T[K]): IFunctionVerifier<T, 'getter', C>;
    /**
     * Defines a single property and allows getters and setters to be defined.
     * Returns a function verifier to verify get and set operations
     *
     * @param propertyname
     * @param getter
     * @param setter
     */
    defineProperty<K extends keyof T>(
        propertyname: K,
        getter?: () => T[K],
        setter?: (value: T[K]) => void,
    ): IFunctionVerifier<T, 'getter', C>;

    /**
     * Sets up a single static function and returns a function verifier to verify calls made and parameters passed.
     *
     * @param functionName
     * @param mockFunction
     */
    setupStaticFunction<K extends keyof FunctionsOnly<C>>(
        functionName: K,
        mockFunction?: C[K],
    ): IFunctionWithParametersVerification<FunctionParams<C[K]>, T, 'staticFunction', C>;
    /**
     * Sets up a single static property and returns a function verifier to verify value get or set operations.
     *
     * @param propertyname
     * @param value
     */
    setupStaticProperty<K extends keyof C>(propertyname: K, value?: C[K]): IFunctionVerifier<T, 'staticGetter', C>;
    /**
     * Defines a single static property and allows getters and setters to be defined.
     * Returns a function verifier to verify get and set operations
     *
     * @param propertyname
     * @param getter
     * @param setter
     */
    defineStaticProperty<K extends keyof C>(
        propertyname: K,
        getter?: () => C[K],
        setter?: (value: C[K]) => void,
    ): IFunctionVerifier<T, 'staticGetter', C>;

    /**
     * Verifies calls to a previously setup function.
     * expect(myMock.withFunction("functionName")).wasNotCalled():
     * expect(myMock.withFunction("functionName")).wasCalledOnce():
     * expect(myMock.withFunction("functionName").withParameters("one", 2)).wasCalledOnce():
     *
     * @param functionName
     */
    withFunction<K extends keyof FunctionsOnly<T>>(
        functionName: K,
    ): IFunctionWithParametersVerification<FunctionParams<T[K]>, T, 'function', C>;
    /**
     * Verifies calls to a previously setup getter.
     * expect(myMock.withGetter("propertyName")).wasNotCalled():
     * expect(myMock.withGetter("propertyName")).wasCalledOnce():
     *
     * @param functionName
     */
    withGetter<K extends keyof T>(propertyname: K): IFunctionVerifier<T, 'getter', C>;
    /**
     * Verifies calls to a previously setup setter.
     * expect(myMock.withSetter("propertyName")).wasNotCalled():
     * expect(myMock.withSetter("propertyName")).wasCalledOnce():
     * expect(myMock.withSetter("propertyName").withParameters("one")).wasCalledOnce():
     *
     * @param functionName
     */
    withSetter<K extends keyof T>(propertyname: K): IFunctionWithParametersVerification<[T[K]], T, 'setter', C>;

    /**
     * Verifies calls to a previously setup static function.
     * expect(myMock.withStaticFunction("functionName")).wasNotCalled():
     * expect(myMock.withStaticFunction("functionName")).wasCalledOnce():
     * expect(myMock.withStaticFunction("functionName").withParameters("one", 2)).wasCalledOnce():
     *
     * @param functionName
     */
    withStaticFunction<K extends keyof FunctionsOnly<C>>(
        functionName: K,
    ): IFunctionWithParametersVerification<FunctionParams<C[K]>, T, 'staticFunction', C>;
    /**
     * Verifies calls to a previously setup static getter.
     * expect(myMock.withStaticGetter("functionName")).wasNotCalled():
     * expect(myMock.withStaticGetter("functionName")).wasCalledOnce():
     *
     * @param functionName
     */
    withStaticGetter<K extends keyof C>(propertyname: K): IFunctionVerifier<T, 'staticGetter', C>;
    /**
     * Verifies calls to a previously setup static setter.
     * expect(myMock.withStaticSetter("functionName")).wasNotCalled():
     * expect(myMock.withStaticSetter("functionName")).wasCalledOnce():
     * expect(myMock.withStaticSetter("functionName").withParameters("one")).wasCalledOnce():
     *
     * @param functionName
     */
    withStaticSetter<K extends keyof C>(
        propertyname: K,
    ): IFunctionWithParametersVerification<[C[K]], T, 'staticSetter', C>;
}

/**
 * Copied here from Jest types to avoid the need for consuming projects to install Jest types
 */
export interface IJestCustomMatcherResult {
    pass: boolean;
    message: () => string;
}

/**
 * Copied here from Jasmine types to avoid the need for consuming projects to install Jasmine types
 */
export interface IJasmineCustomMatcherResult {
    pass: boolean;
    message?: string;
}

export interface ICustomMatcher {
    compare<T>(actual: T, expected: T, ...args: any[]): IJasmineCustomMatcherResult;
    compare(actual: any, ...expected: any[]): IJasmineCustomMatcherResult;
    negativeCompare?<T>(actual: T, expected: T, ...args: any[]): IJasmineCustomMatcherResult;
    negativeCompare?(actual: any, ...expected: any[]): IJasmineCustomMatcherResult;
}
