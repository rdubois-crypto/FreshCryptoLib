import { isEqual } from 'lodash';
import { IParameterMatcher } from './contracts';

export function isParameterMatcher(matcher: unknown): matcher is IParameterMatcher<any> {
    const typedMatcher = matcher as IParameterMatcher<any>;

    return (
        matcher != null &&
        typeof typedMatcher.expectedDisplayValue === 'string' &&
        typeof typedMatcher.isExpectedValue === 'function'
    );
}

/**
 * Checks that the value is not undefined.
 * This will pass if the value is null (like jasmine does).
 * actualValue !== undefined
 */
export function toBeDefined(): IParameterMatcher<any> {
    return {
        isExpectedValue: (actualValue) => actualValue !== undefined,
        expectedDisplayValue: '<mustBeDefined>',
    };
}

/**
 * Checks that the value is not undefined and not null.
 * actualValue != null
 */
export function hasValue(): IParameterMatcher<any> {
    return {
        isExpectedValue: (actualValue) => actualValue != null,
        expectedDisplayValue: '<hasValue>',
    };
}

/**
 * Compares the expected value to the actual value with a strict equality check.
 * actualValue === expectedValue
 *
 * @param expectedValue
 */
export function toBe(expectedValue: any): IParameterMatcher<any> {
    return {
        expectedDisplayValue: mapItemToString(expectedValue),
        isExpectedValue: (actualValue) => actualValue === expectedValue,
    };
}

/**
 * Checks that the expected value is equal to the actual value using deep object comparison.
 * actualValue => isEqual(actualValue, expectedValue)
 *
 * @param expectedValue
 */
export function toEqual(expectedValue: any): IParameterMatcher<any> {
    return {
        expectedDisplayValue: mapItemToString(expectedValue),
        isExpectedValue: (actualValue) => isEqual(actualValue, expectedValue),
    };
}

/**
 * Allows any value.
 */
export function any(): IParameterMatcher<any> {
    return {
        expectedDisplayValue: `<matchAny>`,
        isExpectedValue: () => true,
    };
}

/**
 * Returns a string representation of a value.
 * @param item
 */
export function mapItemToString(item: any): string {
    if (typeof item === 'string') {
        return `"${item}"`;
    }

    if (item === undefined) {
        return 'undefined';
    }

    return JSON.stringify(item, replaceValue).replace(ReplaceUndefinedRegExp, 'undefined');
}

// https://regex101.com/r/BIvQJG/2
const functionStringRegExp = /function\s*\([^)]*\)/;

/**
 * string value that is VERY unlikely to be a genuine value
 */
const UndefinedReplacementString = `___@morgan-stanley/ts-mocking-bird_undefined_@morgan-stanley/ts-mocking-bird___`;

const ReplaceUndefinedRegExp = new RegExp(`"${UndefinedReplacementString}"`, 'g');

function replaceValue(_key: any, value: any) {
    switch (typeof value) {
        case 'function':
            // eslint-disable-next-line no-case-declarations
            const functionString = String(value);
            // eslint-disable-next-line no-case-declarations
            const regexpResult = functionStringRegExp.exec(functionString);

            return regexpResult != null ? regexpResult[0] : 'FUNCTION_BODY_REMOVED';

        case 'undefined':
            return UndefinedReplacementString;
    }

    return value;
}
