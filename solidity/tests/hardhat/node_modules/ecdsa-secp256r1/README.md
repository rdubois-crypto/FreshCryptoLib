# ECDSA secp256r1

<a href="https://www.npmjs.com/package/ecdsa-secp256r1"><img alt="npm-status" src="https://img.shields.io/npm/v/ecdsa-secp256r1.svg?style=flat" /></a>
<a href="https://github.com/forevertz/ecdsa-secp256r1/blob/master/LICENSE"><img alt="license" src="https://img.shields.io/badge/license-MIT_License-blue.svg?style=flat" /></a>

## Getting Started

```shell
$ yarn add ecdsa-secp256r1
```

## How to use

### Node.js

#### Create key

```javascript
const ECDSA = require('ecdsa-secp256r1')

const privateKey = ECDSA.generateKey()

privateKey.toJWK()
/*
{ kty: 'EC',
  crv: 'P-256',
  x: '4YdUIhIDncVu5tScgjxthiXOO_el11FWb56gR3qnhVQ',
  y: 'UyEvWOJbMZa9PtggGeRC9iQcAzOZZsyXpFE1qaF6jFk',
  d: 'TYVI2fW-nHSPGCx0MhWasg2Ggiyl1E_Kq4D1A5LmkxU' }
*/

privateKey.asPublic().toJWK()
/*
{ kty: 'EC',
  crv: 'P-256',
  x: '4YdUIhIDncVu5tScgjxthiXOO_el11FWb56gR3qnhVQ',
  y: 'UyEvWOJbMZa9PtggGeRC9iQcAzOZZsyXpFE1qaF6jFk' }
*/

privateKey.toPEM()
/*
-----BEGIN PRIVATE KEY-----
MIGHAgEAMBMGByqGSM49AgEGCCqGSM49AwEHBG0wawIBAQQgTYVI2fW+nHSPGCx0
MhWasg2Ggiyl1E/Kq4D1A5LmkxWhRANCAAThh1QiEgOdxW7m1JyCPG2GJc4796XX
UVZvnqBHeqeFVFMhL1jiWzGWvT7YIBnkQvYkHAMzmWbMl6RRNamheoxZ
-----END PRIVATE KEY-----
*/

privateKey.asPublic().toPEM()
/*
-----BEGIN PUBLIC KEY-----
MFkwEwYHKoZIzj0CAQYIKoZIzj0DAQcDQgAE4YdUIhIDncVu5tScgjxthiXOO/el
11FWb56gR3qnhVRTIS9Y4lsxlr0+2CAZ5EL2JBwDM5lmzJekUTWpoXqMWQ==
-----END PUBLIC KEY-----
*/

privateKey.toCompressedPublicKey()
/*
A+GHVCISA53FbubUnII8bYYlzjv3pddRVm+eoEd6p4VU
*/
```

#### Sign

```javascript
const message = { text: 'hello' }
privateKey.sign(JSON.stringify(message))
/*
lY3Lf9xDtcsqom5IKu+ZyikxeYHlEuxnPfme4lMxp76NMkIm5BiLxVjbqBSo4itfT/LEuBCzMXl11cB0w/X8dA==
*/
```

#### Verify

```javascript
const key = 'A+GHVCISA53FbubUnII8bYYlzjv3pddRVm+eoEd6p4VU'
const publicKey = ECDSA.fromCompressedPublicKey(key) // or ECDSA.fromJWK
const message = { text: 'hello' }
const signature = 'lY3Lf9xDtcsqom5IKu+ZyikxeYHlEuxnPfme4lMxp76NMkIm5BiLxVjbqBSo4itfT/LEuBCzMXl11cB0w/X8dA=='

publicKey.verify(JSON.stringify(message), signature)
/*
true
*/
```

### Browsers

Support: https://caniuse.com/#feat=cryptography

```javascript
const ECDSA = require('ecdsa-secp256r1/browser')

const privateKey = await ECDSA.generateKey()

// API is the same as Node.js, except everything returns Promises
```

or

```html
<script src="/PATH/TO/browser.js"></script>
<script>
;(async function() {
  const privateKey = await ECDSA.generateKey()
  const message = { text: 'hello' }
  const signature = await privateKey.sign(JSON.stringify(message))
})()
</script>
```

### Thanks

Inspired from https://github.com/relocately/ec-key
