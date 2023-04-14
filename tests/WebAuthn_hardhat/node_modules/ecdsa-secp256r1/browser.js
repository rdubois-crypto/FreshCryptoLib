/* global crypto, btoa, atob */

/* ========================================================================== *
 * From RFC-4492 (Appendix A) Equivalent Curves (Informative)                 *
 * ========================================================================== *
 *                                                                            *
 * +------------------------------------------------------------------------+ *
 * |                         Curve names chosen by                          | *
 * |                   different standards organizations                    | *
 * +-----------+------------+------------+------------+---------------------+ *
 * |   SECG    | ANSI X9.62 |    NIST    |  OpenSSL   |      ASN.1 OID      | *
 * +-----------+------------+------------+------------+---------------------+ *
 * | secp256r1 | prime256v1 | NIST P-256 | prime256v1 | 1.2.840.10045.3.1.7 | *
 * +-----------+------------+------------+------------+---------------------+ *
 * ========================================================================== */

;(function() {
  var ALGO = { name: 'ECDSA', namedCurve: 'P-256' }
  var SIGN_ALGO = { name: 'ECDSA', hash: { name: 'SHA-256' } }

  /* ========================================================================== *
   * CLASS DEFINITION                                                           *
   * ========================================================================== */

  function ECDSA(keys /* { publicKey: CryptoKey, privateKey?: CryptoKey } */) {
    if (keys.publicKey) this.publicKey = keys.publicKey
    if (keys.privateKey) this.privateKey = keys.privateKey
  }

  /* ========================================================================== *
   * UTILS                                                                      *
   * ========================================================================== */

  function toBase64(buffer) {
    return btoa(buffer)
  }

  function fromBase64(string) {
    return atob(string)
  }

  function arrayBufferToString(buffer) {
    return String.fromCharCode.apply(null, new Uint8Array(buffer))
  }

  function stringToArrayBuffer(string) {
    if (window.TextEncoder) { // Chrome, Firefox, Opera
      return new TextEncoder('utf-8').encode(string)
    } else {
      // TextEncoder polyfill (https://developer.mozilla.org/en-US/docs/Web/API/TextEncoder)
      var stringLength = string.length
      var buffer = new Uint8Array(stringLength * 3)
      var resPos = -1
      for (var point = 0, nextcode = 0, i = 0; i !== stringLength;) {
        ;(point = string.charCodeAt(i)), (i += 1)
        if (point >= 0xd800 && point <= 0xdbff) {
          if (i === stringLength) {
            buffer[(resPos += 1)] = 0xef /* 0b11101111 */
            buffer[(resPos += 1)] = 0xbf /* 0b10111111 */
            buffer[(resPos += 1)] = 0xbd /* 0b10111101 */
            break
          }
          nextcode = string.charCodeAt(i)
          if (nextcode >= 0xdc00 && nextcode <= 0xdfff) {
            point = (point - 0xd800) * 0x400 + nextcode - 0xdc00 + 0x10000
            i += 1
            if (point > 0xffff) {
              buffer[(resPos += 1)] = (0x1e /* 0b11110 */ << 3) | (point >>> 18)
              buffer[(resPos += 1)] = (0x2 /* 0b10 */ << 6) | ((point >>> 12) & 0x3f) /* 0b00111111 */
              buffer[(resPos += 1)] = (0x2 /* 0b10 */ << 6) | ((point >>> 6) & 0x3f) /* 0b00111111 */
              buffer[(resPos += 1)] = (0x2 /* 0b10 */ << 6) | (point & 0x3f) /* 0b00111111 */
              continue
            }
          } else {
            buffer[(resPos += 1)] = 0xef /* 0b11101111 */
            buffer[(resPos += 1)] = 0xbf /* 0b10111111 */
            buffer[(resPos += 1)] = 0xbd /* 0b10111101 */
            continue
          }
        }
        if (point <= 0x007f) {
          buffer[(resPos += 1)] = (0x0 /* 0b0 */ << 7) | point
        } else if (point <= 0x07ff) {
          buffer[(resPos += 1)] = (0x6 /* 0b110 */ << 5) | (point >>> 6)
          buffer[(resPos += 1)] = (0x2 /* 0b10 */ << 6) | (point & 0x3f) /* 0b00111111 */
        } else {
          buffer[(resPos += 1)] = (0xe /* 0b1110 */ << 4) | (point >>> 12)
          buffer[(resPos += 1)] = (0x2 /* 0b10 */ << 6) | ((point >>> 6) & 0x3f) /* 0b00111111 */
          buffer[(resPos += 1)] = (0x2 /* 0b10 */ << 6) | (point & 0x3f) /* 0b00111111 */
        }
      }
      buffer = new Uint8Array(buffer.buffer.slice(0, resPos + 1))
      return buffer
    }
  }

  function hash(object) {
    return new Promise(resolve => {
      var buffer = stringToArrayBuffer(typeof object === 'string' ? object : JSON.stringify(object))
      crypto.subtle.digest('SHA-256', buffer).then(sha256 => {
        resolve(toBase64(arrayBufferToString(sha256)))
      })
    })
  }

  /* ========================================================================== *
   * FACTORIES                                                                  *
   * ========================================================================== */

  ECDSA.generateKey = () => /*: Promise<ECDSA> */ {
    return new Promise(resolve => {
      crypto.subtle.generateKey(ALGO, true, ['sign', 'verify']).then((key) => {
        resolve(new ECDSA(key))
      })
    })
  }

  ECDSA.fromJWK = (jwk /*: Object */) => /*: Promise<ECDSA> */ {
    return new Promise(resolve => {
      var publicJwk = { kty: jwk.kty, crv: jwk.crv, x: jwk.x, y: jwk.y }
      var keyPromises = [
        crypto.subtle.importKey('jwk', publicJwk, ALGO, true, ['verify'])
      ]
      if (jwk.d) {
        var privateJwk = { kty: jwk.kty, crv: jwk.crv, d: jwk.d, x: jwk.x, y: jwk.y }
        keyPromises.push(crypto.subtle.importKey('jwk', privateJwk, ALGO, true, ['sign']))
      }
      Promise.all(keyPromises).then(keys => {
        resolve(new ECDSA({ publicKey: keys[0], privateKey: keys[1] }))
      })
    })
  }

  ECDSA.fromCompressedPublicKey = (base64Key /*: string */) => /*: Promise<ECDSA> */ {
    return new Promise(resolve => {
      var rawCompressedKey = stringToArrayBuffer(fromBase64(base64Key))
      crypto.subtle.importKey('raw', rawCompressedKey, ALGO, true, ['verify']).then((key) => {
        resolve(new ECDSA({ publicKey: key }))
      })
    })
  }

  ECDSA.fromBase64PrivateKey = (base64Key /*: string */) => /*: Promise<ECDSA> */ {
    return new Promise(resolve => {
      var pkcs8Key = stringToArrayBuffer(fromBase64(base64Key))
      crypto.subtle.importKey('pkcs8', pkcs8Key, ALGO, true, ['sign']).then((key) => {
        resolve(new ECDSA({ privateKey: key }))
      })
    })
  }

  /* ========================================================================== *
   * SIGNING / VALIDATION                                                       *
   * ========================================================================== */

  ECDSA.prototype.sign = function sign(message /*: string */) /*: Promise<string> */ {
    return new Promise(resolve => {
      crypto.subtle.sign(SIGN_ALGO, this.privateKey, stringToArrayBuffer(message)).then(signature => {
        resolve(toBase64(arrayBufferToString(signature)))
      })
    })
  }

  ECDSA.prototype.hashAndSign = function hashAndSign(message /*: string | Object */) /*: Promise<string> */ {
    return new Promise(resolve => {
      hash(message).then(hashed => {
        resolve(this.sign(hashed))
      })
    })
  }

  ECDSA.prototype.verify = function verify(message /*: string */, signature /*: string */) /*: Promise<Boolean> */ {
    var signatureBuffer = stringToArrayBuffer(fromBase64(signature))
    return crypto.subtle.verify(
      SIGN_ALGO,
      this.publicKey,
      signatureBuffer,
      stringToArrayBuffer(message)
    )
  }

  ECDSA.prototype.hashAndVerify = function hashAndVerify(message /*: string | Object */, signature /*: string */) /*: Promise<Boolean> */ {
    return new Promise(resolve => {
      hash(message).then(hashed => {
        resolve(this.verify(hashed , signature))
      })
    })
  }

  /* ========================================================================== *
   * CONVERSION                                                                 *
   * ========================================================================== */

  ECDSA.prototype.asPublic = function asPublic() {
    if (!this.privateKey) return this
    return new ECDSA({ publicKey: this.publicKey })
  }

  ECDSA.prototype.toJWK = function toJWK() /*: Promise<Object> */ {
    return new Promise(resolve => {
      var key = this.privateKey ? this.privateKey : this.publicKey
      crypto.subtle.exportKey('jwk', key).then(jwk => {
        resolve({ kty: jwk.kty, crv: jwk.crv, d: jwk.d, x: jwk.x, y: jwk.y })
      })
    })
  }

  ECDSA.prototype.toBase64PrivateKey = function toBase64PrivateKey() /*: Promise<string> */ {
    return new Promise(resolve => {
      crypto.subtle.exportKey('pkcs8', this.privateKey).then(key => {
        resolve(toBase64(arrayBufferToString(key)))
      })
    })
  }

  ECDSA.prototype.toCompressedPublicKey = function toCompressedPublicKey() /*: Promise<Uint8Array> */ {
    return new Promise(resolve => {
      crypto.subtle.exportKey('raw', this.publicKey).then(key => {
        var rawKey = new Uint8Array(key)
        var x = new Uint8Array(rawKey.slice(1, rawKey.length / 2 + 1))
        var y = new Uint8Array(rawKey.slice(rawKey.length / 2 + 1))
        var compressedKey = new Uint8Array(x.length + 1)
        compressedKey[0] = 2 + (y[y.length - 1] & 1)
        compressedKey.set(x, 1)
        resolve(compressedKey)
      })
    })

  }

  ECDSA.prototype.toBase64CompressedPublicKey = function toBase64CompressedPublicKey() /*: Promise<string> */ {
    return new Promise(resolve => {
      this.toCompressedPublicKey().then(compressedKey => {
        resolve(toBase64(arrayBufferToString(compressedKey)))
      })
    })
  }

  /* ========================================================================== *
   * EXPORTS                                                                    *
   * ========================================================================== */

  if (typeof module !== 'undefined' && typeof module.exports !== 'undefined') {
    module.exports = ECDSA
  } else {
    if (typeof define === 'function' && define.amd) {
      define([], function() {
        return ECDSA
      })
    } else {
      window.ECDSA = ECDSA
    }
  }
})()
