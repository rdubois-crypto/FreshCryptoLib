const crypto = require('crypto')

const asn = require('asn1.js')
const BN = require('bn.js')

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

const curve = 'prime256v1' /* OpenSSL curve name */
const jwkCurve = 'P-256' /* JWK curve name */
const curveLength = Math.ceil(256 / 8) /* Byte length for validation */

/* ========================================================================== *
 * CLASS DEFINITION                                                           *
 * ========================================================================== */

function ECDSA({ x, y, d }) {
  if (x.length !== curveLength) throw new TypeError('Public EC Key point X of wrong length')
  if (y.length !== curveLength) throw new TypeError('Public EC Key point Y of wrong length')

  this.x = x
  this.y = y
  this.publicCodePoint = Buffer.concat([Buffer.from([0x04]), x, y])
  if (d) {
    this.d = d
    this.isPrivate = true
  } else {
    this.isPrivate = false
  }
}

/* ========================================================================== *
 * FACTORIES                                                                  *
 * ========================================================================== */

ECDSA.generateKey = function generateKeys() {
  const ecdh = crypto.createECDH(curve)
  ecdh.generateKeys()

  return new ECDSA({
    d: ecdh.getPrivateKey(),
    x: ecdh.getPublicKey().slice(1, curveLength + 1),
    y: ecdh.getPublicKey().slice(curveLength + 1)
  })
}

ECDSA.fromJWK = function fromJWK(jwk) {
  return new ECDSA({
    x: Buffer.from(jwk.x, 'base64'),
    y: Buffer.from(jwk.y, 'base64'),
    d: jwk.d ? Buffer.from(jwk.d, 'base64') : false
  })
}

ECDSA.fromCompressedPublicKey = function fromCompressedPublicKey(compressedKey, format = 'base64') {
  const key = crypto.ECDH.convertKey(compressedKey, curve, format, 'base64', 'uncompressed')
  const keyBuffer = Buffer.from(key, 'base64')
  return new ECDSA({
    x: keyBuffer.slice(1, curveLength + 1),
    y: keyBuffer.slice(curveLength + 1)
  })
}

/* ========================================================================== *
 * ASN.1                                                                      *
 * ========================================================================== */

const ASN1 = {
  Rfc5915: asn.define('Rfc5915Key', function() {
    this.seq().obj(
      this.key('version').int(),
      this.key('privateKey').octstr(),
      this.key('parameters')
        .optional()
        .explicit(0)
        .objid({
          '1 2 840 10045 3 1 7': 'prime256v1',
          '1 3 132 0 10': 'secp256k1',
          '1 3 132 0 34': 'secp384r1',
          '1 3 132 0 35': 'secp521r1'
        }),
      this.key('publicKey')
        .optional()
        .explicit(1)
        .bitstr()
    )
  }),
  Pkcs8: asn.define('Pkcs8Key', function() {
    this.seq().obj(
      this.key('version').int(),
      this.key('algorithmIdentifier')
        .seq()
        .obj(
          this.key('privateKeyType').objid({
            '1 2 840 10045 2 1': 'EC'
          }),
          this.key('parameters').objid({
            '1 2 840 10045 3 1 7': 'prime256v1',
            '1 3 132 0 10': 'secp256k1',
            '1 3 132 0 34': 'secp384r1',
            '1 3 132 0 35': 'secp521r1'
          })
        ),
      this.key('privateKey').octstr()
    )
  }),
  Spki: asn.define('SpkiKey', function() {
    this.seq().obj(
      this.key('algorithmIdentifier')
        .seq()
        .obj(
          this.key('publicKeyType').objid({
            '1 2 840 10045 2 1': 'EC'
          }),
          this.key('parameters').objid({
            '1 2 840 10045 3 1 7': 'prime256v1',
            '1 3 132 0 10': 'secp256k1',
            '1 3 132 0 34': 'secp384r1',
            '1 3 132 0 35': 'secp521r1'
          })
        ),
      this.key('publicKey').bitstr()
    )
  }),
  EcdsaDerSig: asn.define('ECPrivateKey', function() {
    return this.seq().obj(this.key('r').int(), this.key('s').int())
  })
}

/* ========================================================================== *
 * SIGNING / VALIDATION                                                       *
 * ========================================================================== */

function hash(object) {
  return crypto
    .createHash('sha256')
    .update(typeof object === 'string' ? object : JSON.stringify(object))
    .digest('base64')
}

ECDSA.prototype.sign = function sign(message, format = 'base64') {
  function removeDerEncoding(signatureBuffer) {
    const { r, s } = ASN1.EcdsaDerSig.decode(signatureBuffer, 'der')
    return Buffer.concat([r.toBuffer(), s.toBuffer()]).toString(format)
  }

  if (!this.isPrivate) throw new Error('EC Private Key needed to sign')
  const sign = crypto.createSign('RSA-SHA256') // RSA works with EC keys, too
  sign.write(message)
  sign.end()
  return removeDerEncoding(sign.sign(this.toPEM()))
}

ECDSA.prototype.hashAndSign = async function hashAndSign(message, format = 'base64') {
  return this.sign(await hash(message), format)
}

ECDSA.prototype.verify = function verify(message, signature, format = 'base64') {
  function signatureToDer(signatureBuffer) {
    const r = new BN(signatureBuffer.slice(0, curveLength).toString('hex'), 16, 'be')
    const s = new BN(signatureBuffer.slice(curveLength).toString('hex'), 16, 'be')
    return ASN1.EcdsaDerSig.encode({ r, s }, 'der')
  }

  const verify = crypto.createVerify('RSA-SHA256') // RSA works with EC keys, too
  verify.write(message)
  verify.end()
  const key = this.isPrivate ? this.asPublic() : this
  const signatureBuffer = Buffer.from(signature, format)
  return verify.verify(
    key.toPEM(),
    signatureBuffer.length <= 2 * curveLength ? signatureToDer(signatureBuffer) : signatureBuffer,
    format
  )
}

ECDSA.prototype.hashAndVerify = async function hashAndVerify(
  message,
  signature,
  format = 'base64'
) {
  return this.verify(await hash(message), signature, format)
}

/* ========================================================================== *
 * CONVERSION                                                                 *
 * ========================================================================== */

ECDSA.prototype.asPublic = function asPublic() {
  if (!this.isPrivate) return this
  return new ECDSA({ x: this.x, y: this.y })
}

ECDSA.prototype.toCompressedPublicKey = function toCompressedPublicKey(format = 'base64') {
  return crypto.ECDH.convertKey(this.publicCodePoint, curve, 'base64', format, 'compressed')
}

ECDSA.prototype.toBuffer = function toBuffer(format) {
  if (format === 'pem') {
    return Buffer.from(this.toPEM(), 'ascii')
  } else if (this.isPrivate) {
    // Strip leading zeroes from private key
    let d = this.d
    while (d[0] === 0) d = d.slice(1)

    // Encode in PKCS8
    return ASN1.Pkcs8.encode(
      {
        version: 0,
        algorithmIdentifier: {
          privateKeyType: 'EC',
          parameters: curve
        },
        // Private key is RFC5915 minus curve
        privateKey: ASN1.Rfc5915.encode(
          {
            version: 1,
            privateKey: d,
            publicKey: { data: this.publicCodePoint }
          },
          'der'
        )
      },
      'der'
    )
  } else {
    return ASN1.Spki.encode(
      {
        algorithmIdentifier: {
          publicKeyType: 'EC',
          parameters: curve
        },
        publicKey: { data: this.publicCodePoint }
      },
      'der'
    )
  }
}

ECDSA.prototype.toPEM = function toPEM() {
  return this.isPrivate
    ? '-----BEGIN PRIVATE KEY-----\n' +
        this.toBuffer('pkcs8')
          .toString('base64')
          .match(/.{1,64}/g)
          .join('\n') +
        '\n-----END PRIVATE KEY-----\n'
    : '-----BEGIN PUBLIC KEY-----\n' +
        this.toBuffer('spki')
          .toString('base64')
          .match(/.{1,64}/g)
          .join('\n') +
        '\n-----END PUBLIC KEY-----\n'
}

ECDSA.prototype.toJWK = function toJWK() {
  function urlsafe(buffer) {
    return buffer
      .toString('base64')
      .replace(/\+/g, '-')
      .replace(/\//g, '_')
      .replace(/=/g, '')
  }

  const jwk = {
    kty: 'EC',
    crv: jwkCurve,
    x: urlsafe(this.x),
    y: urlsafe(this.y)
  }

  if (this.isPrivate) {
    let d = this.d
    if (d.length < curveLength) {
      const remaining = curveLength - d.length
      d = Buffer.concat([Buffer.alloc(remaining).fill(0), d])
    }
    jwk.d = urlsafe(d)
  }

  return jwk
}

/* ========================================================================== *
 * EXPORTS                                                                    *
 * ========================================================================== */

module.exports = ECDSA
