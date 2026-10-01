/**
 * @file sha1.h
 * @brief Self-contained public-domain SHA-1 (FIPS 180-1).
 *
 * Used to emit the PSI mzML spec's <fileChecksum> element during raw->mzML
 * extraction. No external dependency (no OpenSSL). SHA-1 is cryptographically
 * broken; it is used only because the mzML schema mandates it for that
 * element.
 */

#ifndef SHA1_H
#define SHA1_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
   uint32_t state[5];
   uint64_t bitlen;
   uint8_t data[64];
} sha1_ctx;

void sha1_init(sha1_ctx* c);
void sha1_update(sha1_ctx* c, const void* data, size_t len);
void sha1_final(sha1_ctx* c, uint8_t out[20]); /* big-endian digest */
/* convenience — one-shot, returns lowercase hex in out (41 bytes incl. NUL) */
void sha1_hex(const void* data, size_t len, char out[41]);

#endif /* SHA1_H */
