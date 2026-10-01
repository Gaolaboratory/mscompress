/**
 * @file sha1.c
 * @brief Public-domain SHA-1 (FIPS 180-1), in the style of Steve Reid's
 *        reference implementation: an 80-round compression function over
 *        32-bit big-endian words with a left-rotate.
 *
 * The context stores only the chaining state, the total bit length and the
 * current partial 64-byte block; the partial-block fill level is derived as
 * (bitlen >> 3) & 63, so no separate counter is needed.
 */

#include "sha1.h"

#include <string.h>

#define rol(v, bits) (((v) << (bits)) | ((v) >> (32 - (bits))))

static uint32_t
be32_load(const uint8_t* p)
{
   return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) |
          ((uint32_t)p[2] << 8) | (uint32_t)p[3];
}

static void
be32_store(uint8_t* p, uint32_t v)
{
   p[0] = (uint8_t)(v >> 24);
   p[1] = (uint8_t)(v >> 16);
   p[2] = (uint8_t)(v >> 8);
   p[3] = (uint8_t)v;
}

/* Hash one 64-byte block. */
static void
sha1_transform(uint32_t state[5], const uint8_t block[64])
{
   uint32_t w[80];
   uint32_t a, b, c, d, e;

   for (int i = 0; i < 16; i++)
      w[i] = be32_load(block + i * 4);
   for (int i = 16; i < 80; i++)
      w[i] = rol(w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16], 1);

   a = state[0];
   b = state[1];
   c = state[2];
   d = state[3];
   e = state[4];

   for (int i = 0; i < 80; i++) {
      uint32_t f, k;
      if (i < 20) {
         f = (b & c) | ((~b) & d);
         k = 0x5A827999;
      } else if (i < 40) {
         f = b ^ c ^ d;
         k = 0x6ED9EBA1;
      } else if (i < 60) {
         f = (b & c) | (b & d) | (c & d);
         k = 0x8F1BBCDC;
      } else {
         f = b ^ c ^ d;
         k = 0xCA62C1D6;
      }
      uint32_t t = rol(a, 5) + f + e + k + w[i];
      e = d;
      d = c;
      c = rol(b, 30);
      b = a;
      a = t;
   }

   state[0] += a;
   state[1] += b;
   state[2] += c;
   state[3] += d;
   state[4] += e;
}

void
sha1_init(sha1_ctx* c)
{
   c->state[0] = 0x67452301;
   c->state[1] = 0xEFCDAB89;
   c->state[2] = 0x98BADCFE;
   c->state[3] = 0x10325476;
   c->state[4] = 0xC3D2E1F0;
   c->bitlen = 0;
}

void
sha1_update(sha1_ctx* c, const void* data, size_t len)
{
   const uint8_t* p = (const uint8_t*)data;
   size_t fill = (size_t)(c->bitlen >> 3) & 63;

   c->bitlen += (uint64_t)len << 3;

   /* Top up a partial block first, then whole blocks straight from input. */
   if (fill && fill + len >= 64) {
      memcpy(c->data + fill, p, 64 - fill);
      sha1_transform(c->state, c->data);
      p += 64 - fill;
      len -= 64 - fill;
      fill = 0;
   }
   while (len >= 64) {
      sha1_transform(c->state, p);
      p += 64;
      len -= 64;
   }
   if (len)
      memcpy(c->data + fill, p, len);
}

void
sha1_final(sha1_ctx* c, uint8_t out[20])
{
   size_t fill = (size_t)(c->bitlen >> 3) & 63;

   c->data[fill++] = 0x80;
   if (fill > 56) { /* not enough room for the 8-byte length: flush a block */
      memset(c->data + fill, 0, 64 - fill);
      sha1_transform(c->state, c->data);
      fill = 0;
   }
   memset(c->data + fill, 0, 56 - fill);
   for (int i = 0; i < 8; i++) /* 64-bit big-endian bit length */
      c->data[63 - i] = (uint8_t)(c->bitlen >> (8 * i));
   sha1_transform(c->state, c->data);

   for (int i = 0; i < 5; i++)
      be32_store(out + i * 4, c->state[i]);
}

void
sha1_hex(const void* data, size_t len, char out[41])
{
   static const char digits[] = "0123456789abcdef";
   sha1_ctx c;
   uint8_t digest[20];

   sha1_init(&c);
   sha1_update(&c, data, len);
   sha1_final(&c, digest);

   for (int i = 0; i < 20; i++) {
      out[i * 2] = digits[digest[i] >> 4];
      out[i * 2 + 1] = digits[digest[i] & 0x0F];
   }
   out[40] = '\0';
}
