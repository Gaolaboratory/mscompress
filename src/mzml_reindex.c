/**
 * @file mzml_reindex.c
 * @brief Streaming repair of an indexedmzML document's index on decompress.
 *
 * mscompress stores binary arrays decoded and re-encodes them on the way out.
 * When the source arrays were zlib-compressed, the re-deflated streams rarely
 * have the source's exact length, so every byte after the first array moves:
 * the <offset> entries of <indexList>, <indexListOffset> and the
 * <fileChecksum> SHA-1 the source carried all go stale.
 *
 * The reindexer sits between decompress_msz and the output fd. The decoded
 * divisions arrive in document order; it writes them through, recording the
 * true output offset of each <spectrum>/<chromatogram> element by id, and
 * holds back everything after </mzML> (the index tail). At the end it patches
 * that tail in place — each <offset>'s number by its idRef, <indexListOffset>,
 * and <fileChecksum> (SHA-1 from byte 0 through the "<fileChecksum>" open tag,
 * per mzML1.1.x_idx.xsd) — and writes it. Patching rather than regenerating
 * keeps the tail's layout (whitespace, extra attributes) as the source had it.
 *
 * The repair runs only when the caller reports that some array came out a
 * different length than the source's encodedLength — i.e. the output bytes
 * provably moved. Otherwise the tail is written verbatim, so a document whose
 * arrays round-trip unchanged stays byte-identical to its source even when
 * that source's own index or checksum was wrong.
 *
 * A document that does not open with <indexedmzML is passed through verbatim
 * with no scanning or hashing.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mscompress.h"

#ifdef _WIN32
#include <io.h>
#define write _write
#else
#include <unistd.h>
#endif

enum { RX_UNKNOWN = 0, RX_PASSTHROUGH, RX_INDEXED, RX_TAIL };
enum { RX_SPECTRUM = 0, RX_CHROMATOGRAM = 1 };

typedef struct {
   char* id;      /* raw (still XML-escaped) id attribute value */
   uint64_t off;  /* byte offset of the element's '<' in the output */
   int type;      /* RX_SPECTRUM / RX_CHROMATOGRAM */
} rx_entry_t;

struct mzml_reindex {
   int fd;
   int mode;
   uint64_t written; /* bytes emitted to fd so far */
   sha1_ctx sha;

   /* Bytes from an incomplete trailing tag ("<..." with no '>' yet), held
    * until the next chunk completes it. */
   char* pend;
   size_t pend_len, pend_cap;

   /* Everything after </mzML>, held for patching. */
   char* tail;
   size_t tail_len, tail_cap;

   rx_entry_t* ents;
   size_t n_ents, cap_ents;

   /* open-addressing hash: slot -> ents index + 1 (0 = empty) */
   size_t* slots;
   size_t n_slots;

   long n_patched, n_missing; /* stats for the report */
};

static int rx_append(char** buf, size_t* len, size_t* cap, const char* src,
                     size_t n) {
   if (*len + n > *cap) {
      size_t nc = *cap ? *cap : 4096;
      while (nc < *len + n) nc *= 2;
      char* nb = realloc(*buf, nc);
      if (!nb) return -1;
      *buf = nb;
      *cap = nc;
   }
   memcpy(*buf + *len, src, n);
   *len += n;
   return 0;
}

static int rx_emit(mzml_reindex_t* rx, const char* p, size_t n) {
   const char* q = p;
   size_t left = n;
   while (left > 0) {
#ifdef _WIN32
      int w = write(rx->fd, q, (unsigned int)left);
#else
      ssize_t w = write(rx->fd, q, left);
#endif
      if (w <= 0) {
         error("mzml_reindex: write of %zu bytes failed.\n", n);
         return -1;
      }
      q += w;
      left -= (size_t)w;
   }
   if (rx->mode != RX_PASSTHROUGH) sha1_update(&rx->sha, p, n);
   rx->written += n;
   return 0;
}

static uint64_t rx_hash(const char* s, size_t n, int type) {
   uint64_t h = 1469598103934665603ULL ^ (uint64_t)type;
   for (size_t i = 0; i < n; i++) {
      h ^= (unsigned char)s[i];
      h *= 1099511628211ULL;
   }
   return h;
}

static int rx_rehash(mzml_reindex_t* rx, size_t n_slots) {
   size_t* s = calloc(n_slots, sizeof(size_t));
   if (!s) return -1;
   for (size_t i = 0; i < rx->n_ents; i++) {
      rx_entry_t* e = &rx->ents[i];
      size_t k = rx_hash(e->id, strlen(e->id), e->type) & (n_slots - 1);
      while (s[k]) k = (k + 1) & (n_slots - 1);
      s[k] = i + 1;
   }
   free(rx->slots);
   rx->slots = s;
   rx->n_slots = n_slots;
   return 0;
}

static int rx_record(mzml_reindex_t* rx, const char* id, size_t id_len,
                     int type, uint64_t off) {
   if (rx->n_ents == rx->cap_ents) {
      size_t nc = rx->cap_ents ? rx->cap_ents * 2 : 1024;
      rx_entry_t* ne = realloc(rx->ents, nc * sizeof(rx_entry_t));
      if (!ne) return -1;
      rx->ents = ne;
      rx->cap_ents = nc;
   }
   char* dup = malloc(id_len + 1);
   if (!dup) return -1;
   memcpy(dup, id, id_len);
   dup[id_len] = '\0';
   rx->ents[rx->n_ents].id = dup;
   rx->ents[rx->n_ents].off = off;
   rx->ents[rx->n_ents].type = type;
   rx->n_ents++;
   /* keep load factor <= 1/2; first id wins on duplicates (lookup order) */
   if (rx->n_ents * 2 > rx->n_slots)
      return rx_rehash(rx, rx->n_slots ? rx->n_slots * 2 : 2048);
   size_t k = rx_hash(dup, id_len, type) & (rx->n_slots - 1);
   while (rx->slots[k]) k = (k + 1) & (rx->n_slots - 1);
   rx->slots[k] = rx->n_ents;
   return 0;
}

static int rx_lookup(const mzml_reindex_t* rx, const char* id, size_t id_len,
                     int type, uint64_t* off) {
   if (!rx->n_slots) return 0;
   size_t k = rx_hash(id, id_len, type) & (rx->n_slots - 1);
   while (rx->slots[k]) {
      const rx_entry_t* e = &rx->ents[rx->slots[k] - 1];
      if (e->type == type && strncmp(e->id, id, id_len) == 0 &&
          e->id[id_len] == '\0') {
         *off = e->off;
         return 1;
      }
      k = (k + 1) & (rx->n_slots - 1);
   }
   return 0;
}

static int rx_is_ws(char c) {
   return c == ' ' || c == '\t' || c == '\n' || c == '\r';
}

/* Finds attribute `name` (e.g. "id") inside the tag [p, end) — the name must
 * be preceded by whitespace so "id" does not match "nativeID"/"idRef".
 * Returns the value start and sets *vlen, or NULL. */
static const char* rx_attr(const char* p, const char* end, const char* name,
                           size_t* vlen) {
   size_t nl = strlen(name);
   for (const char* q = p; q + nl + 2 < end; q++) {
      if (!rx_is_ws(*q) || memcmp(q + 1, name, nl) != 0) continue;
      const char* r = q + 1 + nl;
      while (r < end && rx_is_ws(*r)) r++;
      if (r >= end || *r != '=') continue;
      r++;
      while (r < end && rx_is_ws(*r)) r++;
      if (r >= end || (*r != '"' && *r != '\'')) continue;
      char quote = *r++;
      const char* v = r;
      while (r < end && *r != quote) r++;
      if (r >= end) return NULL;
      *vlen = (size_t)(r - v);
      return v;
   }
   return NULL;
}

/* Does the tag at p (pointing at '<') start with `name` followed by a tag
 * delimiter? */
static int rx_tag_is(const char* p, const char* end, const char* name) {
   size_t nl = strlen(name);
   if ((size_t)(end - p) < nl + 2) return 0;
   if (memcmp(p + 1, name, nl) != 0) return 0;
   char c = p[1 + nl];
   return rx_is_ws(c) || c == '>' || c == '/';
}

/* Scans and emits a run of complete tags [p, p+n) that starts at output
 * offset rx->written. Switches to RX_TAIL after </mzML>; returns the number of
 * bytes consumed (the rest, if any, belongs to the tail), or -1 on error. */
static long rx_scan_emit(mzml_reindex_t* rx, const char* p, size_t n) {
   const char* end = p + n;
   const char* q = p;
   while (q < end) {
      const char* lt = memchr(q, '<', (size_t)(end - q));
      if (!lt) break;
      const char* gt = memchr(lt, '>', (size_t)(end - lt));
      if (!gt) break; /* caller guarantees completeness; be safe */
      int type = -1;
      if (rx_tag_is(lt, end, "spectrum"))
         type = RX_SPECTRUM;
      else if (rx_tag_is(lt, end, "chromatogram"))
         type = RX_CHROMATOGRAM;
      if (type >= 0) {
         size_t vlen;
         const char* v = rx_attr(lt, gt, "id", &vlen);
         if (v && rx_record(rx, v, vlen, type,
                            rx->written + (uint64_t)(lt - p)) != 0)
            return -1;
      } else if (gt - lt == 6 && memcmp(lt, "</mzML", 6) == 0) {
         size_t upto = (size_t)(gt + 1 - p);
         if (rx_emit(rx, p, upto) != 0) return -1;
         rx->mode = RX_TAIL;
         return (long)upto;
      }
      q = gt + 1;
   }
   if (rx_emit(rx, p, n) != 0) return -1;
   return (long)n;
}

mzml_reindex_t* mzml_reindex_new(int fd) {
   mzml_reindex_t* rx = calloc(1, sizeof(mzml_reindex_t));
   if (!rx) return NULL;
   rx->fd = fd;
   rx->mode = RX_UNKNOWN;
   sha1_init(&rx->sha);
   return rx;
}

/* Decides on the first chunk: an <indexedmzML root before <mzML means the
 * document carries an index worth repairing. */
static void rx_detect(mzml_reindex_t* rx, const char* buf, size_t len) {
   size_t lim = len < 65536 ? len : 65536;
   const char* q = buf;
   const char* end = buf + lim;
   while (q < end) {
      const char* lt = memchr(q, '<', (size_t)(end - q));
      if (!lt) break;
      if (rx_tag_is(lt, end, "indexedmzML")) {
         rx->mode = RX_INDEXED;
         return;
      }
      if (rx_tag_is(lt, end, "mzML")) break;
      q = lt + 1;
   }
   rx->mode = RX_PASSTHROUGH;
}

int mzml_reindex_write(mzml_reindex_t* rx, const char* buf, size_t len) {
   if (len == 0) return 0;
   if (rx->mode == RX_UNKNOWN) rx_detect(rx, buf, len);
   if (rx->mode == RX_PASSTHROUGH) return rx_emit(rx, buf, len);
   if (rx->mode == RX_TAIL)
      return rx_append(&rx->tail, &rx->tail_len, &rx->tail_cap, buf, len);

   /* RX_INDEXED. First complete a tag left open by the previous chunk. */
   if (rx->pend_len) {
      const char* gt = memchr(buf, '>', len);
      if (!gt)
         return rx_append(&rx->pend, &rx->pend_len, &rx->pend_cap, buf, len);
      size_t k = (size_t)(gt + 1 - buf);
      if (rx_append(&rx->pend, &rx->pend_len, &rx->pend_cap, buf, k) != 0)
         return -1;
      size_t plen = rx->pend_len;
      rx->pend_len = 0;
      long used = rx_scan_emit(rx, rx->pend, plen);
      if (used < 0) return -1;
      if ((size_t)used < plen &&
          rx_append(&rx->tail, &rx->tail_len, &rx->tail_cap,
                    rx->pend + used, plen - (size_t)used) != 0)
         return -1;
      buf += k;
      len -= k;
      if (rx->mode == RX_TAIL)
         return rx_append(&rx->tail, &rx->tail_len, &rx->tail_cap, buf, len);
   }

   /* Hold back a trailing "<..." whose '>' lies in the next chunk. */
   size_t safe = len;
   for (size_t i = len; i > 0; i--) {
      char c = buf[i - 1];
      if (c == '>') break;
      if (c == '<') {
         safe = i - 1;
         break;
      }
   }

   long used = rx_scan_emit(rx, buf, safe);
   if (used < 0) return -1;
   if (rx->mode == RX_TAIL)
      return rx_append(&rx->tail, &rx->tail_len, &rx->tail_cap, buf + used,
                       len - (size_t)used);
   if (safe < len)
      return rx_append(&rx->pend, &rx->pend_len, &rx->pend_cap, buf + safe,
                       len - safe);
   return 0;
}

/* Patches the held tail into a new buffer: <offset> numbers by idRef,
 * <indexListOffset>, and the <fileChecksum> digest. `base` is the output
 * offset of the tail's first byte. */
static int rx_patch_tail(mzml_reindex_t* rx, char** out, size_t* out_len) {
   const char* t = rx->tail;
   const char* end = t + rx->tail_len;
   char* o = NULL;
   size_t ol = 0, oc = 0;
   uint64_t base = rx->written;
   int64_t index_list_off = -1;
   int cur_type = RX_SPECTRUM;
   char num[48]; /* offsets, or a 40-char SHA-1 hex + NUL */

   const char* q = t;
   while (q < end) {
      const char* lt = memchr(q, '<', (size_t)(end - q));
      if (!lt) break;
      const char* gt = memchr(lt, '>', (size_t)(end - lt));
      if (!gt) break;

      if (rx_tag_is(lt, end, "indexList")) {
         if (rx_append(&o, &ol, &oc, q, (size_t)(lt - q)) != 0) goto fail;
         index_list_off = (int64_t)(base + ol);
         if (rx_append(&o, &ol, &oc, lt, (size_t)(gt + 1 - lt)) != 0)
            goto fail;
         q = gt + 1;
         continue;
      }
      if (rx_tag_is(lt, end, "index")) {
         size_t vlen;
         const char* v = rx_attr(lt, gt, "name", &vlen);
         if (v && vlen == 12 && memcmp(v, "chromatogram", 12) == 0)
            cur_type = RX_CHROMATOGRAM;
         else
            cur_type = RX_SPECTRUM;
      } else if (rx_tag_is(lt, end, "offset") ||
                 rx_tag_is(lt, end, "indexListOffset") ||
                 rx_tag_is(lt, end, "fileChecksum")) {
         /* value runs from gt+1 to the next '<' */
         const char* vs = gt + 1;
         const char* ve = memchr(vs, '<', (size_t)(end - vs));
         if (!ve) break;
         int n = -1;
         if (lt[1] == 'o') {
            size_t ilen;
            const char* id = rx_attr(lt, gt, "idRef", &ilen);
            uint64_t off;
            if (id && rx_lookup(rx, id, ilen, cur_type, &off)) {
               n = snprintf(num, sizeof(num), "%llu", (unsigned long long)off);
               rx->n_patched++;
            } else {
               rx->n_missing++;
            }
         } else if (lt[1] == 'i') {
            if (index_list_off >= 0)
               n = snprintf(num, sizeof(num), "%lld",
                            (long long)index_list_off);
         } else {
            /* SHA-1 over output byte 0 .. end of "<fileChecksum>" */
            if (rx_append(&o, &ol, &oc, q, (size_t)(vs - q)) != 0) goto fail;
            sha1_ctx c = rx->sha;
            sha1_update(&c, o, ol);
            uint8_t dg[20];
            sha1_final(&c, dg);
            for (int i = 0; i < 20; i++)
               snprintf(num + 2 * i, 3, "%02x", dg[i]);
            if (rx_append(&o, &ol, &oc, num, 40) != 0) goto fail;
            q = ve;
            continue;
         }
         if (n > 0) {
            if (rx_append(&o, &ol, &oc, q, (size_t)(vs - q)) != 0 ||
                rx_append(&o, &ol, &oc, num, (size_t)n) != 0)
               goto fail;
            q = ve;
            continue;
         }
      }
      if (rx_append(&o, &ol, &oc, q, (size_t)(gt + 1 - q)) != 0) goto fail;
      q = gt + 1;
   }
   if (q < end && rx_append(&o, &ol, &oc, q, (size_t)(end - q)) != 0)
      goto fail;
   *out = o;
   *out_len = ol;
   return 0;
fail:
   free(o);
   return -1;
}

int mzml_reindex_finish(mzml_reindex_t* rx, int repair) {
   int rv = 0;
   if (rx->pend_len) {
      /* unterminated tag at EOF: nothing to complete it, emit as is */
      if (rx->mode == RX_TAIL)
         rv |= rx_append(&rx->tail, &rx->tail_len, &rx->tail_cap, rx->pend,
                         rx->pend_len);
      else
         rv |= rx_emit(rx, rx->pend, rx->pend_len);
      rx->pend_len = 0;
   }
   if (rx->mode == RX_TAIL && rx->tail_len && !repair) {
      rv |= rx_emit(rx, rx->tail, rx->tail_len);
   } else if (rx->mode == RX_TAIL && rx->tail_len) {
      char* patched = NULL;
      size_t plen = 0;
      if (rx_patch_tail(rx, &patched, &plen) != 0) {
         error("mzml_reindex: out of memory patching the index.\n");
         return -1;
      }
      rv |= rx_emit(rx, patched, plen);
      free(patched);
      if (rx->n_missing)
         warning("mzml_reindex: %ld index entr%s had no matching element; "
                 "left unchanged.\n",
                 rx->n_missing, rx->n_missing == 1 ? "y" : "ies");
      print("\tRe-indexed mzML: %ld offsets, indexListOffset and fileChecksum "
            "recomputed.\n",
            rx->n_patched);
   }
   return rv ? -1 : 0;
}

void mzml_reindex_free(mzml_reindex_t* rx) {
   if (!rx) return;
   for (size_t i = 0; i < rx->n_ents; i++) free(rx->ents[i].id);
   free(rx->ents);
   free(rx->slots);
   free(rx->pend);
   free(rx->tail);
   free(rx);
}

/* === encodedLength repair for an assembled chunk ========================= */

long mzml_fix_encoded_lengths(char** buf, size_t* len) {
   static const char tag[] = "<binaryDataArray";
   static const char attr[] = "encodedLength";
   const size_t tl = sizeof(tag) - 1, al = sizeof(attr) - 1;
   const char* s = *buf;
   const char* end = s + *len;
   const char* copied = s; /* source bytes before this are already in out */
   char* out = NULL;
   size_t ol = 0, oc = 0;
   long n_fixed = 0;
   char num[32];

   const char* q = s;
   while ((size_t)(end - q) > tl) {
      const char* t = memchr(q, '<', (size_t)(end - q));
      if (!t || (size_t)(end - t) <= tl) break;
      if (memcmp(t, tag, tl) != 0 || !(rx_is_ws(t[tl]) || t[tl] == '>')) {
         q = t + 1;
         continue;
      }
      const char* gt = memchr(t, '>', (size_t)(end - t));
      if (!gt) break;
      q = gt + 1;

      size_t vlen;
      const char* v = rx_attr(t, gt, attr, &vlen);
      if (!v) continue;

      /* this array's <binary> text */
      const char* b0 = NULL;
      for (const char* r = gt; r + 8 <= end; r++) {
         r = memchr(r, '<', (size_t)(end - r));
         if (!r || r + 8 > end) break;
         if (memcmp(r, "<binary>", 8) == 0) {
            b0 = r + 8;
            break;
         }
         if (end - r >= 17 && memcmp(r, "</binaryDataArray", 17) == 0)
            break; /* array without <binary>: leave it */
      }
      if (!b0) continue;
      const char* b1 = memchr(b0, '<', (size_t)(end - b0));
      if (!b1) break;

      int n = snprintf(num, sizeof(num), "%zu", (size_t)(b1 - b0));
      if ((size_t)n == vlen && memcmp(v, num, vlen) == 0) {
         q = b1;
         continue;
      }
      if (rx_append(&out, &ol, &oc, copied, (size_t)(v - copied)) != 0 ||
          rx_append(&out, &ol, &oc, num, (size_t)n) != 0) {
         free(out);
         return -1;
      }
      copied = v + vlen;
      n_fixed++;
      q = b1;
   }

   if (n_fixed == 0) return 0;
   if (rx_append(&out, &ol, &oc, copied, (size_t)(end - copied)) != 0) {
      free(out);
      return -1;
   }
   free(*buf);
   *buf = out;
   *len = ol;
   return n_fixed;
}
