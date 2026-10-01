/**
 * @file test_idx_writer.c
 * @brief Unit tests for the index-tracking mzML output writer (idxout_*) in
 *        src/extract.c.
 *
 * Runs standalone (no msz file needed): a small two-spectrum one-line
 * document in the raw-path layout is driven through
 * idxout_begin/mark/write/finish into a temp file, then read back and
 * verified against real bytes: index offsets must seek onto "<spectrum" in
 * the produced file, the offset count must equal the spectra emitted, and
 * <fileChecksum> must equal SHA-1 recomputed over <indexedmzML .. </indexList>
 * (that recomputation is also the safety net for the incremental
 * sha1_init/update/final path used by the writer).
 */

#include "../sha1.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "mscompress.h"

static int failures = 0;
static int checks = 0;

#define CHECK(cond, msg)                                             \
   do {                                                              \
      checks++;                                                      \
      if (!(cond)) {                                                 \
         fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, msg); \
         failures++;                                                 \
      }                                                              \
   } while (0)

#define RUN_TEST(fn)                                                   \
   do {                                                                \
      int before = failures;                                           \
      fn();                                                            \
      printf("%-58s %s\n", #fn, failures == before ? "PASS" : "FAIL");  \
   } while (0)

/* Max spectra driven through the writer in these tests. */
#define MAX_SPEC 4
/* Longest id used (escaped attribute bytes are re-emitted verbatim). */
#define MAX_ID 128

/* --- shared verification -------------------------------------------------- */

/* Parses the <offset idRef="ID">N</offset> entries of the produced file.
   Fills offs[] and ids[] (ids[i] is malloc'd) and returns the entry count. */
static size_t parse_index_entries(const char* buf, uint64_t* offs,
                                  char* ids[], size_t max) {
   size_t n = 0;
   const char* p = strstr(buf, "<indexList");
   if (!p)
      return 0;
   while (n < max && (p = strstr(p, "<offset idRef=\"")) != NULL) {
      p += strlen("<offset idRef=\"");
      const char* vend = strchr(p, '"');
      if (!vend)
         return n;
      size_t id_len = (size_t)(vend - p);
      if (id_len >= MAX_ID)
         return n;
      ids[n] = malloc(id_len + 1);
      memcpy(ids[n], p, id_len);
      ids[n][id_len] = '\0';
      offs[n] = strtoull(vend + 2, NULL, 10); /* skips ">N */
      n++;
      p = vend;
   }
   return n;
}

/* Extracts the 40-hex-char value of <fileChecksum>..</fileChecksum>. */
static int get_file_checksum(const char* buf, char hex[41]) {
   const char* s = strstr(buf, "<fileChecksum>");
   if (!s)
      return -1;
   s += strlen("<fileChecksum>");
   const char* e = strstr(s, "</fileChecksum>");
   if (!e || (size_t)(e - s) != 40)
      return -1;
   memcpy(hex, s, 40);
   hex[40] = '\0';
   return 0;
}

/* Verifies a wrapped document the way a consumer (pyteomics) would:
   - opens with our declaration, then <indexedmzML
   - ends with "</indexedmzML>\n"
   - indexList carries exactly n_exp offsets, strictly increasing, whose
     idRefs match the expected (escaped) ids in order
   - seeking each offset in the file (fd) lands on "<spectrum", and the
     spectrum found there carries the matching id="..."
   - <fileChecksum> equals SHA-1 over [first byte of <indexedmzML .. end of
     </indexList>] */
static void verify_wrapped_output(int fd, const char* buf, size_t len,
                                  const char* const* exp_ids, size_t n_exp) {
   CHECK(len > 40, "output is non-trivial");
   if (len <= 40)
      return;

   CHECK(strncmp(buf, "<?xml ", 6) == 0, "output starts with <?xml");
   const char* wrap_start = strstr(buf, "<indexedmzML");
   CHECK(wrap_start != NULL, "output contains <indexedmzML");
   CHECK(strncmp(buf, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n", 39) == 0,
         "our XML declaration is first");
   CHECK(wrap_start == buf + 39, "<indexedmzML follows the declaration");
   CHECK(strcmp(buf + len - strlen("</indexedmzML>\n"), "</indexedmzML>\n") == 0,
         "output ends with </indexedmzML>\\n");

   CHECK(strstr(buf, "<indexList count=\"1\">") != NULL,
         "indexList declares one index");
   CHECK(strstr(buf, "<index name=\"spectrum\">") != NULL,
         "spectrum index is named");

   uint64_t offs[MAX_SPEC];
   char* ids[MAX_SPEC];
   size_t n = parse_index_entries(buf, offs, ids, MAX_SPEC);
   CHECK(n == n_exp, "index count equals emitted spectra");
   if (n != n_exp)
      return;

   uint64_t prev = 0;
   for (size_t i = 0; i < n; i++) {
      CHECK(offs[i] > prev, "offsets strictly increasing");
      CHECK(offs[i] < len, "offset inside file");
      prev = offs[i];

      CHECK(strcmp(ids[i], exp_ids[i]) == 0, "idRef matches spectrum id");

      /* Real seek: this is what an index consumer does. */
      char window[512];
      if (lseek(fd, (off_t)offs[i], SEEK_SET) == (off_t)-1) {
         CHECK(0, "seek to spectrum offset");
         continue;
      }
      ssize_t got = read(fd, window, sizeof(window) - 1);
      CHECK(got >= 9, "read at spectrum offset");
      if (got < 9)
         continue;
      window[got] = '\0';
      CHECK(strncmp(window, "<spectrum", 9) == 0,
            "offset lands on '<' of <spectrum");

      /* The spectrum reached through the offset is the one the idRef
         promises: its id=" attribute must carry the same escaped bytes. */
      char want[MAX_ID + 8];
      snprintf(want, sizeof(want), "id=\"%s\"", exp_ids[i]);
      CHECK(strstr(window, want) != NULL,
            "seeking the offset finds the idRef'd spectrum");
   }
   for (size_t i = 0; i < n; i++)
      free(ids[i]);

   /* Checksum: PSI convention — SHA-1 from the first byte of <indexedmzML
      through the last byte of </indexList>, recomputed here in one shot. */
   const char* il_end = strstr(buf, "</indexList>");
   CHECK(il_end != NULL, "found </indexList>");
   if (!il_end)
      return;
   size_t span = (size_t)(il_end + strlen("</indexList>") - wrap_start);
   char hex[41];
   sha1_hex(wrap_start, span, hex);
   char emitted[41];
   CHECK(get_file_checksum(buf, emitted) == 0, "fileChecksum present");
   CHECK(strcmp(hex, emitted) == 0,
         "fileChecksum equals SHA-1 over <indexedmzML .. </indexList>");

   /* <indexListOffset> is required by the indexedmzML schema and must name
      the absolute position where our <indexList begins. */
   const char* ilo = strstr(buf, "<indexListOffset>");
   CHECK(ilo != NULL, "indexListOffset present");
   const char* il_start = strstr(buf, "<indexList count=\"1\">");
   CHECK(il_start != NULL, "found our <indexList");
   if (ilo && il_start) {
      uint64_t want = (uint64_t)(il_start - buf);
      CHECK(strtoull(ilo + strlen("<indexListOffset>"), NULL, 10) == want,
            "indexListOffset equals position of <indexList");
      CHECK(ilo > il_end, "indexListOffset emitted after </indexList>");
   }
}

/* --- driver ---------------------------------------------------------------- */

/* Pushes a document through the writer exactly the way extract_msz feeds
   it: the header chunk through idxout_begin, one mark+write per spectrum
   (marking after any leading whitespace so the offset points at '<'),
   the footer chunk, then idxout_finish. Reads the file back and returns it
   malloc'd (+1 byte, NUL-terminated for string parsing); *out_len is the
   real file size. fd stays open for verify_wrapped_output when the caller
   passes a non-NULL pfd. */
static char* run_idxout(const char* header, const char* const* spec_xml,
                        const char* const* spec_ids, size_t n_spec,
                        const char* footer, size_t* out_len, int* pfd) {
   *out_len = 0;
   if (pfd)
      *pfd = -1;

   char tmpl[] = "/tmp/idx_writer_test_XXXXXX";
   int fd = mkstemp(tmpl);
   CHECK(fd >= 0, "mkstemp for idxout sink");
   if (fd < 0)
      return NULL;
   unlink(tmpl); /* fd stays usable until close(); no temp files left behind */

   idxout_t io;
   CHECK(idxout_begin(&io, fd, header, strlen(header)) == 0,
         "idxout_begin succeeds");

   for (size_t i = 0; i < n_spec; i++) {
      const char* s = spec_xml[i];
      size_t lead = 0;
      while (s[lead] == ' ' || s[lead] == '\t' || s[lead] == '\r' ||
             s[lead] == '\n')
         lead++;
      if (lead)
         CHECK(idxout_write(&io, s, lead) == 0, "idxout_write leading ws");
      CHECK(idxout_mark_spectrum(&io, spec_ids[i]) == 0,
            "idxout_mark_spectrum");
      CHECK(idxout_write(&io, s + lead, strlen(s + lead)) == 0,
            "idxout_write spectrum");
   }
   CHECK(idxout_write(&io, footer, strlen(footer)) == 0,
         "idxout_write footer");
   CHECK(idxout_finish(&io) == 0, "idxout_finish");

   struct stat st;
   fstat(fd, &st);
   char* buf = malloc((size_t)st.st_size + 1);
   CHECK(buf != NULL, "malloc output readback");
   if (!buf) {
      close(fd);
      return NULL;
   }
   lseek(fd, 0, SEEK_SET);
   ssize_t got = read(fd, buf, (size_t)st.st_size);
   CHECK(got == st.st_size, "read back full output");
   buf[got > 0 ? got : 0] = '\0';
   *out_len = (size_t)got;

   if (pfd)
      *pfd = fd; /* caller verifies and closes */
   else
      close(fd);
   return buf;
}

/* --- tests ------------------------------------------------------------------ */

/* Raw-path layout: bare <mzML> one-liner (Task-1 newline-terminated stream)
   must come out <indexedmzML>-wrapped with a true spectrum index. */
static void test_wrap_two_spectra(void) {
   const char* header =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<spectrumList count=\"2\">";
   const char* specs[] = {
      "<spectrum index=\"0\" id=\"scan=1\" defaultArrayLength=\"4\">"
      "<binary>AQIDBA==</binary></spectrum>\n",
      "<spectrum index=\"1\" id=\"scan=2\" defaultArrayLength=\"2\">"
      "<binary>BQY=</binary></spectrum>\n",
   };
   const char* ids[] = {"scan=1", "scan=2"};
   const char* footer = "</spectrumList></mzML>\n";

   size_t len = 0;
   int fd = -1;
   char* buf = run_idxout(header, specs, ids, 2, footer, &len, &fd);
   if (!buf)
      return;

   /* The stored declaration must not survive: we emit our own. */
   const char* decl = buf;
   int decl_count = 0;
   while ((decl = strstr(decl, "<?xml")) != NULL) {
      decl_count++;
      decl++;
   }
   CHECK(decl_count == 1, "exactly one <?xml declaration in output");

   /* The index block starts on its own line after </mzML>. */
   CHECK(strstr(buf, "</mzML>\n<indexList") != NULL,
         "indexList follows </mzML> on a new line");

   verify_wrapped_output(fd, buf, len, ids, 2);
   close(fd);
   free(buf);
}

/* Filtered extraction: only spectra actually emitted are indexed. */
static void test_wrap_filtered_single_spectrum(void) {
   const char* header =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<spectrumList count=\"1\">";
   const char* specs[] = {
      "<spectrum index=\"1\" id=\"scan=7\">"
      "<binary>BQY=</binary></spectrum>\n",
   };
   const char* ids[] = {"scan=7"};
   const char* footer = "</spectrumList></mzML>\n";

   size_t len = 0;
   int fd = -1;
   char* buf = run_idxout(header, specs, ids, 1, footer, &len, &fd);
   if (!buf)
      return;
   verify_wrapped_output(fd, buf, len, ids, 1);
   close(fd);
   free(buf);
}

/* Escaped id attribute bytes are re-emitted verbatim inside idRef="...". */
static void test_wrap_escaped_id_roundtrips(void) {
   const char* header =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<spectrumList count=\"1\">";
   const char* specs[] = {
      "<spectrum index=\"0\" id=\"scan=1&amp;controllerType=0\">"
      "<binary>AQ==</binary></spectrum>\n",
   };
   const char* ids[] = {"scan=1&amp;controllerType=0"};
   const char* footer = "</spectrumList></mzML>\n";

   size_t len = 0;
   int fd = -1;
   char* buf = run_idxout(header, specs, ids, 1, footer, &len, &fd);
   if (!buf)
      return;
   CHECK(strstr(buf, "<offset idRef=\"scan=1&amp;controllerType=0\">") !=
         NULL, "escaped id re-emitted verbatim");
   verify_wrapped_output(fd, buf, len, ids, 1);
   close(fd);
   free(buf);
}

/* Source document without a trailing newline before the indexList: the
   writer must put the index on its own line and the checksum must still
   span exactly <indexedmzML .. </indexList>. */
static void test_wrap_without_footer_newline(void) {
   const char* header =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"><spectrumList count=\"1\">";
   const char* specs[] = {
      "<spectrum index=\"0\" id=\"scan=1\"><binary>AQ==</binary></spectrum>",
   };
   const char* ids[] = {"scan=1"};
   const char* footer = "</spectrumList></mzML>"; /* no newline */

   size_t len = 0;
   int fd = -1;
   char* buf = run_idxout(header, specs, ids, 1, footer, &len, &fd);
   if (!buf)
      return;
   CHECK(strstr(buf, "</mzML>\n<indexList") != NULL,
         "newline inserted before indexList");
   verify_wrapped_output(fd, buf, len, ids, 1);
   close(fd);
   free(buf);
}

/* Header with no <?xml declaration at all still wraps. */
static void test_wrap_declless_header(void) {
   const char* header = "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
                        "<spectrumList count=\"1\">";
   const char* specs[] = {
      "<spectrum index=\"0\" id=\"scan=1\"><binary>AQ==</binary></spectrum>\n",
   };
   const char* ids[] = {"scan=1"};
   const char* footer = "</spectrumList></mzML>\n";

   size_t len = 0;
   int fd = -1;
   char* buf = run_idxout(header, specs, ids, 1, footer, &len, &fd);
   if (!buf)
      return;
   verify_wrapped_output(fd, buf, len, ids, 1);
   close(fd);
   free(buf);
}

/* Thermo-style nativeIDs contain spaces (raw2ms emits them verbatim); the
   id must travel through idRef="…" byte-exact. Seen in the wild as
   id="controllerType=0 controllerNumber=1 scan=N". */
static void test_wrap_spaced_thermo_id(void) {
   const char* header =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<spectrumList count=\"1\">";
   const char* specs[] = {
      "<spectrum index=\"4\" id=\"controllerType=0 controllerNumber=1 scan=5\""
      " defaultArrayLength=\"3\"><binary>AQID</binary></spectrum>\n",
   };
   const char* ids[] = {"controllerType=0 controllerNumber=1 scan=5"};
   const char* footer = "</spectrumList></mzML>\n";

   size_t len = 0;
   int fd = -1;
   char* buf = run_idxout(header, specs, ids, 1, footer, &len, &fd);
   if (!buf)
      return;
   CHECK(strstr(buf,
                "<offset idRef=\"controllerType=0 controllerNumber=1 "
                "scan=5\">") != NULL,
         "spaced nativeID re-emitted verbatim in idRef");
   verify_wrapped_output(fd, buf, len, ids, 1);
   close(fd);
   free(buf);
}

/* Leading whitespace before <spectrum must not shift the recorded offset. */
static void test_wrap_leading_whitespace_chunk(void) {
   const char* header =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<spectrumList count=\"1\">";
   const char* specs[] = {
      "\n  <spectrum index=\"0\" id=\"scan=1\">"
      "<binary>AQ==</binary></spectrum>\n",
   };
   const char* ids[] = {"scan=1"};
   const char* footer = "</spectrumList></mzML>\n";

   size_t len = 0;
   int fd = -1;
   char* buf = run_idxout(header, specs, ids, 1, footer, &len, &fd);
   if (!buf)
      return;
   verify_wrapped_output(fd, buf, len, ids, 1);
   close(fd);
   free(buf);
}

/* mzML-path msz: the stored document is already <indexedmzML>-wrapped, so
   the writer must pass it through byte-identically and never double-wrap.
   finish() only guarantees the trailing newline. */
static void test_passthrough_already_indexed(void) {
   const char* doc =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<spectrumList count=\"1\">"
      "<spectrum index=\"0\" id=\"scan=1\"><binary>AQ==</binary></spectrum>"
      "</spectrumList></mzML>"
      "<indexList count=\"1\"><index name=\"spectrum\">"
      "<offset idRef=\"scan=1\">123</offset>"
      "</index></indexList>"
      "<fileChecksum>0123456789abcdef0123456789abcdef01234567</fileChecksum>"
      "</indexedmzML>"; /* no trailing newline */

   char tmpl[] = "/tmp/idx_writer_pass_XXXXXX";
   int fd = mkstemp(tmpl);
   CHECK(fd >= 0, "mkstemp for passthrough");
   if (fd < 0)
      return;

   idxout_t io;
   CHECK(idxout_begin(&io, fd, doc, strlen(doc)) == 0,
         "idxout_begin passthrough");
   CHECK(io.wrap == 0, "already-indexed source selects passthrough");
   CHECK(idxout_finish(&io) == 0, "idxout_finish passthrough");

   struct stat st;
   fstat(fd, &st);
   CHECK((size_t)st.st_size == strlen(doc) + 1,
         "passthrough output is input plus exactly one newline");
   char* buf = malloc((size_t)st.st_size + 1);
   CHECK(buf != NULL, "malloc passthrough readback");
   if (!buf) {
      close(fd);
      unlink(tmpl);
      return;
   }
   lseek(fd, 0, SEEK_SET);
   ssize_t got = read(fd, buf, (size_t)st.st_size);
   CHECK(got == st.st_size, "read back passthrough output");
   buf[st.st_size] = '\0';

   CHECK(memcmp(buf, doc, strlen(doc)) == 0,
         "passthrough emits the source verbatim");
   CHECK(buf[strlen(doc)] == '\n', "finish appends the missing newline");

   /* No second wrapper: exactly one <indexedmzML and one <?xml in output. */
   const char* p = buf;
   int wraps = 0, decls = 0;
   while ((p = strstr(p, "<indexedmzML")) != NULL) {
      wraps++;
      p++;
   }
   p = buf;
   while ((p = strstr(p, "<?xml")) != NULL) {
      decls++;
      p++;
   }
   CHECK(wraps == 1, "no double wrap in passthrough");
   CHECK(decls == 1, "no injected declaration in passthrough");
   CHECK(strstr(buf, "<indexListOffset") == NULL,
         "no indexListOffset injected in passthrough");

   free(buf);
   close(fd);
   unlink(tmpl);
}

/* Same, but the source already ends with '\n': nothing must be appended. */
static void test_passthrough_trailing_newline_kept(void) {
   const char* doc =
      "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
      "<indexedmzML xmlns=\"http://psi.hupo.org/ms/mzml\">"
      "<mzML><spectrumList count=\"0\"></spectrumList></mzML>"
      "</indexedmzML>\n";

   char tmpl[] = "/tmp/idx_writer_pass2_XXXXXX";
   int fd = mkstemp(tmpl);
   CHECK(fd >= 0, "mkstemp for passthrough newline case");
   if (fd < 0)
      return;

   idxout_t io;
   CHECK(idxout_begin(&io, fd, doc, strlen(doc)) == 0, "idxout_begin");
   CHECK(io.wrap == 0, "indexed source is passthrough");
   CHECK(idxout_finish(&io) == 0, "idxout_finish");

   struct stat st;
   fstat(fd, &st);
   CHECK((size_t)st.st_size == strlen(doc),
         "no extra byte appended when input ends in newline");
   lseek(fd, 0, SEEK_SET);
   char* buf = malloc(strlen(doc) + 1);
   CHECK(buf != NULL, "malloc");
   if (buf) {
      ssize_t got = read(fd, buf, strlen(doc));
      CHECK(got == (ssize_t)strlen(doc), "read back");
      CHECK(memcmp(buf, doc, strlen(doc)) == 0, "byte-identical passthrough");
      CHECK(strstr(buf, "<indexListOffset") == NULL,
            "no indexListOffset injected in passthrough");
      free(buf);
   }
   close(fd);
   unlink(tmpl);
}

int main(void) {
   RUN_TEST(test_wrap_two_spectra);
   RUN_TEST(test_wrap_filtered_single_spectrum);
   RUN_TEST(test_wrap_escaped_id_roundtrips);
   RUN_TEST(test_wrap_spaced_thermo_id);
   RUN_TEST(test_wrap_without_footer_newline);
   RUN_TEST(test_wrap_declless_header);
   RUN_TEST(test_wrap_leading_whitespace_chunk);
   RUN_TEST(test_passthrough_already_indexed);
   RUN_TEST(test_passthrough_trailing_newline_kept);

   printf("%d checks, %d failures\n", checks, failures);
   return failures ? 1 : 0;
}
