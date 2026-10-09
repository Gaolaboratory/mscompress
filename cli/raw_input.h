/**
 * @file raw_input.h
 * @brief Direct vendor raw file support: raw -> mzML (in memory) -> msz.
 *
 * Reads Thermo (.raw), Bruker timsTOF (TDF .d, TSF .d, BAF .d), SCIEX
 * (.wiff/.wiff2/.t2d) and Waters (.raw directory) files through the
 * compiled raw2ms C library
 * (libraw2ms_capi) loaded at runtime with dlopen()/LoadLibraryA(), then
 * synthesizes an mzML document and feeds it into the existing compression
 * pipeline (preprocess_mzml + compress_mzml).
 *
 * The mzML is materialized either fully in memory (when an OOM check says it
 * fits) or streamed spectrum-by-spectrum into a memfd/temp file that is then
 * mmapped, so peak heap usage stays bounded on large runs.
 */

#ifndef RAW_INPUT_H
#define RAW_INPUT_H

#include <stddef.h>
#include <stdint.h>

#include "mscompress.h"

/* Vendor codes mirrored from the raw2ms C ABI (RAW2MS_VENDOR_*). */
#define RAW_VENDOR_UNKNOWN 0
#define RAW_VENDOR_THERMO 1
#define RAW_VENDOR_BRUKER_TIMS 2
#define RAW_VENDOR_SCIEX 3
#define RAW_VENDOR_WATERS 4
#define RAW_VENDOR_SCIEX_T2D 5
#define RAW_VENDOR_BRUKER_TSF 6  /* Bruker .d without mobility (TSF) */
#define RAW_VENDOR_BRUKER_BAF 7  /* Bruker QTOF .d (BAF): micrOTOF, impact, maXis */

/* Polarity / peak-mode codes mirrored from the raw2ms C ABI. */
#define RAW_POLARITY_POSITIVE 1
#define RAW_POLARITY_NEGATIVE (-1)
#define RAW_PEAKS_CENTROID 1
#define RAW_PEAKS_PROFILE 2

/**
 * @brief One decoded spectrum.
 *
 * Field order and types mirror Raw2msSpectrum in the raw2ms C ABI — keep in
 * step. Absent floating-point values are NaN; absent ids are -1.
 */
typedef struct {
   uint64_t index;
   uint8_t ms_level;           /* 0 when the vendor file does not say */
   int8_t polarity;            /* RAW_POLARITY_* */
   int8_t peak_mode;           /* RAW_PEAKS_* */
   double rt_seconds;
   double drift_time_ms;       /* NaN unless per-spectrum IM (Waters) */
   double injection_time_ms;   /* NaN when unavailable */
   double precursor_mz;        /* NaN when not MS/MS or no single ion */
   double isolation_lower;     /* NaN when there is no window */
   double isolation_upper;
   int64_t precursor_id;       /* -1 when absent */
   int32_t charge;             /* -1 when absent */
   uint64_t n_peaks;
   const double* mz;           /* n_peaks entries, ascending */
   const float* intensity;     /* n_peaks entries */
   const float* mobility;      /* n_peaks entries or NULL (not emitted) */
   const char* native_id;      /* vendor-native id, NUL-terminated */
} raw_spec_t;

/**
 * @brief A spectrum's own precursor stage.
 *
 * Mirrors Raw2msPrecursor in the raw2ms C ABI (raw2ms_spectrum_precursor)
 * field-for-field -- keep in step. It carries what raw_spec_t cannot: the
 * isolation *target* (the window centre the quadrupole was set to), which
 * differs from the selected ion on a large minority of Thermo DDA scans
 * (monoisotopic correction) and is what mzML's MS:1000827 means.
 */
typedef struct {
   double mz;                /* selected ion; NaN when no single ion */
   double isolation_target;  /* NaN when unrecorded, as are the edges */
   double isolation_lower;
   double isolation_upper;
   double collision_energy;  /* eV; NaN when unrecorded */
   int64_t parent_index;     /* -1 when absent */
   int64_t precursor_id;     /* -1 when absent */
   int32_t charge;           /* -1 when absent */
   int32_t activation;       /* RAW2MS_ACTIVATION_*; 0 unknown */
} raw_prec_t;

/**
 * @brief Spectra are pulled one at a time so neither the in-memory nor the
 *        chunked sink ever holds more than one spectrum's mzML text.
 */
typedef struct {
   void* ctx;
   const raw_spec_t* (*fetch)(void* ctx, uint64_t index); /* NULL on error */
   void (*release)(void* ctx, const raw_spec_t* spec);    /* NULL if static */
   /* Detail for the most recent failed fetch (provider-specific), or NULL.
    * May be NULL itself when a provider has nothing to add. */
   const char* (*last_error)(void* ctx);
   /* Optional: fill `out` with `spec`'s own precursor stage and return 0, or
    * nonzero when the spectrum has none or the provider cannot say. NULL
    * when the provider has no such notion; the writer then centres the
    * window on the selected ion. */
   int (*precursor)(void* ctx, const raw_spec_t* spec, raw_prec_t* out);
} raw_provider_t;

/**
 * @brief A vendor run plus the metadata the mzML header needs.
 */
typedef struct {
   raw_provider_t provider;
   uint64_t n_spectra;
   const char* source_path;      /* vendor file the run was read from */
   const char* run_id;           /* may be NULL */
   const char* instrument_model; /* may be NULL or "" */
   const char* raw2ms_version;   /* version of the reader, may be NULL */
   int32_t vendor;               /* RAW_VENDOR_* */
   /* Set when the vendor file was corrupt and raw2ms recovered what it could
    * (--salvage): its account of what was lost. NULL for a healthy file. */
   const char* recovery;
} raw_run_t;

/**
 * @brief In-memory mzML produced by raw_write_mzml_mem().
 *
 * buf is malloc'd and NUL-terminated at buf[len]; pattern_detect() and
 * scan_mzml() rely on that terminator (mmap'd inputs get it from the
 * zero page past EOF, a malloc'd buffer must provide it). Free with free().
 */
typedef struct {
   char* buf;
   size_t len;
} mzml_buf_t;

/**
 * @brief Write `run` as a complete mzML 1.1.0 document into a growable
 *        memory buffer.
 * @param run The run to serialize.
 * @param salvage When nonzero, spectra whose fetch fails are skipped (with a
 *                summary warning) instead of aborting; an error is still
 *                returned when nothing at all is readable.
 * @param out Receives the malloc'd buffer and length on success.
 * @return 0 on success, -1 on allocation/provider failure.
 */
int raw_write_mzml_mem(const raw_run_t* run, int salvage, mzml_buf_t* out);

/**
 * @brief Write `run` as a complete mzML 1.1.0 document to an open file
 *        descriptor (memfd or temp file), one spectrum at a time.
 *
 * Used by the chunked (OOM-avoiding) path: the caller mmaps the fd
 * afterwards and runs the standard pipeline over the mapping.
 *
 * @param run The run to serialize.
 * @param salvage Skip unreadable spectra instead of failing (see
 *                raw_write_mzml_mem).
 * @param fd Open, writable file descriptor positioned at the start.
 * @return 0 on success, -1 on failure.
 */
int raw_write_mzml_fd(const raw_run_t* run, int salvage, int fd);

/**
 * @brief Estimate the mzML size a raw input will expand to.
 *
 * Peak data dominates: vendors store ~4-8 bytes/peak, the mzML carries
 * ~16 bytes/peak of base64 (12 raw bytes for double m/z + float intensity),
 * so ~3x the on-disk size is a conservative figure. XML metadata adds a
 * per-spectrum constant.
 *
 * @param raw_bytes Size of the vendor file (or the sum of a vendor
 *                  directory's files) on disk.
 * @param n_spectra Spectrum count of the run.
 * @return The estimate in bytes.
 */
uint64_t raw_estimate_mzml_bytes(uint64_t raw_bytes, uint64_t n_spectra);

/**
 * @brief Bytes currently available for allocation without swapping hard.
 *
 * Reads MemAvailable from /proc/meminfo on Linux; falls back to
 * sysconf(_SC_AVPHYS_PAGES) * page size elsewhere. Returns 0 when it cannot
 * be determined (callers then default to the chunked path).
 */
uint64_t raw_available_bytes(void);

/**
 * @brief OOM check deciding between the in-memory and the chunked path.
 * @param estimate raw_estimate_mzml_bytes() result.
 * @param available_bytes raw_available_bytes() result.
 * @return 1 when the mzML should be built in chunks (estimate exceeds 60% of
 *         what is available), 0 when it is safe to build in memory.
 */
int raw_should_chunk(uint64_t estimate, uint64_t available_bytes);

/**
 * @brief True when `path` looks like a supported vendor raw file by
 *        extension: .raw (Thermo file / Waters directory), .d (Bruker),
 *        .wiff, .wiff2, .t2d. Case-insensitive.
 */
int is_raw_vendor_path(const char* path);

/**
 * @brief Compress a vendor raw file directly to msz on `output_fd`.
 *
 * Loads libraw2ms_capi at runtime, converts run `run_index` of the file to
 * an mzML (in memory when the OOM check passes, memfd/temp-file chunked
 * otherwise) and runs the existing preprocess_mzml + compress_mzml pipeline.
 * All standard compression flags in `arguments` apply.
 *
 * @param input_path Path to the vendor file.
 * @param run_index Run to convert inside multi-run containers (.wiff).
 * @param salvage Skip unreadable spectra instead of failing (see
 *                raw_write_mzml_mem), and ask raw2ms to recover a corrupt
 *                file (truncated, zero-filled tail, shifted data blocks)
 *                instead of refusing it.
 * @param arguments Populated CLI arguments.
 * @param output_fd Open output .msz file descriptor.
 * @param recovered Out, may be NULL: set to 1 when the file was corrupt and
 *                  raw2ms recovered it (the caller names the output
 *                  *_recovered.msz), 0 otherwise.
 * @return 0 on success, -1 on failure (message already printed).
 */
int compress_raw(const char* input_path, uint64_t run_index, int salvage,
                 Arguments* arguments, int output_fd, int* recovered);

#endif /* RAW_INPUT_H */
