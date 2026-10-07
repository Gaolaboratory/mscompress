#!/usr/bin/env python3
"""End-to-end extraction gate: convert a fixture mzML -> msz -> extract, then
assert the extracted document is a newline-terminated indexedmzML whose
spectrum index really lands on <spectrum elements.

Usage: test_extract_index.py <mscompress> <workdir>

Leg 1 (always): test/data/test.mzML -> msz -> --extract. The fixture is an
indexed mzML, so extraction passes the stored document through verbatim: the
test asserts no double wrapping (exactly one <indexedmzML) and that the
source's own index block survived byte-for-byte. The fixture was made by
truncating a 431 MB run, so its inherited indexList is stale by construction
(41,608 offsets for 50 spectra, many past EOF); offset-landing is therefore
asserted on the raw leg only, where extraction builds a fresh true index.

Leg 2 (only when $MSC_RAW_FIXTURE points at a vendor raw): convert the raw to
msz via raw2ms, extract in full and filtered, and assert the full gate:
newline termination, single indexedmzML wrapper, every spectrum offset
seek-lands on b'<spectrum', PSI fileChecksum recomputes over the
<indexedmzML .. </indexList> span, sane isolation windows (low <= target <=
high in offset terms, width in (0, 1000], no zero offsets around a set
target), exactly one peak-mode label per spectrum, selected-ion m/z in every
precursor, and a pyteomics indexed parse that yields the same count as the
index. Needs $RAW2MS_LIBRARY as well; skipped cleanly (exit 0) when
$MSC_RAW_FIXTURE is unset so CI without vendor files still runs leg 1.
"""

import hashlib
import os
import re
import shutil
import subprocess
import sys

CENTROID = b'"MS:1000127"'  # centroid spectrum
PROFILE = b'"MS:1000128"'  # profile spectrum
SELECTED_ION_MZ = b'"MS:1000744"'  # selected ion m/z

FAILURES = []


def check(name, cond, detail=""):
    """Record and print one assertion; the process exits non-zero if any
    check failed."""
    if cond:
        print("  [ok] " + name)
    else:
        print("  [FAIL] " + name + (": " + detail if detail else ""))
        FAILURES.append(name)
    return cond


def info(msg):
    print("  " + msg)


def run(mscompress, args, env=None, name=None):
    """Run mscompress; on failure record it and print the captured output."""
    proc = subprocess.run(
        [mscompress] + args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
    )
    if proc.returncode != 0:
        out = proc.stdout.decode("utf-8", "replace")
        check(
            "mscompress %s exits 0" % (name or args[0]),
            False,
            "return code %d\n%s" % (proc.returncode, out[-2000:]),
        )
        return None
    return proc


def num_spectra(mscompress, msz):
    """Spectrum count recorded in the msz footer (--json -d)."""
    proc = subprocess.run(
        [mscompress, "--json", "-d", msz],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0:
        return None
    m = re.search(rb'"num_spectra":\s*(\d+)', proc.stdout)
    return int(m.group(1)) if m else None


def parse_spectrum_index(data):
    """(idRefs, offsets) of the spectrum <index>, or (None, None)."""
    m = re.search(rb'<index name="spectrum">(.*?)</index>', data, re.S)
    if not m:
        return None, None
    pairs = re.findall(rb'<offset idRef="([^"]*)">(\d+)</offset>', m.group(1))
    return [p[0] for p in pairs], [int(p[1]) for p in pairs]


def indexed_root_after_decl(data):
    """The bytes opening the document root that follows any <?xml
    declaration (leading whitespace tolerated), compared exactly."""
    i = 0
    n = len(data)
    while i < n and data[i:i+1] in (b" ", b"\t", b"\r", b"\n"):
        i += 1
    if data[i:i+5] == b"<?xml":
        j = data.find(b"?>", i)
        if j == -1:
            return None
        i = j + 2
        while i < n and data[i:i+1] in (b" ", b"\t", b"\r", b"\n"):
            i += 1
    return data[i:i + len(b"<indexedmzML")]


def common_document_checks(tag, data, n_body_spectra):
    """Newline + wrapper checks shared by both extraction modes."""
    check("%s: starts with <?xml declaration" % tag, data[:5] == b"<?xml",
          "first bytes are %r" % data[:20])
    head = indexed_root_after_decl(data)
    check("%s: document root is <indexedmzML" % tag,
          head == b"<indexedmzML", "root after decl is %r" % (head or b""))
    check("%s: exactly one <indexedmzML (no double wrap)" % tag,
          data.count(b"<indexedmzML") == 1,
          "found %d occurrences" % data.count(b"<indexedmzML"))
    check("%s: file ends with a newline" % tag, data[-1:] == b"\n",
          "last byte is %r" % data[-1:])
    newlines = data.count(b"\n")
    check("%s: newline-terminated multi-line mzML" % tag,
          newlines >= max(n_body_spectra, 1) and b"\n" in data[:200],
          "%d newlines total, none in first 200 bytes" % newlines
          if b"\n" not in data[:200] else "only %d newlines" % newlines)
    idrefs, offsets = parse_spectrum_index(data)
    check("%s: spectrum index present" % tag, offsets is not None and len(offsets) > 0,
          "no <index name=\"spectrum\"> with <offset> entries found")
    return idrefs, offsets


def check_fresh_index(tag, data, offsets, expected_count):
    """Wrap-mode assertions: true offsets + PSI fileChecksum (extraction
    built this index itself, so every value must be verifiable)."""
    bad = sum(1 for o in offsets if data[o:o + 9] != b"<spectrum")
    check("%s: every spectrum offset seek-lands on b'<spectrum'" % tag,
          bad == 0, "%d of %d offsets miss" % (bad, len(offsets)))
    check("%s: index size == extracted spectrum count" % tag,
          len(offsets) == expected_count,
          "%d offsets vs %d spectra" % (len(offsets), expected_count))
    m = re.search(rb"<fileChecksum>([0-9a-f]+)</fileChecksum>", data)
    if not check("%s: <fileChecksum> present" % tag, m is not None):
        return
    start = data.index(b"<indexedmzML")  # span starts at the wrapper itself
    end = data.rindex(b"</indexList>") + len(b"</indexList>")
    digest = hashlib.sha1(data[start:end]).hexdigest()
    check("%s: fileChecksum == SHA-1 over <indexedmzML .. </indexList>" % tag,
          digest == m.group(1).decode(),
          "recomputed %s, file says %s" % (digest, m.group(1).decode()))
    # <indexListOffset> is a required indexedmzML child: it must name the
    # byte position where the <indexList we just verified begins.
    ilo = re.search(rb"<indexListOffset>(\d+)</indexListOffset>", data)
    if check("%s: <indexListOffset> present" % tag, ilo is not None):
        check("%s: indexListOffset == position of <indexList" % tag,
              int(ilo.group(1)) == data.index(b"<indexList"),
              "file says %s, actual %d"
              % (ilo.group(1).decode(), data.index(b"<indexList")))


def has_msn(data):
    """Whether any spectrum is MS2 or deeper. An MS1-only run (a TSF MALDI
    imaging run, a survey-only file) has no isolation windows or precursors
    to check."""
    return any(int(v) >= 2 for v in
               re.findall(rb'name="ms level" value="([0-9]+)"', data))


def check_isolation_windows(tag, data):
    """T1 window rules over every isolationWindow in the document."""
    if not has_msn(data):
        info("%s: MS1-only document: isolation-window checks not applicable" % tag)
        return
    blocks = re.findall(rb"<isolationWindow>(.*?)</isolationWindow>", data, re.S)
    if not check("%s: isolation windows present" % tag, len(blocks) > 0):
        return
    bad_missing = bad_order = bad_width = bad_zero = 0
    first = None
    for b in blocks:
        t = re.search(rb'isolation window target m/z" value="([-0-9.eE+]+)"', b)
        lo = re.search(rb'isolation window lower offset" value="([-0-9.eE+]+)"', b)
        up = re.search(rb'isolation window upper offset" value="([-0-9.eE+]+)"', b)
        if not (t and lo and up):
            bad_missing += 1
            continue
        target, lower, upper = (float(t.group(1)), float(lo.group(1)),
                                float(up.group(1)))
        if first is None:
            first = (target, lower, upper)
        # low = target-lower <= target <= target+upper = high
        if lower < 0 or upper < 0:
            bad_order += 1
        if not (0 < lower + upper <= 1000):
            bad_width += 1
        if target != 0 and lower == 0 and upper == 0:
            bad_zero += 1  # the original redo bug signature
    check("%s: every window has target/lower/upper cvParams" % tag,
          bad_missing == 0, "%d of %d windows incomplete" % (bad_missing, len(blocks)))
    check("%s: every window has low <= target <= high" % tag,
          bad_order == 0, "%d of %d windows malformed" % (bad_order, len(blocks)))
    check("%s: every window width in (0, 1000] Da" % tag,
          bad_width == 0, "%d of %d windows out of range" % (bad_width, len(blocks)))
    check("%s: no set target bracketed by 0.0 offsets" % tag,
          bad_zero == 0, "%d windows carry the old zero-offset signature" % bad_zero)
    if first:
        info("first isolation window: target %.6f lower %.6f upper %.6f"
             % first)


def check_peak_mode_and_dia(tag, data, n_body_spectra):
    """C4 peak-mode labels and C5 DIA params (selected ion)."""
    centroid = data.count(CENTROID)
    profile = data.count(PROFILE)
    check("%s: exactly one peak-mode label per spectrum" % tag,
          centroid + profile == n_body_spectra,
          "%d centroid + %d profile vs %d spectra"
          % (centroid, profile, n_body_spectra))
    if not has_msn(data):
        info("%s: MS1-only document: precursor checks not applicable" % tag)
        return
    precursors = re.findall(rb"<precursorList.*?</precursorList>", data, re.S)
    if check("%s: precursorList elements present" % tag, len(precursors) > 0):
        bad_window = sum(1 for p in precursors if p.count(b"<isolationWindow>") != 1)
        check("%s: every precursor has exactly one isolation window" % tag,
              bad_window == 0, "%d of %d precursors" % (bad_window, len(precursors)))
        # DDA selects one ion per precursor; DIA (diaPASEF, SWATH) isolates a
        # window and selects none. Decided per document, so a DDA run that
        # loses a selected ion still fails here.
        with_ion = sum(1 for p in precursors
                       if b"selectedIonList" in p and SELECTED_ION_MZ in p)
        if with_ion == 0:
            info("%s: no selected ions: DIA document, windows only" % tag)
        else:
            check("%s: every precursor has a selected ion m/z (DDA)" % tag,
                  with_ion == len(precursors),
                  "%d of %d precursors lack one" % (len(precursors) - with_ion,
                                                   len(precursors)))


def pyteomics_count(path):
    from pyteomics import mzml  # deferred: only needed on legs that use it
    return sum(1 for _ in mzml.read(path, use_index=True))


def scan_numbers(idrefs):
    scans = []
    for idr in idrefs:
        m = re.search(rb"scan=(\d+)$", idr)
        if m:
            scans.append(int(m.group(1)))
    return scans


def leg_mzml(mscompress, workdir):
    print("leg 1: mzML path (repo fixture, indexed source -> passthrough)")
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "..", "test", "data", "test.mzML")
    src = os.path.abspath(src)
    if not check("leg 1: fixture %s exists" % src, os.path.isfile(src)):
        return
    with open(src, "rb") as f:
        source = f.read()
    src_index_at = source.find(b"<indexList count=")
    check("leg 1: fixture opens as indexed mzML (test expects passthrough)",
          src_index_at != -1 and b"<indexedmzML" in source[:200],
          "fixture layout changed; update this leg")

    msz = os.path.join(workdir, "fixture.msz")
    out = os.path.join(workdir, "fixture_full.mzML")
    if run(mscompress, [src, msz], name="compress") is None:
        return
    total = num_spectra(mscompress, msz)
    if not check("leg 1: msz reports num_spectra", total is not None and total > 0,
                 "got %r" % total):
        return
    if run(mscompress, ["--extract", "--extract-indices", "0-%d" % (total - 1),
                        msz, out], name="extract full") is None:
        return
    with open(out, "rb") as f:
        data = f.read()

    n_body = data.count(b"<spectrum ")
    src_body = source.count(b"<spectrum ")
    src_scans = [int(s) for s in
                 re.findall(rb'<spectrum [^>]*?scan=(\d+)', source)]
    common_document_checks("leg 1", data, n_body)
    check("leg 1: extracted every source spectrum", n_body == src_body,
          "%d extracted vs %d in source" % (n_body, src_body))
    out_scans = [int(s) for s in
                 re.findall(rb'<spectrum [^>]*?scan=(\d+)', data)]
    check("leg 1: extracted scan ids match the source, in order",
          out_scans == src_scans, "body scans differ from the source fixture")
    # The brief's no-double-wrap + surviving-index assertions for an
    # already-indexed source. Offsets themselves are the fixture's own
    # (stale by construction), so they are not seek-verified here.
    out_index_at = data.find(b"<indexList count=")
    check("leg 1: original <indexList survived extraction byte-for-byte",
          out_index_at != -1 and source[src_index_at:] == data[out_index_at:],
          "trailing index block differs from the source fixture")


def leg_raw(mscompress, workdir):
    fixture = os.environ.get("MSC_RAW_FIXTURE")
    if not fixture:
        info("MSC_RAW_FIXTURE unset: skipping the vendor raw leg "
             "(set it, plus RAW2MS_LIBRARY, to enable)")
        return
    print("leg 2: vendor raw path (%s)" % fixture)
    raw_lib = os.environ.get("RAW2MS_LIBRARY")
    if not check("leg 2: RAW2MS_LIBRARY is set", bool(raw_lib)):
        return
    # A vendor raw is a file (Thermo .raw, SCIEX .wiff/.wiff2/.t2d) or a
    # directory (Bruker .d, Waters .raw); isfile() alone refused every .d.
    if not check("leg 2: raw fixture exists",
                 os.path.isfile(fixture) or os.path.isdir(fixture), fixture):
        return
    if not check("leg 2: RAW2MS_LIBRARY file exists",
                 os.path.isfile(raw_lib), raw_lib):
        return
    env = dict(os.environ)
    env["RAW2MS_LIBRARY"] = raw_lib
    env["MSCOMPRESS_FORCE_CHUNKED"] = "1"

    msz = os.path.join(workdir, "fixture.msz")
    full = os.path.join(workdir, "fixture_full.mzML")
    filtered = os.path.join(workdir, "fixture_filtered.mzML")
    if run(mscompress, [fixture, msz], env=env, name="raw convert") is None:
        return
    total = num_spectra(mscompress, msz)
    if not check("leg 2: msz reports num_spectra", total is not None and total > 0,
                 "got %r" % total):
        return
    info("num_spectra: %d" % total)
    if run(mscompress, ["--extract", "--extract-indices", "0-%d" % (total - 1),
                        msz, full], name="raw extract full") is None:
        return
    with open(full, "rb") as f:
        data = f.read()
    n_body = data.count(b"<spectrum ")
    idrefs, offsets = common_document_checks("leg 2 full", data, n_body)
    if offsets:
        check_fresh_index("leg 2 full", data, offsets, n_body)
        check("leg 2 full: body spectrum count == msz num_spectra",
              n_body == total, "%d vs %d" % (n_body, total))
        check_isolation_windows("leg 2 full", data)
        check_peak_mode_and_dia("leg 2 full", data, n_body)

    # Filtered extraction mirrors the manual smoke: --extract-scans 2-25 —
    # where scan numbers mean that. timsTOF ids are "frame=N scan=M" with M a
    # mobility scan, SCIEX ids count cycles, so a scan range selects nothing
    # meaningful there; those fixtures select by position instead.
    scan_ids = bool(idrefs) and all(
        b"scan=" in i and b"frame=" not in i for i in idrefs)
    if scan_ids:
        sel, how = ["--extract-scans", "2-25"], "scans 2-25"
    else:
        sel, how = ["--extract-indices", "1-24"], "indices 1-24"
    if run(mscompress, ["--extract"] + sel + [msz, filtered],
           name="raw extract %s" % how) is None:
        return
    with open(filtered, "rb") as f:
        fdata = f.read()
    fn_body = fdata.count(b"<spectrum ")
    check("leg 2 filtered: the %s selection matched spectra" % how, fn_body > 0,
          "extracted an empty document")
    fidrefs, foffsets = common_document_checks("leg 2 filtered", fdata, fn_body)
    if foffsets:
        check_fresh_index("leg 2 filtered", fdata, foffsets, fn_body)
        if scan_ids:
            wanted = [s for s in range(2, 26) if s in set(scan_numbers(idrefs or []))]
            check("leg 2 filtered: exactly scans 2-25 extracted",
                  sorted(scan_numbers(fidrefs)) == wanted,
                  "got %s, expected %s" % (sorted(scan_numbers(fidrefs))[:8], wanted[:8]))
        else:
            wanted = list((idrefs or [])[1:25])
            check("leg 2 filtered: exactly indices 1-24 extracted",
                  list(fidrefs) == wanted,
                  "got %s, expected %s" % (list(fidrefs)[:3], wanted[:3]))
        check_isolation_windows("leg 2 filtered", fdata)
        check_peak_mode_and_dia("leg 2 filtered", fdata, fn_body)

    try:
        n_py = pyteomics_count(full)
        check("leg 2: pyteomics indexed parse count == index size",
              n_py == len(offsets or []), "pyteomics %d vs %d offsets"
              % (n_py, len(offsets or [])))
        n_pyf = pyteomics_count(filtered)
        check("leg 2: pyteomics parses the filtered document",
              n_pyf == fn_body, "pyteomics %d vs %d spectra" % (n_pyf, fn_body))
    except Exception as e:  # noqa: BLE001 - report, do not crash the gate
        check("leg 2: pyteomics indexed parse", False,
              "%s: %s (pyteomics needs the psims package for mzML)"
              % (type(e).__name__, e))


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    mscompress, workdir = os.path.abspath(sys.argv[1]), sys.argv[2]
    if not check("mscompress binary exists", os.path.isfile(mscompress),
                 mscompress):
        return 1
    if os.path.isdir(workdir):
        shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir, exist_ok=True)

    leg_mzml(mscompress, workdir)
    leg_raw(mscompress, workdir)

    if FAILURES:
        print("FAILED (%d): %s" % (len(FAILURES), ", ".join(FAILURES)))
        print("artifacts kept in %s" % workdir)
        return 1
    print("all extraction gate checks passed")
    shutil.rmtree(workdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
