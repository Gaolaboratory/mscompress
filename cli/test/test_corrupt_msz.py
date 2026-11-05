#!/usr/bin/env python3
"""Corrupt .msz files must be rejected with an error, never crash.

usage: test_corrupt_msz.py <mscompress> <valid.msz> <tmpdir>
"""
import os, subprocess, sys

exe, good, tmp = sys.argv[1:4]
os.makedirs(tmp, exist_ok=True)
data = open(good, "rb").read()
cases = {
    # bytes removed from the middle: footer at EOF is intact, offsets are stale
    "cut_middle": data[:600] + data[1600:],
    # stale tail: a whole second copy appended after the footer region
    "doubled": data + data[512:],
    # truncated: footer lost
    "truncated": data[: len(data) // 2],
    # just the header
    "header_only": data[:512],
}

def run(path):
    out = os.path.join(tmp, "out.mzML")
    r = subprocess.run([exe, path, out], capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr

rc, _ = run(good)
assert rc == 0, f"valid file failed to decompress (rc={rc})"

failed = False
for name, blob in cases.items():
    p = os.path.join(tmp, name + ".msz")
    open(p, "wb").write(blob)
    rc, msg = run(p)
    if rc < 0 or rc >= 128:
        print(f"FAIL {name}: crashed (rc={rc})"); failed = True
    elif rc == 0:
        print(f"FAIL {name}: corrupt file was accepted"); failed = True
    else:
        print(f"ok   {name}: rc={rc}")
sys.exit(1 if failed else 0)
