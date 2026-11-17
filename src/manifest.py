#!/usr/bin/env python3
"""Write MANIFEST.sha256 for the release: every file under the given directories (sorted), one line
'<sha256>  <relative path>'.  usage: manifest.py OUT DIR [DIR ...]   (run from the release root);
    manifest.py --check MANIFEST   verifies every listed file."""
import hashlib, os, sys


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if __name__ == "__main__":
    if sys.argv[1] == "--check":
        bad = 0
        n = 0
        for line in open(sys.argv[2]):
            h, p = line.rstrip("\n").split("  ", 1)
            n += 1
            if not os.path.exists(p) or sha(p) != h:
                print("MISMATCH/MISSING:", p)
                bad += 1
        print(f"manifest check: {n} files, {bad} problems")
        sys.exit(1 if bad else 0)
    out, dirs = sys.argv[1], sys.argv[2:]
    files = []
    for d in dirs:
        for root, _, names in os.walk(d):
            for nm in names:
                p = os.path.join(root, nm)
                if "__pycache__" in p or p.endswith(".pyc") or os.path.islink(p):
                    continue
                files.append(p)
    with open(out, "w") as fo:
        for p in sorted(files):
            fo.write(f"{sha(p)}  {p}\n")
    print(f"wrote {out}: {len(files)} files")
