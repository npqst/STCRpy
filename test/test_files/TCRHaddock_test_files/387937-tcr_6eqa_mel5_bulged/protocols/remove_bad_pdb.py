import os
import sys
import time

pdbfile = sys.argv[1]
pdbfile0 = pdbfile + "0"

time.sleep(5)
err = 0
pdblines = None
while err < 10:
    try:
        pdblines = open(pdbfile0).read()
        break
    except Exception:
        err += 1
        time.sleep(5)

ok = True
if pdblines:
    for line in pdblines.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            ok = False
            break
else:
    ok = False

if os.path.exists(pdbfile0):
    if ok:
        os.rename(pdbfile0, pdbfile)
    else:
        os.remove(pdbfile0)
