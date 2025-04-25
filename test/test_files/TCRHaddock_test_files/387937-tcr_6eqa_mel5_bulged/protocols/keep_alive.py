import os
import sys
import time

cmd = "touch " + sys.argv[1]
while 1:
    os.system(cmd)
    time.sleep(20)
