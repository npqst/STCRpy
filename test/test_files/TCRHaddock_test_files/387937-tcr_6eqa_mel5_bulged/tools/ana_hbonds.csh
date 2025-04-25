#!/bin/csh
#
# ana_hbonds.csh: analysis the output of the CNS script print_hbonds.inp
#                 for hbond analysis with as output for each detected 
#                 H-bond the average distance and the number of times detected
# count_hbonds.awk is needed
#
# A. Bonvin, Utrecht University Feb. 03

if ($#argv < 1) goto usage

awk -f count_hbonds.awk $1 | sort -k7
goto exit

usage:
echo 'Usage: ana_hbonds.csh hbond.disp'

exit:

