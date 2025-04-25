#!/bin/csh
#
# ana_dihed_viol: analysis the output of the xplor script accept.inp
#               for diehdral angle violations with as output for each violated
#               restraints, the average dihedral, average violation
#               and number of time it is violated in the ensemble of 
#               structures
# ana_dihed_viol.awk and count_dihed_viol.awk are needed
#
# A. Bonvin, Utrecht University Nov. 99
#
set found=`which awk |grep -v found |wc -l`
if ($found == 0) then
  echo 'awk not found'
  echo '==> no NOE violation analysis'
  goto exit
endif

if ($#argv < 1) goto usage

if (`wc -l $1 |awk '{print $1}'` == 0) goto exit

egrep -B 4 -A 0 "Dihedral=" $1 | awk '{if (NF > 1 && $1 != "Dihedral=") {printf "%4s %4s %4s ",$1,$2,$3} else if (NF >1 && $1 == "Dihedral=") {print $0} }' | sort -k2 -n >ana_dihed_tmp1

echo 'DONE' >>ana_dihed_tmp1
if (`wc -l ana_dihed_tmp1 |awk '{print $1}'` == 1) then
  \rm ana_dihed_tmp1
  goto exit
endif
awk -f count_dihed_viol.awk ana_dihed_tmp1 | sort -n -k6 
\rm ana_dihed_tmp1
goto exit

usage:
echo 'Usage: ana_dihed_viol accept.out'
echo ' '
echo '       accept.out contains the dihedral angle violations analysis'
echo '                  from the X-Plor script accept.inp '

exit:
