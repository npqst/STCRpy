#!/bin/csh
#
# ana_noe_viol: analysis the output of the xplor script accept.inp
#               for NOE violations with as output for each violated
#               NOE restraints, the average distance, average violation
#               and number of time it is violated in the ensemble of 
#               structures
# ana_noe_viol.awk and count_noe_viol.awk are needed
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
sed -e 's/=1/=\ 1/g' -e 's/=2/=\ 2/g' -e 's/=3/=\ 3/g' $1 | awk -f ana_noe_viol.awk |sort -k2 -k1 >ana_noe_tmp1
echo 'DONE' >>ana_noe_tmp1
if (`wc -l ana_noe_tmp1 |awk '{print $1}'` == 1) then
  \rm ana_noe_tmp1
  goto exit
endif
awk -f count_noe_viol.awk ana_noe_tmp1 | sort -n -k6 
\rm ana_noe_tmp1
goto exit

usage:
echo 'Usage: ana_noe_viol accept.out'
echo ' '
echo '       accept.out contains the NOE violations analysis'
echo '                  from the X-Plor script accept.inp '

exit:
