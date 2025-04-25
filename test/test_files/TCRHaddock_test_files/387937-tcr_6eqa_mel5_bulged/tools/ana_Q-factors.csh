#!/bin/csh
#
if ($#argv < 1) goto usage
#
# Check for awk 
#
set found=`which awk |grep -v found |wc -l`
if ($found == 0) then
  echo 'awk not found'
  echo '==> no NOE violation analysis'
  goto exit
endif
#
# Define the location of pales
#
if ( `printenv |grep PALES | wc -l` == 0) then
  set found=`which pales |wc -l`
  if ($found == 0) then
     echo 'PALES environment variable not defined'
     echo '==> stopping'
     goto exit
  else
     setenv PALES `which pales`
  endif
endif
#
\rm tmp.pdb tmp1 >&/dev/null
foreach iclu (file.nam)
  echo "#struc Q-factor Corr-R" >$iclu'_Q-factor'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_Q-factor'
    $HADDOCKTOOLS/pdb_blank_segid $i >tmp.pdb
    $PALES -bestFit -inD $1 -pdb tmp.pdb >tmp.qfactor
    grep CORNILESCU tmp.qfactor |awk '{printf "%6.3f ",$4}' >>$iclu'_Q-factor'
    grep "CORR R" tmp.qfactor |awk '{print $4}' >>$iclu'_Q-factor'
    \rm tmp.pdb
    \mv tmp.qfactor $i:r.pales
  end
end
goto exit

usage:


echo 'ana_Q-factor.csh: calculates Q-factor for residual dipolar couplings using pales'
echo ' '
echo ' usage: ana_Q-factor.csh pales_input_file'

exit:

