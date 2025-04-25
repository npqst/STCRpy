#!/bin/csh
#
if ($#argv < 1) goto usage
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

\rm tmp.pdb tmp1 >&/dev/null
foreach i (*.pdb)
    $HADDOCKTOOLS/pdb_blank_segid $i |pdb_chain |pdb_reres >tmp.pdb
    pales -bestFit -inD $1 -pdb tmp.pdb >$i:r.pales
    \rm tmp.pdb
end
goto exit

usage:


echo 'ana_pdb_Q-factor.csh: calculates Q-factor for residual dipolar couplings using pales'
echo ' '
echo ' usage: ana_Q-factor.csh pales_input_file'

exit:

