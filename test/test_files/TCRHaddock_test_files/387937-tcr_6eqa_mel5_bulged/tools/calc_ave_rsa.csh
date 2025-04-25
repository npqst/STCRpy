#!/bin/tcsh
#
if ( `printenv |grep NACCESS | wc -l` == 0) then
  set found=`which naccess |wc -l`
  if ($found == 0) then
     echo 'NACCESS environment variable not defined'
     echo '==> exiting'
     goto exit
  else
     setenv NACCESS `which naccess`
  endif
endif
#
# Check for awk 
#
set found=`which awk |grep -v found |wc -l`
if ($found == 0) then
  echo 'awk or gawk not found'
  echo '==> no NOE violation analysis'
  goto exit
endif

foreach i ($argv)
  cp $i .
  $NACCESS `basename $i`
end

set rsaname=`ls -al *.rsa|head -2 |tail -1 |awk '{if (NF==10) {print $10} else if (NF==9) {print $9} else {print $8}}'`
echo $rsaname
#set ires=`grep RES $rsaname | grep -v atoms | tail -1 | awk '{if (NF==14){print $4} else {print $3}}'`
set ires=`grep RES $rsaname | grep -v atoms | wc -l |awk '{print $1}'`
#set i=`grep RES $rsaname | grep -v atoms | head -1 | awk '{if (NF==14){print $4} else {print $3}}'`
#@ i-=1
set i=0

\rm rsa_ave.lis
echo '# resnam resnum <rsa_all> (sd) <rsa_back> (sd) <rsa_side> (sd)' >rsa_ave.lis
while ($i < $ires)
 @ i+=1
 cat /dev/null > tmpall
 cat /dev/null > tmpback
 cat /dev/null > tmpside
 foreach str (*.rsa)
   grep RES $str |grep -v atoms | head -$i |tail -1 |awk '{if (NF == 13) {print $2,$3,$5}; if (NF == 14){print $2,$3,$4,$6}}' >>tmpall
   grep RES $str |grep -v atoms | head -$i |tail -1 |awk '{if (NF == 13) {print $2,$3,$9}; if (NF == 14){print $2,$3,$4,$10}}' >>tmpback
   grep RES $str |grep -v atoms | head -$i |tail -1 |awk '{if (NF == 13) {print $2,$3,$7}; if (NF == 14){print $2,$3,$4,$8}}' >>tmpside
 end
 head -1 tmpall | awk '{if (NF==3) {printf "%4s %6d",$1,$2} else {printf "%4s %2s %6d",$1,$2,$3}}' >>rsa_ave.lis
 awk '{if (NF==4) {print $4} else {print $3}}' tmpall | $HADDOCKTOOLS/average.sh | awk '{printf "%8.3f %8.3f ",$1,$2}' >> rsa_ave.lis
 awk '{if (NF==4) {print $4} else {print $3}}' tmpback | $HADDOCKTOOLS/average.sh | awk '{printf "%8.3f %8.3f ",$1,$2}' >> rsa_ave.lis
 awk '{if (NF==4) {print $4} else {print $3}}' tmpside | $HADDOCKTOOLS/average.sh | awk '{printf "%8.3f %8.3f\n",$1,$2}' >> rsa_ave.lis
# set oname=`head -1 tmpall | awk '{if (NF==4){print $1"_"$2"_"$3".rsa_all"} else {print $1"_"$2".rsa_all"}}' `
# set onameb=`head -1 tmpback | awk '{if (NF==4){print $1"_"$2"_"$3".rsa_back"} else {print $1"_"$2".rsa_back"}}' `
# set onames=`head -1 tmpside | awk '{if (NF==4){print $1"_"$2"_"$3".rsa_side"} else {print $1"_"$2".rsa_side"}}' `
# \mv tmpall $oname
# \mv tmpback $onameb
# \mv tmpside $onames
end

exit:
