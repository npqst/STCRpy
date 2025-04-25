#!/bin/csh
#
# Check for awk 
#
set compresspdb=false
set found=`which awk |grep -v found |wc -l`
if ($found == 0) then
  echo 'awk not found'
  echo '==> no analysis possible'
  goto exit
endif
#
# Define the location of profit
#
if ( `printenv |grep PROFIT | wc -l` == 0) then
  set found=`which profit |wc -l`
  if ($found == 0) then
     echo 'PROFIT environment variable not defined'
     echo '==> no rmsd calculations '
  else
     setenv PROFIT `which profit`
  endif
endif
#
# extract HADDOCK score
#
foreach iclu (file.nam)
  echo "#struc haddock-score" >$iclu'_haddock-score'
  foreach i (`cat $iclu`)
    if ( -e $iclu.gz && ! -e $iclu ) then
      set compresspdb=true
      gunzip $iclu.gz
    endif
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_haddock-score'
    grep $i file.list | awk '{print $3}'>>$iclu'_haddock-score'
  end
end
#
# Extract buried surface area
#
foreach iclu (file.nam)
  echo "#struc bsa" >$iclu'_bsa'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_bsa'
    head -n 50 $i | grep -i buried | awk '{print $5}'>>$iclu'_bsa'
  end
end
#
# Extract energies
#
foreach iclu (file.nam)
  echo "#struc Einter  Enb Evdw+0.1Eelec Evdw  Eelec Eair Ecdih Ecoup Esani Evean Edani" >$iclu'_ener'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_ener'
    head -n 50 $i | grep -i energies | sed -e's/\,//g' | awk '{if(NF==15){print $3,$8+$9,$8+0.1*$9,$8,$9,$10,$11,$12,$13,$14,$15} else {print $3,$8+$9,$8+0.1*$9,$8,$9,$10,$11,$12,$13,$14,0.0}}'>>$iclu'_ener'
  end
end
#
# extract violations
#
foreach iclu (file.nam)
  echo "#struc #NOEviol #Dihedviol #Jviol #Saniviol #veanviol #Daniviol" >$iclu'_viol'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_viol'
    head -n 50 $i | grep -i " violation" | sed -e's/\,//g' | awk '{if (NF==8){print $3,$4,$5,$6,$7,$8} else {print $3,$4,$5,$6,$7,0}}'>>$iclu'_viol'
  end
end
#
# extract binding energy
#
foreach iclu (file.nam)
  echo "#struc dH" >$iclu'_dH'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_dH'
    head -n 50 $i | grep -i " Binding" | awk '{print $4}'>>$iclu'_dH'
  end
end
#
# extract desolvation energy
#
foreach iclu (file.nam)
  echo "#struc Edesolv" >$iclu'_Edesolv'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_Edesolv'
    head -n 50 $i | grep -i " desolv" | awk '{print $4}'>>$iclu'_Edesolv'
  end
end
if (-e l-RMSD.dat) then
   \cp l-RMSD.dat file.nam_rmsd
   set RMSD=true
else
#
# RMSD calculations only if PROFITDIR is defined
#
if (-e $PROFIT) then
   set RMSD=true
#
# calculate rmsd from lowest energy model (bakcbone -start and end residues)
#
set refe=`head -n 1 file.nam`
set atoms='CA,C2,P'
set refe=`head -1 file.nam`
set numCA=`grep -c ' CA ' $refe `
set numC2=`grep -c ' C2 ' $refe `
set numP=`grep -c ' P ' $refe `
set numatom=`echo $numCA $numC2 $numP | awk '{print $1+$2+$3}'`

if ($numatom < 3) then
  set atoms='*'
endif

$HADDOCKTOOLS/pdb_segid-to-chain $refe >refe.pdb

cat /dev/null > rmsd_best.disp
foreach iclu (file.nam)
  foreach i (`cat $iclu`)
    echo $i >>rmsd_best.disp
    $HADDOCKTOOLS/pdb_segid-to-chain $i >mobi.pdb
    $PROFIT <<_Eod_ |grep RMS |tail -1 >>rmsd_best.disp
      refe refe.pdb
      mobi mobi.pdb
      atom $atoms
      zone A*
      fit
      ratoms *
      rzone B*
      quit
_Eod_
  end
  echo ' ' >>rmsd_best.disp
  echo "#structure rmsd_all" >$iclu'_rmsd'
  awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else if  ($2 == "RMS:") {printf "%8.3f ",$3} else {printf "\n %s ",$1}}' rmsd_best.disp |grep pdb >>$iclu'_rmsd'
  echo ' ' >>$iclu'_rmsd'
 \rm rmsd_best.disp refe.pdb mobi.pdb
end
endif
endif
if ($RMSD == true) then
#
# Combine results in one file
# 
foreach iclu (file.nam)
  paste $iclu'_haddock-score' $iclu'_rmsd' $iclu'_ener' $iclu'_viol' $iclu'_bsa' $iclu'_dH' $iclu'_Edesolv'\
 | awk '{print $1,$2,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$18,$19,$20,$21,$22,$23,$25,$27,$29}' >structures_haddock-sorted.stat
end
grep "#struc" structures_haddock-sorted.stat > structures_rmsd-sorted.stat
sort -gk3 structures_haddock-sorted.stat |grep pdb >> structures_rmsd-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_ene-sorted.stat
sort -gk4 structures_haddock-sorted.stat |grep pdb >> structures_ene-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nb-sorted.stat
sort -gk5 structures_haddock-sorted.stat |grep pdb >> structures_nb-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nbw-sorted.stat
sort -gk6 structures_haddock-sorted.stat |grep pdb >> structures_nbw-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_air-sorted.stat
sort -gk9 structures_haddock-sorted.stat |grep pdb >> structures_air-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_airviol-sorted.stat
sort -gk15 structures_haddock-sorted.stat |grep pdb >> structures_airviol-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_bsa-sorted.stat
sort -rgk21 structures_haddock-sorted.stat |grep pdb >> structures_bsa-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_dH-sorted.stat
sort -gk22 structures_haddock-sorted.stat |grep pdb >> structures_dH-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_Edesolv-sorted.stat
sort -gk23 structures_haddock-sorted.stat |grep pdb >> structures_Edesolv-sorted.stat
#
# esle no RMSD calculation
#
else
#
# Combine results in one file
# 
foreach iclu (file.nam)
  paste $iclu'_haddock-score' $iclu'_ener' $iclu'_viol' $iclu'_bsa'  $iclu'_dH' $iclu'_Edesolv'\
 | awk '{print $1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$17,$18,$19,$20,$21,$23,$25,$27}' >structures_haddock-sorted.stat
end
grep "#struc" structures_haddock-sorted.stat > structures_ene-sorted.stat
sort -gk3 structures_haddock-sorted.stat |grep pdb >> structures_ene-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nb-sorted.stat
sort -gk4 structures_haddock-sorted.stat |grep pdb >> structures_nb-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nbw-sorted.stat
sort -gk5 structures_haddock-sorted.stat |grep pdb >> structures_nbw-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_air-sorted.stat
sort -gk8 structures_haddock-sorted.stat |grep pdb >> structures_air-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_airviol-sorted.stat
sort -gk14 structures_haddock-sorted.stat |grep pdb >> structures_airviol-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_bsa-sorted.stat
sort -rgk20 structures_haddock-sorted.stat |grep pdb >> structures_bsa-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_dH-sorted.stat
sort -gk21 structures_haddock-sorted.stat |grep pdb >> structures_dH-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_Edesolv-sorted.stat
sort -gk22 structures_haddock-sorted.stat |grep pdb >> structures_Edesolv-sorted.stat

endif

if ($compresspdb == true) then
  gzip *.pdb
endif

