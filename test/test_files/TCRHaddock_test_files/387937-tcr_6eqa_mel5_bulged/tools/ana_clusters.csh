#!/bin/csh
if ($#argv < 1) goto usage

if ($#argv > 1) then
  if ($#argv != 3) goto usage
  if (`echo $argv[1]| sed -e "1,$ s/-//"` != "best") goto usage
  set nmb=$2
  set clusfile=$3
  set best=on
else
  set clusfile=$1
  set best=off
endif

set compresspdb=false
#
# Check for awk 
#
set found=`which awk |grep -v found |wc -l`
if ($found == 0) then
  echo 'awk not found'
  echo '==> analysis not possible'
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
# Start the analysis of the clusters
#
set numclu=`wc -l $clusfile | awk '{print $1}'`

set iclu=0

while ($iclu < $numclu)

  @ iclu += 1
  set clu=`head -n $iclu $clusfile |tail -1 | awk '{$1=" ";$2=" ";$3=" "; print $0}'`
  cat /dev/null >file.list_clust$iclu
  cat /dev/null >file.list_tmp
  set num=0
  foreach i ($clu)
    @ num +=1
    head -n $i file.list |tail -1 >>file.list_tmp
  end
  sort -gk3 file.list_tmp >file.list_clust$iclu

end

set iclu=0

while ($iclu < $numclu)

  @ iclu += 1
  sed -e 's/"//g' -e 's/PREVIT\://' file.list_clust$iclu| awk '{print $1}' >file.list_tmp
  set pdbs=`cat file.list_tmp`
  cat /dev/null >file.nam_clust$iclu
  head -n 8 file.cns >file.cns_clust$iclu
  set num=0
  foreach i ($pdbs)
    if ( -e $i.gz && ! -e $i ) then
      set compresspdb=true
      gunzip $i.gz
    endif
    @ num +=1 
    echo $i >>file.nam_clust$iclu
    echo 'evaluate (&filenames.bestfile_'$num'="PREVIT:'$i'")'>>file.cns_clust$iclu
  end

end
#
# calculate average cluster energies
#
echo "#Cluster Nstruc Einter sd Enb sd Evdw+0.1Eelec sd Evdw sd Eelec sd Eair sd Ecdih sd Ecoup sd Esani sd Evean sd Edani sd" >cluster_ener.txt
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo "#struc Einter Enb Evdw+0.1Eelec Evdw Eelec Eair Ecdih Ecoup Esani Evean Edani" >$iclu'_ener'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_ener'
    head -n 50 $i | grep -i energies | sed -e's/\,//g' | awk '{print $3,$8+$9,$8+0.1*$9,$8,$9,$10,$11,$12,$13,$14,$15}'>>$iclu'_ener'
  end
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_ener.txt
  grep pdb $iclu | wc -l |awk '{printf "%5d ",$1}' >>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $2}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $3}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $4}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $5}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $6}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $7}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $8}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $9}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $10}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $11}' |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt
  grep -v '#' $iclu'_ener' | awk '{print $12}' |$HADDOCKTOOLS/average.sh >>cluster_ener.txt
end
#
# calculate average cluster violations
#
echo "#Cluster #AIRviol sd #dihedviol sd #Coupviol sd #Saniviol sd #Veanviol sd #Daniviol sd" >cluster_viol.txt
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo "#struc #NOEviol #Dihedviol #Coupviol #Saniviol #Veanviol #Daniviol" >$iclu'_viol'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_viol'
    head -n 50 $i |grep -i " violation" | sed -e's/\,//g' | awk '{print $3,$4,$5,$6,$7,$8}'>>$iclu'_viol'
  end
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_viol.txt
  grep -v "#" $iclu'_viol' | awk '{print $2}' |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt
  grep -v "#" $iclu'_viol' | awk '{print $3}' |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt
  grep -v "#" $iclu'_viol' | awk '{print $4}' |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt
  grep -v "#" $iclu'_viol' | awk '{print $5}' |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt
  grep -v "#" $iclu'_viol' | awk '{print $6}' |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt
  grep -v "#" $iclu'_viol' | awk '{print $7}' |$HADDOCKTOOLS/average.sh >>cluster_viol.txt
end
#
# calculate average HADDOCK score
#
echo "#Cluster haddock-score sd" >cluster_haddock-score.txt
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo "#struc haddock-score" >$iclu'_haddock-score'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_haddock-score'
    grep $i file.list | awk '{print $3}'>>$iclu'_haddock-score'
  end
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_haddock-score.txt
  grep -v "#" $iclu'_haddock-score' | awk '{print $2}' |$HADDOCKTOOLS/average.sh >>cluster_haddock-score.txt
end
#
# calculate average cluster buried surface area
#
echo "#Cluster BSA sd" >cluster_bsa.txt
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo "#struc BSA" >$iclu'_bsa'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_bsa'
    head -n 50 $i |grep -i buried | awk '{print $5}'>>$iclu'_bsa'
  end
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_bsa.txt
  grep -v "#" $iclu'_bsa' | awk '{print $2}' |$HADDOCKTOOLS/average.sh >>cluster_bsa.txt
end
#
# calculate average cluster binding energy
#
echo "#Cluster #dH sd" >cluster_dH.txt
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo "#struc dH" >$iclu'_dH'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_dH'
    head -n 50 $i | grep -i " Binding" | awk '{print $4}'>>$iclu'_dH'
  end
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_dH.txt
  grep -v "#" $iclu'_dH' | awk '{print $2}' |$HADDOCKTOOLS/average.sh >>cluster_dH.txt
end
#
# calculate average cluster desolvation energy
#
echo "#Cluster #Edesolv sd" >cluster_Edesolv.txt
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo "#struc Edesolv" >$iclu'_Edesolv'
  foreach i (`cat $iclu`)
    echo $i |awk '{printf "%s ",$1}'>>$iclu'_Edesolv'
    head -n 50 $i | grep -i " desolv" | awk '{print $4}'>>$iclu'_Edesolv'
  end
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_Edesolv.txt
  grep -v "#" $iclu'_Edesolv' | awk '{print $2}' |$HADDOCKTOOLS/average.sh >>cluster_Edesolv.txt
end
#
# RMSD calculations only if PROFITDIR is defined
#
if ($?PROFIT) then
#
# calculate average cluster rmsd from lowest energy structure inside the cluster
#
set atoms='CA,C2,P'
set refe=`head -1 file.nam`
set numCA=`grep -c ' CA ' $refe `
set numC2=`grep -c ' C2 ' $refe `
set numP=`grep -c ' P ' $refe `
set numatom=`echo $numCA $numC2 $numP | awk '{print $1+$2+$3}'`

if ($numatom < 3) then
  set atoms='*'
endif

echo "#Cluster rmsd sd" >cluster_rmsd.txt
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  cat /dev/null >rmsd_best.disp
  set ifil=`echo $iclu |sed -e 's/nam/list/'`
  set refe=`sort -gk3 $ifil |head -n 1 | awk '{print $1}' |sed -e 's/"PREVIT://' -e 's/"//'`
  foreach i (`cat $iclu`)
    echo $i >>rmsd_best.disp
    $PROFIT <<_Eod_ |grep RMS >>rmsd_best.disp
      refe $refe
      mobi $i
      atom $atoms
      fit
      quit
_Eod_
  end
  echo "#Structure rmsd_all" >$iclu'_rmsd'
  awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' rmsd_best.disp \
      | grep pdb >>$iclu'_rmsd'
  \rm rmsd_best.disp
  echo $iclu |sed -e 's/list/nam/' |awk '{printf "%s ",$1}' >>cluster_rmsd.txt
  grep -v "#" $iclu'_rmsd' | awk '{if (NF == 2) print $2}' |$HADDOCKTOOLS/average.sh >>cluster_rmsd.txt
end
#
# if file present, calculate average cluster i-RMSD from target
#
if (-e i-RMSD.dat) then
  echo "#Cluster i-RMSD sd" >cluster_rmsd-Emin.txt
  foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9])
    echo "#Structure i-RMSD" >$iclu'_rmsd-Emin'
    foreach i (`cat $iclu`)
      grep $i i-RMSD.dat >>$iclu'_rmsd-Emin'
    end
    echo $iclu |sed -e 's/list/nam/' |awk '{printf "%s ",$1}' >>cluster_rmsd-Emin.txt
    grep -v "#" $iclu'_rmsd-Emin' | awk '{if (NF == 2) print $2}' |$HADDOCKTOOLS/average.sh >>cluster_rmsd-Emin.txt
  end
else
#
# calculate average cluster rmsd from lowest energy structure
#
  set atoms='CA,C2,P'
  set refe=`head -1 file.nam`
  set numCA=`grep -c ' CA ' $refe `
  set numC2=`grep -c ' C2 ' $refe `
  set numP=`grep -c ' P ' $refe `
  set numatom=`echo $numCA $numC2 $numP | awk '{print $1+$2+$3}'`

  if ($numatom < 3) then
    set atoms='*'
  endif

  echo "#Cluster rmsd-Emin sd" >cluster_rmsd-Emin.txt
  foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
    cat /dev/null >rmsd_best-Emin.disp
    set refe=`head -n 1 file.nam`
    foreach i (`cat $iclu`)
      echo $i >>rmsd_best-Emin.disp
      $PROFIT <<_Eod_ |grep RMS >>rmsd_best-Emin.disp
        refe $refe
        mobi $i
        atom $atoms
        fit
        quit
_Eod_
    end
    echo "#Structure rmsd_Emin" >$iclu'_rmsd-Emin'
    awk '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' rmsd_best-Emin.disp \
        | grep pdb >>$iclu'_rmsd-Emin'
    \rm rmsd_best-Emin.disp
    echo $iclu |sed -e 's/list/nam/' |awk '{printf "%s ",$1}' >>cluster_rmsd-Emin.txt
    grep -v "#" $iclu'_rmsd-Emin' | awk '{if (NF == 2) print $2}' |$HADDOCKTOOLS/average.sh >>cluster_rmsd-Emin.txt
  end  
endif
#
# Combine cluster results in one file
# 
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  paste $iclu'_rmsd' $iclu'_rmsd-Emin' $iclu'_ener' $iclu'_viol' $iclu'_bsa' \
 | awk '{print $1,$2,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$18,$19,$20,$22,$23,$25}' >$iclu.stat
end
#
# Combine results in one file
# 
paste cluster_haddock-score.txt cluster_rmsd.txt cluster_rmsd-Emin.txt cluster_ener.txt cluster_viol.txt cluster_bsa.txt cluster_dH.txt cluster_Edesolv.txt \
 | awk '{print $1,$2,$3,$5,$6,$8,$9,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$48,$49,$51,$52,$54,$55}' >clusters.stat
head -n 1 clusters.stat >clusters_haddock-sorted.stat
sort -gk2 clusters.stat |grep file >>clusters_haddock-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_ene-sorted.stat
sort -gk9 clusters_haddock-sorted.stat |grep file >>clusters_ene-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_nb-sorted.stat
sort -gk11 clusters_haddock-sorted.stat |grep file >>clusters_nb-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_nbw-sorted.stat
sort -gk13 clusters_haddock-sorted.stat |grep file >>clusters_nbw-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_air-sorted.stat
sort -gk19 clusters_haddock-sorted.stat |grep file >>clusters_air-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_sani-sorted.stat
sort -gk25 clusters_haddock-sorted.stat |grep file >>clusters_sani-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_vean-sorted.stat
sort -gk27 clusters_haddock-sorted.stat |grep file >>clusters_vean-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_dani-sorted.stat
sort -gk29 clusters_haddock-sorted.stat |grep file >>clusters_dani-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_bsa-sorted.stat
sort -grk43 clusters_haddock-sorted.stat |grep file >>clusters_bsa-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_dH-sorted.stat
sort -gk45 clusters_haddock-sorted.stat |grep file >>clusters_dH-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_Edesolv-sorted.stat
sort -gk47 clusters_haddock-sorted.stat |grep file >>clusters_Edesolv-sorted.stat

else
#
# Combine cluster results in one file
# 
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  paste $iclu'_ener' $iclu'_viol' $iclu'_bsa' \
 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14,$15,$16,$17,$18,$19,$21}' >$iclu.stat
end
#
# Combine results in one file
# 
paste cluster_haddock-score.txt cluster_ener.txt cluster_viol.txt cluster_bsa.txt cluster_dH.txt cluster_Edesolv.txt \
| awk '{print $1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$42,$43,$45,$46,$48,$49,$51,$52}' >clusters.stat
head -n 1 clusters.stat >clusters_haddock-sorted.stat
sort -gk2 clusters.stat |grep file >>clusters_haddock-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_ene-sorted.stat
sort -gk5 clusters_haddock-sorted.stat |grep file >>clusters_ene-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_nb-sorted.stat
sort -gk7 clusters_haddock-sorted.stat |grep file >>clusters_nb-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_nbw-sorted.stat
sort -gk11 clusters_haddock-sorted.stat |grep file >>clusters_nbw-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_air-sorted.stat
sort -gk15 clusters_haddock-sorted.stat |grep file >>clusters_air-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_sani-sorted.stat
sort -gk21 clusters_haddock-sorted.stat |grep file >>clusters_sani-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_vean-sorted.stat
sort -gk23 clusters_haddock-sorted.stat |grep file >>clusters_vean-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_dani-sorted.stat
sort -gk25 clusters_haddock-sorted.stat |grep file >>clusters_dani-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_bsa-sorted.stat
sort -grk39 clusters_haddock-sorted.stat |grep file >>clusters_bsa-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_dH-sorted.stat
sort -gk41 clusters_haddock-sorted.stat |grep file >>clusters_dH-sorted.stat
head -n 1 clusters_haddock-sorted.stat >clusters_Edesolv-sorted.stat
sort -gk43 clusters_haddock-sorted.stat |grep file >>clusters_Edesolv-sorted.stat


endif

#
# Now repeat for best nmb structures if required
#
if ($best == "on") then
#
# generate file.list file.nam and file.cns
#
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  head -n $nmb $iclu >$iclu'_best'$nmb
  echo ' ' >>$iclu'_best'$nmb
end
foreach iclu (file.list_clust[0-9] file.list_clust[0-9][0-9] file.list_clust[0-9][0-9][0-9])
  head -n $nmb $iclu >$iclu'_best'$nmb
  echo ' ' >>$iclu'_best'$nmb
end
set totline=$nmb
@ totline+=8
foreach iclu (file.cns_clust[0-9] file.cns_clust[0-9][0-9] file.cns_clust[0-9][0-9][0-9])
  head -n $totline $iclu >$iclu'_best'$nmb
  echo ' ' >>$iclu'_best'$nmb
end
#
# calculate average cluster energies on best nmb structures
#
echo "#Cluster Nstruc Einter sd Enb sd Evdw+0.1Eelec sd Evdw sd Eelec sd Eair sd Ecdih sd Ecoup sd Esani sd Evean sd Edani sd" >cluster_ener.txt_best$nmb
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_ener.txt_best$nmb
  grep pdb $iclu | wc -l |awk '{printf "%5d ",$1}' >>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $3}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $4}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $5}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $6}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $7}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $8}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $9}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $10}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $11}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
  | awk '{printf "%8.2f %8.2f ",$1,$2}'>>cluster_ener.txt_best$nmb
  grep -v '#' $iclu'_ener' | awk '{print $12}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_ener.txt_best$nmb
end
#
# calculate average cluster violations
#
echo "#Cluster #AIRviol sd #dihedviol sd #Coupviol sd #Saniviol sd #Veanviol sd #Daniviol sd" >cluster_viol.txt_best$nmb
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_viol.txt_best$nmb
  grep -v "#" $iclu'_viol' | awk '{print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt_best$nmb
  grep -v "#" $iclu'_viol' | awk '{print $3}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt_best$nmb
  grep -v "#" $iclu'_viol' | awk '{print $4}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt_best$nmb
  grep -v "#" $iclu'_viol' | awk '{print $5}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt_best$nmb
  grep -v "#" $iclu'_viol' | awk '{print $6}' | head -n $nmb |$HADDOCKTOOLS/average.sh \
   | awk '{printf "%8.2f %8.2f ",$1,$2}' >>cluster_viol.txt_best$nmb
  grep -v "#" $iclu'_viol' | awk '{print $7}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_viol.txt_best$nmb
end
#
# calculate average cluster buried surface area
#
echo "#Cluster BSA sd" >cluster_bsa.txt_best$nmb
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_bsa.txt_best$nmb
  grep -v "#" $iclu'_bsa' | awk '{print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_bsa.txt_best$nmb
end
#
# calculate average HADDOCK score
#
echo "#Cluster haddock-score sd" >cluster_haddock-score.txt_best$nmb
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_haddock-score.txt_best$nmb
  grep -v "#" $iclu'_haddock-score' | awk '{print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_haddock-score.txt_best$nmb
end
#
# calculate average cluster binding energy
#
echo "#Cluster #dH sd" >cluster_dH.txt_best$nmb
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_dH.txt_best$nmb
  grep -v "#" $iclu'_dH' | awk '{print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_dH.txt_best$nmb
end
#
# calculate average cluster desolvation energy
#
echo "#Cluster #Edesolv sd" >cluster_Edesolv.txt_best$nmb
foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |awk '{printf "%s ",$1}' >>cluster_Edesolv.txt_best$nmb
  grep -v "#" $iclu'_Edesolv' | awk '{print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_Edesolv.txt_best$nmb
end
#
# RMSD calculations only if PROFITDIR is defined
#
if (-e $PROFIT) then
#
# calculate average cluster rmsd from lowest energy structure inside the cluster
#
echo "#Cluster rmsd sd" >cluster_rmsd.txt_best$nmb

foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |sed -e 's/list/nam/' |awk '{printf "%s ",$1}' >>cluster_rmsd.txt_best$nmb
  grep -v "#" $iclu'_rmsd' | awk '{if (NF == 2) print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_rmsd.txt_best$nmb
end
#
# calculate average cluster rmsd from lowest energy structure
#
echo "#Cluster rmsd-Emin sd" >cluster_rmsd-Emin.txt_best$nmb

foreach iclu (file.nam_clust[0-9] file.nam_clust[0-9][0-9] file.nam_clust[0-9][0-9][0-9])
  echo $iclu |sed -e 's/list/nam/' |awk '{printf "%s ",$1}' >>cluster_rmsd-Emin.txt_best$nmb
  grep -v "#" $iclu'_rmsd-Emin' | awk '{if (NF == 2) print $2}' | head -n $nmb |$HADDOCKTOOLS/average.sh >>cluster_rmsd-Emin.txt_best$nmb
end  
#
# Combine results in one file
# 
paste cluster_haddock-score.txt_best$nmb cluster_rmsd.txt_best$nmb cluster_rmsd-Emin.txt_best$nmb cluster_ener.txt_best$nmb cluster_viol.txt_best$nmb cluster_bsa.txt_best$nmb cluster_dH.txt_best$nmb cluster_Edesolv.txt_best$nmb \
 | awk '{print $1,$2,$3,$5,$6,$8,$9,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$35,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$48,$49,$51,$52,$54,$55}' >clusters.stat_best$nmb
head -n 1 clusters.stat_best$nmb >clusters_haddock-sorted.stat_best$nmb
sort -gk2 clusters.stat_best$nmb |grep file >>clusters_haddock-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_ene-sorted.stat_best$nmb
sort -gk9 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_ene-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_nb-sorted.stat_best$nmb
sort -gk11 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_nb-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_nbw-sorted.stat_best$nmb
sort -gk13 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_nbw-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_air-sorted.stat_best$nmb
sort -gk19 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_air-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_sani-sorted.stat_best$nmb
sort -gk25 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_sani-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_vean-sorted.stat_best$nmb
sort -gk27 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_vean-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_dani-sorted.stat_best$nmb
sort -gk29 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_dani-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_bsa-sorted.stat_best$nmb
sort -grk43 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_bsa-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_dH-sorted.stat_best$nmb
sort -gk45 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_dH-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_Edesolv-sorted.stat_best$nmb
sort -gk47 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_Edesolv-sorted.stat_best$nmb

else
#
# Combine results in one file
# 
paste cluster_haddock-score.txt_best$nmb cluster_ener.txt_best$nmb cluster_viol.txt_best$nmb cluster_bsa.txt_best$nmb cluster_dH.txt_best$nmb cluster_Edesolv.txt_best$nmb \ 
| awk '{print $1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$29,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$42,$43,$45,$46,$48,$49,$51,$52}' >clusters.stat_best$nmb
head -n 1 clusters.stat_best$nmb >clusters_haddock-sorted.stat_best$nmb
sort -gk2 clusters.stat_best$nmb |grep file >>clusters_haddock-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_ene-sorted.stat_best$nmb
sort -gk5 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_ene-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_nb-sorted.stat_best$nmb
sort -gk7 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_nb-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_nbw-sorted.stat_best$nmb
sort -gk9 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_nbw-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_air-sorted.stat_best$nmb
sort -gk15 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_air-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_sani-sorted.stat_best$nmb
sort -gk21 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_sani-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_vean-sorted.stat_best$nmb
sort -gk23 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_vean-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_dani-sorted.stat_best$nmb
sort -gk25 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_dani-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_bsa-sorted.stat_best$nmb
sort -grk35 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_bsa-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_dH-sorted.stat_best$nmb
sort -gk37 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_dH-sorted.stat_best$nmb
head -n 1 clusters_haddock-sorted.stat_best$nmb >clusters_Edesolv-sorted.stat_best$nmb
sort -gk39 clusters_haddock-sorted.stat_best$nmb |grep file >>clusters_Edesolv-sorted.stat_best$nmb


endif

endif

if ($compresspdb == true) then
  gzip *.pdb
endif

touch CLUSTANA_DONE

goto exit

usage:

echo 'ana_cluster.csh: calculates average energies and RMSDs on a cluster basis'
echo ' '
echo ' usage: ana_cluster.csh [-best num] clustering_file_output'
echo ' '
echo ' The -best option allows to specify how many structures are being used to'
echo ' to calculate the average cluster values. If specified, additional files with'
echo ' extension .stat_best_num will be created.'


exit:
