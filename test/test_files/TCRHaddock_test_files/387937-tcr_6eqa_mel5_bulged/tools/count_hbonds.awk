#
# usage: 
#    awk -f chb.awk hbonds.disp | sort -n -K7
#
#
# replacement for ana_hbond.csh & count_hbonds.awk
#
# bugs in those: 
#   - possibly infinite while loops, infinite memory
#   - gln HE2[12] translated to HH2#
#   - asn HD2[12] not translated 
#

BEGIN { stderr="/dev/stderr" }

$1=="atoms" && $NF == "apart" {
    
    # substitute atoms
#    gsub(/HZ[123]/,"HZ#")
#    gsub(/OD[12]/,"OD#")  
#    gsub(/OE[12]/,"OE#")  
#    gsub(/HT[123]/,"HT#")
#    gsub(/HH1[12]/,"HH1#")
#    gsub(/HH2[12]/,"HH2#")
#   #gsub(/HD2[12]/,"HH2#")  # ERROR ? 
#   gsub(/HD2[12]/,"HD2#")   # !
#    gsub(/HE2[12]/,"HE2#")   # !

    # all substitutions in one line # if your awk doesn't have this use 
    # the gsub()'s above
    $0=gensub(/(HZ|OD|OE|HT|HH[12]|H[ED]2)[123]/,"\\1#","g",$0)

#   for comparison
#   $0=gensub(/(HZ|OD|OE|HT|HH[12])[123]/,"\\1#","g",$0)

    # get DA: donor-acceptor pair
    split($0,istore,"      ")
    DA=istore[1]

    # get distance 
    distcol=NF-2
    dist = $distcol

    sumdist[DA] += dist
    sum2dist[DA] += dist*dist
    numdist[DA]++
}

END { 
    for (DA in sumdist) { 
        avedist=sumdist[DA]/numdist[DA]
        ave2dist=sum2dist[DA]/numdist[DA]
        sddist=(ave2dist - (avedist*avedist))**(0.5)
        printf "d(H-D)=%8.2f +- %5.2f #hb=%5d %s\n",avedist,sddist,numdist[DA],DA 
            #| "/usr/bin/sort -n -k7" #optional
    }
    # what about the DONE line? why is it necessary
}
