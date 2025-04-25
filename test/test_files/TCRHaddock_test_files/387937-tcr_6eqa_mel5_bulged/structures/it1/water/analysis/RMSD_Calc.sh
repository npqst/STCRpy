#!/bin/sh

# Profit Location
#
# Define the location of profit
#

if [ `printenv |grep PROFIT | wc -l` -eq 0 ] 
then
  found=`which profit |wc -l`
  if [ $found -eq 0 ]
  then
     export PROFIT="/trinity/home/enmr/software/bin/profit-rocky"
  else
     export PROFIT=`which profit`
  fi
fi

# Python location
PYTHONANA=/trinity/home/enmr/software/python/bin/python2.7
# Subset of Atoms on which to Calculate RMSD
ATOMS="CA,C,N,O,P,C1*,C2*,C3*,C4*,C5*,C6*,C7*,C8*,C9*"

if [ $# -lt 1 ]
then
  echo "./iRMSD_calc.sh file.nam_clustX [interface/ligand] [symmetry] [debug]"
  echo "file.nam_clustX - List of paths for PDB files. First structure is taken as reference"
  echo
  echo "If interface option is given, the script automatically calculates an interface"
  echo "using the contact script and creates specific ZONES to be used in the RMSD calculation."
  echo
  echo "If symmetry is on, chain permutation is turned on."
  echo
  echo "If debug on, keeps output files from profit for each structure"
  echo "Files should have both chain ID and seg ID"
  exit
fi

if [ ! -f $PROFIT ]
then
    echo "Profit is not where you said so: $PROFIT"
    exit
fi

if [ $2 == "interface" ]
then
    echo "Interface Mode ON"
    INTERFACE_MODE="ON"
fi

if [ $2 == "ligand" ]
then
    echo "Ligand Mode ON"
    LIGAND_MODE="ON"
fi

SYMMETRY="OFF"
if [ $3 == "symmetry" ]
then
    echo "Symmetry Mode ON"
    SYMMETRY="ON"
fi
    
if [ -f $1 ]
then
  # Grep non empty lines
  structures=`egrep -v "(^$|^#)" $1`
  nstructs=`egrep -v "(^$|^#)" $1 | wc -l | tr -d " "`
  echo "Loaded $nstructs structures"

  # Parse Structure List
  reference=`echo $structures | cut -d" " -f1`
  echo "Using $reference as reference"
  mobi=`echo $structures | cut -d" " -f1-`
else
  echo "File not found: $1"
  exit
fi

# Get chains for permutations
if [ -z $INTERFACE_MODE ]
then
    # Parse chains from reference (Not including HETATM I guess..)
    if [ -f $reference ]
    then
      chains=$(egrep "^ATOM" $reference | cut -c 22 | uniq | tr -d "\n")
      echo "Read ${#chains} chains from reference file: $chains"
    else
      echo "Reference structure not found: $reference"
      exit
    fi
else
    if [ -z $HADDOCKTOOLS ]
    then
        echo "HADDOCKTOOLS is not defined."
        echo "Exiting..."
        exit
    fi
    # Parse chains from contacts
    refe_contacts=${reference%%.pdb}_5A.contacts
    $HADDOCKTOOLS/contact $reference 10.0 > $refe_contacts
    chains=$(awk '{print $2"\n"$5}' $refe_contacts | sort -u | tr -d "\n")
    echo "Read ${#chains} chains from contacts file: $chains"
fi

# Check if receptor chains are the same between reference and models
diff=0
if [ ! -z $LIGAND_MODE ]
then
    chain=${chains:0:1}
    l_ref=$(egrep "^.{21}$chain" $reference | cut -c 23-26 | uniq | wc -l | tr -d " ")
    l_mobi=$(egrep "^.{21}$chain" `head -2 $1 | tail -1` | cut -c 23-26 | uniq | wc -l | tr -d " ")
    if [ $l_ref -ne $l_mobi ]
    then
        echo "Receptor chains have not the same lenght, zones for ligand RMSD will be automatically defined by profit alignment."
        diff=1
        grep ' CA ' $reference | awk '{print $5$6"_"$4}' | grep -v "-" | sort > tmp1
        grep ' CA ' `head -2 $1 | tail -1` | awk '{print $5$6"_"$4}' | grep -v "-" | sort > tmp2
        zone_align=`comm -1 -2 tmp1 tmp2 | sed -e 's/_/\ /' | awk '{print "zone "$1"-"$1}' | grep $chain`
#        \rm tmp1 tmp2
    fi
#    echo $l_ref
#    echo $l_mobi
fi

if [ ! -z $SYMMETRY ]
then
# Check if all chains have the same length
# Only on reference
SYMOFF="false"
echo -e "Symmetry option on. Checking chain lengths: \c"
prev=0
for (( i=0; i<${#chains}; i+=1 ))
do
    chain=${chains:$i:1}
    # Calculates length of chain:
    # egrep "^.{21}$chain" - greps character $chain at 22nd position
    # cut -c 23-26 - cuts residue number
    l=$(egrep "^.{21}$chain" $reference | cut -c 23-26 | uniq | wc -l | tr -d " ")
    if [ $prev -eq 0 ]
    then
        prev=$l
    else
        if [ $l -ne $prev ]
        then
            echo -e "WARNING!! Chains have different lengths! (${chains:$i:1}:$l != ${chains:`expr $i-1`:1}:$prev)"
            echo "You main run into problems later on..."
            SYMOFF="true"
        fi
        prev=$l
    fi
done
echo "OK!"
fi

# Parse specific interface residues and builds PROFIT zones per chain
if [ ! -z $INTERFACE_MODE ]
then
    declare -a INTERFACE
    declare -a RES_PER_CHAIN # To later on match common residues between chains when permutating
    for (( i=0; i<${#chains}; i+=1 )) # iterates chains
    do
        chain=${chains:$i:1}
        sideA=$(awk '{print $1, $2}' $refe_contacts | grep $chain | sort -nuk1) # Get residues from contact (first half)
        sideB=$(awk '{print $4, $5}' $refe_contacts | grep $chain | sort -nuk1) # second half
        residues="${sideA}\n${sideB}\n" # concatenate residues
        # Ugly one liner!
        # Makes a unique sorted list of residues and with sed builds the profit ZONE CHAIN_XRESI-CHAIN_XRESI syntax in one go.
        # CHAIN_X will be replace later when permutating
        # Limited to residue numbers up to 5 characters. Change the 5 in the sed to change that!
        RES_PER_CHAIN[$i]=$(echo -e "$residues" | sort -nuk1 | awk '{print $1}')
        INTERFACE[$i]=$(echo -e "${RES_PER_CHAIN[$i]}" | sed -e "s/^\([0-9]\{1,5\}\)/ZONE CHAIN_${i}\1-CHAIN_${i}\1:CHAIN_${i}\1-CHAIN_${i}\1/") 
    done
fi
# $INTERFACE is an array with one entry per zone definition.
# We just need to change chain IDs with sed according to the permutations
# and we should get the decent combinations. Fingers crossed!

# Echoes the zones
# for (( i=0; i<${#INTERFACE[@]}; i+=1 ))
# do
#     echo "${INTERFACE[$i]}"
# done

echo "Generating ZONE definitions.."

if [ $SYMOFF == "true" ]
then
    echo "Symmetry option was turned on, but chains are not identical: Turning symmetry off"
fi

if [ ! -z $SYMMETRY && $SYMOFF != "true" ]
then
    echo "Symmetry option on. Generation all permutations"
    # Generate Permutations of chains, keeping the first chain fixed
    # E.G. For ABC gives: ABC, ACB
    firstchain=${chains:0:1}
    if [ -z $PYTHONANA ]
    then
        echo "Python is required to generate permutations of the chains"
        echo "Set the PYTHONANA variable at the beginning of the script to your Python installation."
        exit
    fi
    permutations=`$PYTHONANA -c "from itertools import permutations; print '\n'.join([''.join(i) for i in permutations('$chains') if i[0] == '$firstchain'])"`
else
    permutations=$chains
fi

# Define Zones based on permutations
# Store in array
# Zone clear and FIT already written

declare -a ZONES
j=0 # Counter for zones
for p in $permutations # Iterates over permutations
do  
    zone="ZONE CLEAR \n"
    if [ -z $INTERFACE_MODE ] && [ -z $LIGAND_MODE ]
    then # Uses all residues!
        for (( i=0; i<${#chains}; i+=1 )) # Iterates over chains
        do
            pair="ZONE ${chains:$i:1}*:${p:$i:1}*"
            zone="${zone}$pair\n"
        done
    else # Using $INTERFACE
        # Repetition of code to avoid if-in-for-loop performance loss
        if [ ! -z $INTERFACE_MODE ]
        then
            for (( i=0; i<${#chains}; i+=1 )) # Iterates over chains
            do
                chain_zone=$(echo "${INTERFACE[$i]}" | sed -e "s/CHAIN_${i}/${chains:$i:1}/" -e "s/CHAIN_${i}/${chains:$i:1}/" -e "s/CHAIN_${i}/${p:$i:1}/" -e "s/CHAIN_${i}/${p:$i:1}/")
                # Find common residues (used only in symmetrical complexes)
                if [ ! -z $SYMMETRY ]
                then
                    index=$(expr index $p ${chains:$i:1})
                    num=$(expr $index - 1 )
                    res_A=$(echo ${RES_PER_CHAIN[$i]} | tr "\n" " ")
                    res_B=$(echo ${RES_PER_CHAIN[$num]} | tr "\n" " ")
                    difference=`$PYTHONANA -c 'print " ".join(set("'"$res_A"'".split()).symmetric_difference("'"$res_B"'".split()))'`
                    for dif in $difference # removes differences from the zone file
                    do
                        chain_zone=$(echo -e "$chain_zone" | sed -e "/$dif/d")
                    done
                fi
                # Write ZONE
                zone="${zone}$chain_zone\n"
            done
        else # ZONE definition for ligand calculation
            if [ $diff -eq 0 ]
            then
                zone=${zone}" ZONE ${chains:0:1}* \n FIT \n"
            else
                zone=${zone}" "$zone_align" \n FIT \n"
            fi
            for (( i=1; i<${#chains}; i+=1 )) # Iterates over chains
            do
                pair="RZONE ${chains:$i:1}"
                pair="${pair}*\n"
                zone="${zone}$pair"
            done
        fi
    fi
    if [ -z $LIGAND_MODE ]
    then
        ZONES[$j]="$zone FIT" # Adds zone to array
    else
        ZONES[$j]="${zone}"
    fi
    j=`expr $j + 1` # increment counter
done
echo "Defined ${#ZONES[@]} zone(s)"
echo $ZONES
#echo -e "${ZONES[@]}" # Echoes all zones

# Prepare output file
cat /dev/null > rmsd.temp
if [ -z $LIGAND_MODE ]
then
    echo "#Interface RMSD" >> rmsd.temp
else
    echo "#Ligand RMSD" >> rmsd.temp
fi

# Run PROFIT iteratively on each structure
for s in $mobi
do
    echo -e "$s \c" >> rmsd.temp
#    echo ${ZONE[@]}
    IN="refe $reference\nmobi $s\nATOMS $ATOMS\n"${ZONES[@]}"\nquit"
    OUT=$(echo -e "$IN" | $PROFIT) # Capture output
    if [ -z $LIGAND_MODE ]
    then
        RMS=$(echo "$OUT" | grep RMS)
    else
        RMS=$(echo "$OUT" | grep RMS | tail -1)
    fi
    echo $RMS
    echo ${#ZONES[@]}
    # Error checking (number of RMS outputs equals number of FIT commands)
    if [ ${#chains} -gt "2" ]
    then
	echo "$RMS" | sort -nk 2 | awk '{print $NF}' | head -1 >> rmsd.temp
        echo "$OUT" > ${s}_profit.out
    else
	if [ $(echo "$RMS" | wc -l) -eq "${#ZONES[@]}" ]
	then
	    echo "$RMS" | sort -nk 2 | awk '{print $NF}' | head -1 >> rmsd.temp
	    echo "$OUT" > ${s}_profit.out
	else
	    echo "ERROR" >> rmsd.temp
	    echo "$OUT" > ${s}_profit.err
	fi
    fi
done

# Produce sorted files
if [ -z $LIGAND_MODE ]
then
    cat rmsd.temp > i-RMSD.dat
#    sort -k2 -t"_" -n rmsd.temp > i-RMSD.dat
    sort -k2 -n rmsd.temp > i-RMSD-sorted.dat
	else
    if [ ! -z $LIGAND_MODE ]
    then
        cat rmsd.temp > l-RMSD.dat
#        sort -k2 -t"_" -n rmsd.temp > l-RMSD.dat
        sort -k2 -n rmsd.temp > l-RMSD-sorted.dat
    else
        cat rmsd.temp > RMSD.dat
#        sort -k2 -t"_" -n rmsd.temp > RMSD.dat
        sort -k2 -n rmsd.temp > RMSD-sorted.dat
    fi
fi

# Clean up
#\rm rmsd.temp
echo $4
if [ -z $4 ]
then
    \rm *_profit.out *_profit.err &> /dev/null
fi
exit
