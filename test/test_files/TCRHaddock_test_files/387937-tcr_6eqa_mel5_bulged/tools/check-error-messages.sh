#!/usr/bin/env bash
#
export found=`tail -n 500 $1 | grep "CHAIN LENGHT FOR SYMMETRY RESTRAINTS DO NOT MATCH" | grep -v display | wc -l | awk '{print $1}'`

if [ "$found" -ne "0" ]; then
  tail -n 500 $1 |grep -A 1 "CHAIN LENGHT FOR SYMMETRY RESTRAINTS DO NOT MATCH" >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "NCS-restraints error encountered: Improperly defined non-crystallographic symmetry"`

if [ "$found" -ne "0" ]; then
  tail -n 500 $1 |grep "NCS-restraints error encountered: Improperly defined non-crystallographic symmetry" >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "error in SYMMETRY potential, check NOE table"`

if [ "$found" -ne "0"  ]; then
  tail -n 500 $1 |grep "error in SYMMETRY potential, check your symmetry restraints definition" >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "exceeded allocation for NOE-restraints"`

if [ "$found" -ne "0" ]; then
  tail -n 500 $1 |grep "exceeded allocation for NOE-restraints" >>FAILED
  echo "Check your definition of active and passive residues" >>FAILED
  echo "Make sure to filter those for solvent accessibility"  >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "SELRPN error encountered: parsing error"`

if [ "$found" -ne "0" ]; then
  tail -n 500 $1 |grep -B 10 "SELRPN error encountered: parsing error" >>FAILED
  echo "Check your restraint files " >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "PARSER error encountered: Encountered too many parsing errors"`

if [ "$found" -ne "0" ]; then
  tail -n 500 $1 |grep -B 10 "PARSER error encountered: Encountered too many parsing errors" >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "XMREAD error encountered:  sectioning of map incompatible with resolution"`

if [ "$found" -ne "0" ]; then
  tail -n 500 $1 | grep -A 2 -B 2 "XMREAD error encountered:  sectioning of map incompatible with resolution" >>FAILED
  echo "Check your EM map resolution and sectioning" >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "ALLHP error encountered: not enough memory available"`

if [ "$found" -ne "0" ]; then
  head -n 100 $1 | grep "Running on" >>WARNING
  tail -n 500 $1 | grep "ALLHP error encountered: not enough memory available" >>WARNING
  echo "Check your definition of active and passive residues"     >>WARNING
  echo "Make sure to filter those for solvent accessibility"      >>WARNING
  echo "Try to decrease the size of your system where possible"   >>WARNING
fi

export found=`tail -n 500 $1 | grep -c "error encountered: missing SCATter definition for SELEcted atoms"`

if [ "$found" -ne "0" ]; then
  tail -n 500 $1 |grep -B 5 "error encountered: missing SCATter definition for SELEcted atoms" >>FAILED
  echo "Unsupported atoms/molecules for cryo-EM restraints" >>FAILED
fi

export found=`tail -n 500 $1 | grep -c "ROTMAT error encountered: rotation vector has zero length"`

if [ "$found" -ne "0" ]; then
  head -n 100 $1 | grep "Running on" >>WARNING
  tail -n 500 $1 | grep "CNS problem encounter"                         >>WARNING
  echo "Check your input parameters and restraints"                     >>WARNING
  echo "Possibly try turning off the sampling of 180 degrees rotattion" >>WARNING
fi

if [ -e WARNING ]; then
  export warnlength=`wc -l WARNING | awk '{print $1}'`
  # if more than 100 structures fail we will flag the run as failed.
  if [ "$warnlength" -gt "500" ]; then
    \cat WARNING >FAILED
  fi
fi
