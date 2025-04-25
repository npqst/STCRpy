#!/bin/csh
#
if ($#argv != 3) goto usage

cat <<_Eod_ >ene_rmsd.xmgr
# ACE/gr parameter file
#
@ TYPE xy
@ subtitle "HADDOCK run"
@ xaxis  label "Rmsd from lowest energy structure [A]"
@ yaxis  label "HADDOCK score [a.u.]"
@ s0 symbol 2
@ s0 symbol size 0.8
@ s0 symbol fill 0
@ s0 symbol color 1
@ s0 symbol linewidth 2
@ s0 linestyle 0
_Eod_

echo "grep -v '#' "$3" | awk '{print "\$$1","\$$2"}'>>ene_rmsd.xmgr" > cmd.tmp
chmod +x cmd.tmp
source ./cmd.tmp 
\rm cmd.tmp
cat <<_Eod2_ >>&ene_rmsd.xmgr
&
_Eod2_

goto exit

usage:

echo 'make_ene-rmsd_graph.csh:  generates a XMGR graph with data'
echo '                          taken from the user-specified columns'
echo ' '
echo '  usage: make_ene-rmsd_graph.csh  column_number1 column_number2 data_file'

exit:

