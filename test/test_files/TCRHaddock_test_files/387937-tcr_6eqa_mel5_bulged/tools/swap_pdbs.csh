#!/bin/csh
#
awk '{if ($3 < $2){print $1}}' file.nam_rmsd | grep pdb > file.swap

foreach i (`cat file.swap`)

  pdb_xsegchain $i | pdb_selchain -B | pdb_chain -A | grep ATOM > swap
  pdb_xsegchain $i | pdb_selchain -A | pdb_chain -B | grep ATOM >> swap
  pdb_xsegchain $i | pdb_selchain -C | grep ATOM >> swap
  grep REMARK $i >$i:r_swap.pdb
  pdb_reatom swap | pdb_chainxseg |pdb_chain >> $i:r_swap.pdb
  echo END >> $i:r_swap.pdb
  mv $i $i:r.pdb-store
  mv $i:r_swap.pdb $i
  \rm swap
end

