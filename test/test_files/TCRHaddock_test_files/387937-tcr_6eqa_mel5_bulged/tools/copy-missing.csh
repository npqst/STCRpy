#!/bin/csh -f
#
if ($#argv != 2) goto usage
set root=$1
set num=$2

set it0file='./structures/it0/'`head -$num structures/it0/file.nam |tail -1`
set it1file='./structures/it1/'$root'_'$num'.pdb'
set it0seed='./structures/it0/'`head -$num structures/it0/file.nam |tail -1 | sed -e 's/pdb/seed/g'`
set it1seed='./structures/it1/'$root'_'$num'.seed'
if (! -e $it1file) then
  cp $it0file $it1file
endif
goto exit

usage:
echo "Usage: copy-missing.csh rootfilename structure-number"
exit:

