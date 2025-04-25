#!/bin/csh
#
if (! -e file.nam_all) then
  cp file.nam file.nam_all
endif
if (! -e file.cns_all) then
  cp file.cns file.cns_all
endif
if (! -e file.list_all) then
  cp file.list file.list_all
endif
if (! -e analysis_all) then
  mv analysis analysis_all
endif
#mv -f [A-Z]* analysis >&/dev/null

rm -f analysis file.list file.nam file.cns >&/dev/null
mkdir analysis_$1
ln -s file.list_$1 file.list
ln -s file.nam_$1 file.nam
ln -s file.cns_$1 file.cns
ln -s analysis_$1 analysis
