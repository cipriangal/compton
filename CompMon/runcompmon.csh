#!/bin/csh

if ($#argv < 1) then
    echo "Usage: $0 [runnumber] [optional - number of events - default 1M]"
else
    if ($#argv == 2) set nEvent = $2
    else set nEvent = 1000000
    set runnum = $1
endif
	    
mkdir -p rootfiles
set output = "./rootfiles/compmon_"$runnum".root"
rm lnkCompMon.output
ln -s $output lnkCompMon.output
rm $output

rm lnkCompMon.input
set input = "/home/compton/data/raw/run"$runnum".dat"
ln -s $input lnkCompMon.input

./bin/compmon $runnum $nEvent


