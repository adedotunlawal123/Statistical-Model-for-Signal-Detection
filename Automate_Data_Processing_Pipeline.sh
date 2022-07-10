#!/bin/bash
#list="0 1 2"
mkdir {0..2}hr
i=0
while [ $i -lt 3 ]
do
#for i in $list; do
    #if $i<3; then
	    cd ${i}hr
	    cp ~/project/umfs_son_2020/profiles/umfs.2020{09..11}{01..31}0${i}.dat ~/project/umfs_son_2020/svds/${i}hr
	    i=`expr $i + 1`
            cp ~/pp/testing/data_analysis data_analysis
            ./data_analysis
            chmod +x dd.py
            ./dd.py
	    cp ~/pp/testing/eofs_pcs_comparecdfs.py eofs_pcs_comparecdfs.py
	    chmod +x eofs_pcs_comparecdfs.py
	    ./eofs_pcs_comparecdfs.py
	    rm -f umfs*
	    cd ..
    #fi
done


    #when copying this script from a terminal editor make sure to copy from this terminal editor to another terminal editor in linux
    #When copying this script to new directory change the profiles, testing ( or whatever required directory name ) to the new directories name, absolute or relative.
    #Also remember to change/include name of scripts for cases where new scripts are included
    #Also change i, 2 in the mkdir and 3 in the while loop to the required value ( +1 above the directory with the highest number in that of while loop ) for the analysis.
    #For the 10-23hr profiles, you change i to 10 and take away the 0 (zero) before ${i} in umfs.2106{01..30}0${i}.dat
    
    
