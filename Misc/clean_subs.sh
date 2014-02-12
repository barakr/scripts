#!/bin/bash
for f in $@; do
    echo Cleaning $f
    cd $f;
    if ls movie*.rmf &> /dev/null; then
        for i in movie*.rmf; do
            j=S${i%\.rmf}.rmf;
            echo $i $j;
            ~/RMF_bin/rmf_slice -s 20 $i  $j >& LOG.slice_next ;
            if [ -e $j ]; then rm $i; fi;
        done;
    fi
    cd ..;
 done &