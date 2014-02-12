#!/bin/bash
if [ $# < 1 ] then;
    PEP_FASTA = pep.fasta
fi
$ROSETTA/rosetta_tools/perl_tools/getFastaFromCoords.pl -pdbfile ../native.pdb -chain B > $PEP_FASTA;
export RECEPTOR_LENGTH=`awk 'substr($0,13,3)==" CA" && substr($0,22,1)=="A"' ../native.pdb | wc -l`
export PEP_LENGTH=`awk 'substr($0,13,3)==" CA" && substr($0,22,1)=="B"' ../native.pdb | wc -l`
echo RECEPTOR_LENGTH $RECEPTOR_LENGTH
echo PEP_LENGTH $RECEPTOR_LENGTH
echo ROSETTA: ''$ROSETTA''
if [ $PEP_LENGTH -ge 9 ]; then
    export frag_sizes="3,5,9"
    export frag_sizes_list="3 5 9"
else
    export frag_sizes="3,5"
    export frag_sizes_list="3 5"
fi
$ROSETTA/rosetta_tools/fragment_tools/make_fragments.pl -frag_sizes $frag_sizes -n_frags 500 -id fragB $PEP_FASTA
    # shift frags by receptor length
for i in $frag_sizes_list; do
    awk '{if ( substr ( $0,1,3 ) == "pos" ) {print substr ( $0,0,18 ) sprintf ("%4d",substr ( $0,19,4 ) + '"$RECEPTOR_LENGTH"' ) substr ( $0,23,1000 ) ; } else {print ; }}' \
        fragB.500.${i}mers > frags.${i}mers.offset
done
cd ..
