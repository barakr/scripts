#!/bin/csh

# assumes gnuin is present in running dir

# usage: cluster.sh <top x structures> <radius> <score.sc> <native.pdb> <decoys.silent> <stype=column number of score according to which to select the top x structures>
# e.g.
# $ cluster_shuf5000.sh 500 2 ../score_flags_abinitio.1.sc ../run_src/native.pdb ../decoys_flags_abinitio.silent reweighted_sc
#    ---> clusters top 500 by reweighted score with 2A peptide bb rms radius.

set PATH_TO_EXE=$ROSETTABIN
set PATH_TO_DB=$ROSETTADB

set top=$1
set radius=$2
set scorefile=$3
set native=$4
set decoys=$5
set stype=$6

#calculate column for specified score type, and the column for the model description (name)
set stype_column=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="'"$stype"'") print i}}'`
set description_column=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="description") print i}}'`
set rmsBB_if_column=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="rmsBB_if") print i}}'`
set bestRMS_5mer_all_column=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="bestRMS_5mer_all") print i}}'`
set bestRMS_6mer_all_column=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="bestRMS_6mer_all") print i}}'`
set startRMSbb_column=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="startRMSbb") print i}}'`
set startRMSphipsi_column=`head -1 $scorefile | awk '{for (i=1;i<=NF;i++){if ($i=="startRMSphipsi") print i}}'`
echo $stype_column, $description_column, $rmsBB_if_column, $bestRMS_5mer_all_column, $bestRMS_6mer_all_column, $startRMSbb_column, $startRMSphipsi_column

#calculate actual clustering radius taking intoaccount protein lenght.
set len=`grep CA $native | wc -l`
set plen=`awk '$5=="B"' $native | grep CA | wc -l`
set actualR=`date | awk '{print sqrt('$plen'/'$len')*'$radius'}'`

#perform clustering run
if (! -e c.0.0.pdb ) then
	$PATH_TO_EXE/cluster.linuxgccrelease -in:file:silent $decoys -in:file:silent_struct_type binary -database $PATH_TO_DB -cluster:radius $actualR -in:file:fullatom -tags `shuf $scorefile | head -n 5000 | sort -nrk ${stype_column} | tail -$top | awk '{print $'"${description_column}"'}'` -silent_read_through_errors >! clog
endif

set topW5=0
@ topW5 = $top + 5

#report top clusters ranked by score $stype
tail -$topW5 clog | head -$top | awk '{print $(NF-2), $(NF-1) }' | sort -nk2 | sort >! tmp.listA
cat $scorefile | awk '{print $'"$description_column"',$'"$stype_column"',$'"$rmsBB_if_column"',$'"$bestRMS_5mer_all_column"',$'"$startRMSbb_column"',$'"$startRMSphipsi_column"'}' | sort > tmp.listB
set output_file=clusters_by_${stype}.txt
echo description cluster $stype > ${output_file}
# representative select - column by which to select representative from cluster, 3 by energy, 4 by closest to native
set representative_select_col=4
# cluster_sort_col - column by which to sort clusters, 2 by cluster size, 3 by energy
set cluster_sort_col = 3
join -1 1 -2 1 tmp.listA tmp.listB | sort -nk2 -k $representative_select_col | awk 'BEGIN{cur=0}{if(cur==$2){print; cur++;}}' | sort -nk $cluster_sort_col > ${output_file}

#rm -rf tmp.list?
