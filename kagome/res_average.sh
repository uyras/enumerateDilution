#!/bin/bash
cd res
rm -f g_aver_*.dat
holes=( 0 3 7 10 15 )
for i in ${holes[@]}
do 
	for f in g_${i}_*.dat
	do
		cat $f | awk -f ../res_aver.awk >> g_aver_${i}.dat
	done
done
