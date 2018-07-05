#!/bin/bash


dir="res"
#holes=$(seq 1 18) # варианты возможных разбавлений образца, список
holes=( 0 3 7 10 15 )
h=$(LANG=en_US seq 0.0 0.01 5.0)     # варианты возможных полей, список

echo "$i $j"

for i in ${holes[@]}
do
	for j in $h
	do
		cd $dir
		printf "$i\n$j" | .././main.o
		cd ..
		#./calc_heating.sh $i $j $dir
		#FILES="${dir}/g_${i}_*_${j}.dat"
		#rm -f $FILES
	done
done