#!/bin/bash
declare -a seeds=('798245613' '123456789' '25022020' '17042026' '459268694')
for i in "${seeds[@]}"
do
	echo -e "\n$i:\n"
	./silly_season iris_set.dat iris_set_const_10.const 3 $i
	./silly_season newthyroid_set.dat newthyroid_set_const_10.const 3 $i
	./silly_season ecoli_set.dat ecoli_set_const_10.const 8 $i
	./silly_season rand_set.dat rand_set_const_10.const 3 $i
	echo
	./silly_season iris_set.dat iris_set_const_20.const 3 $i
	./silly_season newthyroid_set.dat newthyroid_set_const_10.const 3 $i
	./silly_season ecoli_set.dat ecoli_set_const_20.const 8 $i
	./silly_season rand_set.dat rand_set_const_20.const 3 $i     
done