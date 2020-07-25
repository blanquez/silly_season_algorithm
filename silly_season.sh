#!/bin/bash
declare -a seeds=('798245613' '123456789' '25022020' '17042026' '459268694')
for i in "${seeds[@]}"
do
	echo -e "\n$i:\n"
	./code/silly_season iris_set.dat iris_set_const_10.const 3 $i
	./code/silly_season newthyroid_set.dat newthyroid_set_const_10.const 3 $i
	./code/silly_season rand_set.dat rand_set_const_10.const 3 $i
	echo
	./code/silly_season iris_set.dat iris_set_const_20.const 3 $i
	./code/silly_season newthyroid_set.dat newthyroid_set_const_20.const 3 $i
	./code/silly_season rand_set.dat rand_set_const_20.const 3 $i     
done
