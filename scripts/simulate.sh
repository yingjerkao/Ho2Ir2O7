#!/bin/bash

for i in `ls input`;
do
	../code/sweep.out "$i".raw 0 < input/$i
done
