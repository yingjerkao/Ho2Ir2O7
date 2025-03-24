#!/bin/bash

for i in `ls input`;
do
	../code/sweep "$i".raw 0 < input/$i
done
