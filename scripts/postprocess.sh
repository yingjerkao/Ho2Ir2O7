#!/bin/bash

for i in `ls input`;
do
	python3 pool.py "$i".raw "$i".pool
	python3 nicetable.py "$i".pool $i "$i".table
done
