#!/bin/bash

for i in {0..19}
do
	for j in {0..29}
	do
		qsub 'jobs/ring_'$i'_'$j
	done
done
