#!/bin/bash

for i in {0..29}
do
	for j in {0..29}
	do
		qsub 'jobs/TRF_'$i'_'$j
	done
done
