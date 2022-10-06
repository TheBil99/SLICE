#!/bin/bash

if [ $# -eq 1 ] ; then
	if [ ! -d ./data/"$1" ]; then
	   	mkdir -p ./data/"$1";
	  	mkdir -p ./data/"$1"/beta;
		mkdir -p ./data/"$1"/NPMI;
		mkdir -p ./data/"$1"/PI2;
		mkdir -p ./data/"$1"/PI2_inter_beta_evaluation;
		mkdir -p ./data/"$1"/PI3;
		mkdir -p ./data/"$1"/PI3_viewpoint;
	fi
else 
	echo 'pass only one argument containing the name of the dataset '
fi

