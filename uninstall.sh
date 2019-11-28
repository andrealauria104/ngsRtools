#!/bin/bash
# '''
# ngsRtools: A Suite of R Packages for NGS-based Epigenomics Data Analysis.
# 
# '''
# uninstall suite
# 
# 0. Variables/Functions ---
WD=$PWD
PKGS=$(ls -d packages/*tools | sed -e 's/.*\///g')

# 1. Uninstall suite ---
for P in $PKGS
do
	R CMD REMOVE $P
done
