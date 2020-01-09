#!/bin/bash
# '''
# ngsRtools: A Suite of R Packages for NGS-based Epigenomics Data Analysis.
# 
# '''
# Build and install suite
# 
# 0. Variables/Functions ---
WD=$PWD
PKGS=$(ls packages | grep -v utilsRtools)

makeRtools () {
	
	cd $1
	BNAMES=${1##*/}
	make
	
	if [ $? -ne 0 ]; then
		echo "  [!] Make error."
		echo "  [!] $BNAMES not build/installed."
		rm -rf $BNAMES*.gz $BNAMES.Rcheck
		exit 1
	fi
}


echo ""
echo "=== ngsRtools installation script ==="
echo ""
# 1. Install dependencies ---
echo "[+] Checking and installing dependencies ..."
echo ""
${WD}/scripts/install-dependencies.r
echo "[+] All dependencies are correctly solved."
echo ""
# 2. Build utils-dependency ---
makeRtools ${WD}/packages/utilsRtools

# 3. Build others ---
for P in $PKGS
do
	makeRtools ${WD}/packages/${P}
done
