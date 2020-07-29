#!/bin/sh
# Default directories ---
ROOTDIR=$PWD
OUTPUTDIR="Rmarkdown"

# Command line options ---
usage() { echo "Usage: $0 [-r <path-to-rscript>] [-d <path-to-rootdir>] [-o <path-to-outdir>]" 1>&2; exit 1; }

while getopts ":hr:d:o:" opt; do
	case $opt in
		h)	usage;;
		r)	RSCRIPT=$OPTARG;;
		d)	ROOTDIR=$OPTARG;;
		o)	OUTPUTDIR=$OPTARG;;
		\?)	echo "Invalid option: ${OPTARG}" 1>&2; usage;;
		:)	echo "Invalid option: ${OPTARG} requires an argument" 1>&2; usage;;
	esac
done

shift $((OPTIND -1))

if [ -z "${RSCRIPT}" ]; then
    usage
fi
# Execution ---
echo ""
echo "[*] Render Rmarkdown report from Rscript [*]"
echo ""
echo " -- Root directory: ${ROOTDIR}"
echo " -- Output directory: ${OUTPUTDIR}"
echo ""

echo Rscript --vanilla -e "rmarkdown::render('${RSCRIPT}', output_dir='${OUTPUTDIR}')" ${ROOTDIR}

Rscript --vanilla -e "rmarkdown::render('${RSCRIPT}', output_dir='${OUTPUTDIR}')" ${ROOTDIR}
