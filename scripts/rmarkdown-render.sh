#!/bin/sh
# Default directories ---
OUTPUTDIR="Rmarkdown"
RARGS=""

# Command line options ---
usage() { echo "Usage: $0 [-r <path-to-rscript>] [-o <path-to-outdir>] [-a <args-to-rscript>]" 1>&2; exit 1; }

while getopts ":hr:o:a:" opt; do
	case $opt in
		h)	usage;;
		r)	RSCRIPT=$OPTARG;;
		o)	OUTPUTDIR=$OPTARG;;
		a)	RARGS=$OPTARG;;
		\?)	echo "Invalid option: ${OPTARG}" 1>&2; usage;;
		:)	echo "Invalid option: ${OPTARG} requires an argument" 1>&2; usage;;
	esac
done

shift $((OPTIND -1))

if [ -z "${RSCRIPT}" ]; then
    usage
fi

if [[ ${RARGS} == "-h" ]]; then
    if [ -x ${RSCRIPT} ]; then
      Rscript --vanilla ${RSCRIPT} ${RARGS}
    else
      echo "\n[!] ${RSCRIPT} does not have help message for argument parsing.\n"
    fi
    exit 1
fi
# Execution ---
echo ""
echo "[*] Render Rmarkdown report from Rscript [*]"
echo ""
echo " -- Output directory: ${OUTPUTDIR}"
echo ""

echo Rscript --vanilla -e "rmarkdown::render('${RSCRIPT}', intermediates_dir='${PWD}', output_dir='${OUTPUTDIR}')" ${RARGS}

Rscript --vanilla -e "rmarkdown::render('${RSCRIPT}', intermediates_dir='${PWD}', output_dir='${OUTPUTDIR}')" ${RARGS}
