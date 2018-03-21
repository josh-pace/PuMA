#!/bin/bash

set -u

FILE=""

function USAGE() {
    printf "Usage:\n  %s -f FILE\n\n" "$(basename "$0")"

    echo "Required arguments:"
    echo " -f FILE"
    echo
    exit "${1:-0}"
}

[[ $# -eq 0 ]] && USAGE 1

while getopts :f:h OPT; do
    case $OPT in
        f)
            FILE="$OPTARG"
            ;;
        h)
            USAGE
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument."
            exit 1
            ;;
        \?)
            echo "Error: Invalid option: -${OPTARG:-""}"
            exit 1
    esac
done

[[ -z "$FILE" ]] && USAGE
    
if [[ ! -e "$FILE" ]]; then
    echo "-f \"$FILE\" is not a file."
    exit 1
fi

ORFS="$(basename "$FILE").orfs"
./orf_predictor.py -o "$ORFS" "$FILE"

BLAST_OUT="blast-out"
blastp -query "$ORFS" -subject blast_database/conserved.fa -outfmt 6 -evalue 0.001 > "$BLAST_OUT"

cat "$BLAST_OUT"

# 
# Find URR
# 
URR="urr.fa"
./urr.py "$BLAST_OUT" "$FILE" > "$URR"

#
# Find E2BS
#
#FIMO_OUT_E2BS="fimo-out-e2bs"
fimo --oc "$FIMO_OUT_E2BS" --norc --verbosity 1 --thresh 1.0E-3 meme_3000_TOTAL.txt "$URR" 

#FIMO_RESULT_E2BS=""$FIMO_OUT_E2BS"/fimo.txt"
#E2BS_FINAL="e2bs.txt"
#./e2bs.py "$FILE" "$URR" "$FIMO_RESULT_E2BS" > "$E2BS_FINAL" 

#
#Find E1BS
#
#FIMO_OUT_E1BS="fimo-out-e1bs"
#fimo --oc "$FIMO_OUT_E1BS" --norc --verbosity 1 --thresh 1.0E-4 --bgfile background_model_E1BS.txt meme_E1BS_1motif_18_21.txt "$URR"

#FIMO_RESULTS_E1BS=""$FIMO_OUT_E1BS"/fimo.txt"
#E1BS_FINAL="e1bs.txt"
#./e1bs.py "$FILE" "$URR" "$FIMO_RESULTS_E1BS" > "$E1BS_FINAL"
