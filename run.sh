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
URR="urr.txt"
./urr.py "$BLAST_OUT" "$FILE" # > "$URR"

#
# Find E2BS
#
FIMO_OUT="fimo-out"
fimo "$URR" > "$FIMO_OUT"
