#!/usr/bin/env bash

READER="python EventAnalyzer/trackreader.py"

REFERENCE_SOURCE=results_serial
EXEC_MODE=-serial
for i in "$@"; do
    case $i in
        -g|--gcc)
            REFERENCE_SOURCE=results_serial
            shift
            ;;
        -i|--icc)
            REFERENCE_SOURCE=results_serial_icc
            shift
            ;;
        -o|--openmp)
            EXEC_MODE=-openmp
            shift
            ;;
        -t|--tbb)
            EXEC_MODE=-tbb
            shift
            ;;
        *)
            echo "Unknown option provided: " $i
            exit 1
            ;;
    esac
done

function tracks_in_reference() {
    echo $($READER $REFERENCE_SOURCE/$1_serial_txt.out $REFERENCE_SOURCE/$1_serial_txt.out --tracks)
}

function current_tracks() {
    bin/x86_64/Release/clpixel $EXEC_MODE -tex mcdata/$1.dat >/dev/null
    echo $($READER results/$1_serial_txt.out results/$1_serial_txt.out --tracks)
}

function match() {
   $READER results/$1_serial_txt.out $REFERENCE_SOURCE/$1_serial_txt.out --test-equal > /dev/null
}

# get rid of all the results present locally
rm results/*

set $(seq 0 50) 77
equal="true"

echo   "--------+------+------+-------+"
printf "  Input | Act. | Ref. | Match |\n"
echo   "--------+------+------+-------+"
while [ "$#" -gt 0 ]; do
    cur_tracks=$(current_tracks $1)
    reference_tracks=$(tracks_in_reference $1)
    match $1
    if [ $? -eq 0 ]; then
        match_str='\e[32mY\e[39m'
    else
        match_str='\e[31mN\e[39m'
        equal="false"
    fi
    printf "%3d.dat | %4d | %4d |" $1 $cur_tracks $reference_tracks
    echo -e "   $match_str   |"
    shift
done

echo   "--------+------+------+-------+"
if [ $equal == "true" ]; then
    echo -e "\e[32mCurrent and reference results match.\e[39m"
else
    echo -e "\e[31mCurrent and reference results don\'t match.\e[39m"
fi

echo ""
exit 0
