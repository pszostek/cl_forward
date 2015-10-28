#!/usr/bin/env bash

function tracks_in_reference() {
    reference_run=$1
    cat results_serial/$reference_run.out | grep Track | wc -l
}

function current_tracks() {
    ./clpixel -serial mcdata/$1.dat | grep Found | grep -Eo '[0-9]*'
}

set $(seq 0 50) 77

echo   "--------+------+------+"
printf "  Input | Act. | Ref. |\n"
echo   "--------+------+------+"
while [ "$#" -gt 0 ]; do
    tracks=$(current_tracks $1)
    reference_tracks=$(tracks_in_reference $1)
    printf "%3d.dat | %4d | %4d |\n" $1 $tracks $reference_tracks
    if [ $tracks -ne $reference_tracks ]; then
        echo   "--------+------+------+"
        echo -e "\e[31mTracks differ for the $1.dat dataset\e[39m"
        echo "Reference result: $reference_tracks"
        echo "Current result:   $tracks"
        exit 1
    fi
    shift
done

echo   "--------+------+------+"
echo -e "\e[32mCurrent and reference results match.\e[39m"
echo ""