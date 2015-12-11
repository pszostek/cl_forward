#!/usr/bin/env bash

READER="python analysis/trackreader.py"

function tracks_in_reference() {
    echo $($READER results_serial/$1.out results_serial/$1.out --tracks)
}

function current_tracks() {
    ./clpixel -serial -tex mcdata/$1.dat >/dev/null
    echo $($READER results/$1_serial_txt.out results/$1_serial_txt.out --tracks)
}

function match() {
   $READER results/$1_serial_txt.out results_serial/$1.out --test-equal > /dev/null
}


set $(seq 0 50) 77

echo   "--------+------+------+-------+"
printf "  Input | Act. | Ref. | Match |\n"
echo   "--------+------+------+-------+"
while [ "$#" -gt 0 ]; do
    cur_tracks=$(current_tracks $1)
    reference_tracks=$(tracks_in_reference $1)
    match $1
    if [ $? -eq 0 ]; then
        match_str='\e[32mY\e[39m'
        equal=true
    else
        match_str='\e[31mN\e[39m'
        equal=false
    fi
    printf "%3d.dat | %4d | %4d |" $1 $cur_tracks $reference_tracks
    echo -e "   $match_str   |"
    if [ $cur_tracks -ne $reference_tracks ]; then
        echo   "--------+------+------+-------+"
        echo -e "\e[31mTracks differ for the $1.dat dataset\e[39m"
        echo "Reference result: $reference_tracks"
        echo "Current result:   $cur_tracks"
        exit 1
    elif [ $equal == false ]; then
        echo "\e[31mTracks differ.\e[39m"
        exit 1
    fi

    shift
done

echo   "--------+------+------+-------+"
echo -e "\e[32mCurrent and reference results match.\e[39m"
echo ""
exit 0
