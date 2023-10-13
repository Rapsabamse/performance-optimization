#!/bin/bash
# Author: David Holmqvist <daae19@student.bth.se>

echo "NOTE: this script relies on the binaries blur, blur_par, threshold and threshold_par to exist"

status=0
temp_dir=$(mktemp -d 2> /dev/null || mktemp -d -t 'mytmpdir')
red=$(tput setaf 1)
reset=$(tput sgr0)

for thread in 1 2 4 8 16 32
do
    for image in im1 im2 im3 im4
    do
        ./blur 15 "data/$image.ppm" "$temp_dir/blur_${image}_seq.ppm"
        ./blur_par 15 "data/$image.ppm" "$temp_dir/blur_${image}_par.ppm" $thread

        if ! cmp -s "$temp_dir/blur_${image}_seq.ppm" "$temp_dir/blur_${image}_par.ppm"
        then
            echo "${red}Incongruent output data detected when blurring image $image.ppm with $thread thread(s)${reset}"
            status=1
        fi

        rm "$temp_dir/blur_${image}_seq.ppm" "$temp_dir/blur_${image}_par.ppm"

        ./threshold "data/$image.ppm" "$temp_dir/threshold_${image}_seq.ppm"
        ./threshold_par "data/$image.ppm" "$temp_dir/threshold_${image}_par.ppm" $thread

        if ! cmp -s "$temp_dir/threshold_${image}_seq.ppm" "$temp_dir/threshold_${image}_par.ppm"
        then
            echo "${red}Incongruent output data detected when thresholding image $image.ppm with $thread thread(s)${reset}"
            status=1
        fi

        rm "$temp_dir/threshold_${image}_seq.ppm" "$temp_dir/threshold_${image}_par.ppm"
    done
done

exit $status