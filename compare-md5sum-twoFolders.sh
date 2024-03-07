#!/bin/zsh

# Aim: compare the md5sum of all files under two input directories
# NOTE: This script assumes the directory structure is the same in both directories

# Check if the number of arguments is correct
if [ $# -ne 2 ]; then
    echo "Usage: $0 <directory1> <directory2>"
    exit 1
fi

dir1=$1
dir2=$2

# Find all files in directory1, strip the prefix, and sort them
files1=($(find "$dir1" -type f | sed "s|^$dir1/||" | sort))
# Find all files in directory2, strip the prefix, and sort them
files2=($(find "$dir2" -type f | sed "s|^$dir2/||" | sort))

# count the number of files in directory1 and directory2
count1=${#files1[@]}
count2=${#files2[@]}
echo "Number of files in $1: $count1"
echo "Number of files in $2: $count2"

# Initialize the counter
counter=0
mismatch=0

for rel_path in $files1; do
    file1="$dir1/$rel_path"
    file2="$dir2/$rel_path"

    if [ ! -f "$file2" ]; then
        echo "File $rel_path does not exist in $dir2"
        ((mismatch++))
    else
        md5sum1=$(md5sum "$file1" | cut -d' ' -f1)
        md5sum2=$(md5sum "$file2" | cut -d' ' -f1)
        if [ "$md5sum1" != "$md5sum2" ]; then
            echo "------------------------------"
            echo "Mismatch found:"
            echo "$file1"
            echo "$file2"
            echo "MD5 checksums do not match"
            echo "------------------------------"
            ((mismatch++))
        fi
    fi
    ((counter++))
    echo -ne "Compared files: $counter\r"
done

echo
echo "Total files compared: $counter"
echo "Total mismatches found: $mismatch"
echo
