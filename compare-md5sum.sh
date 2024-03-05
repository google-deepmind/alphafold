#!/bin/zsh

# Aim: compare the md5sum of all files under two input directories

# Check if the number of arguments is correct
if [ $# -ne 2 ]; then
    echo "Usage: $0 <directory1> <directory2>"
    exit 1
fi

# Get the list of files in directory1
files1=($(find $1 -type f))

# Get the list of files in directory2
files2=($(find $2 -type f))

# count the number of files in directory1 and directory2
count1=${#files1[@]}
count2=${#files2[@]}
echo "Number of files in $1: $count1"
echo "Number of files in $2: $count2"

# Initialize the counter
counter=0

# Compare the MD5 checksums of the files
for file1 in $files1; do
    file2=$(find $2 -type f -name "$(basename $file1)" -print -quit)
    if [ -z "$file2" ]; then
        echo "File $(basename $file1) does not exist in directory2"
    else
        md5sum1=$(md5sum $file1 | cut -d' ' -f1)
        md5sum2=$(md5sum $file2 | cut -d' ' -f1)
        if [ "$md5sum1" != "$md5sum2" ]; then
            echo "------------------------------"
            echo $file1
            echo $file2
            echo "MD5 checksums of $(basename $file1) and $(basename $file2) do not match"
            echo "------------------------------"
        fi
    fi
    # Increment the counter
    ((counter++))
    # Echo the number of compared files, overwriting the line each time
    echo -ne "Compared files: $counter\r"
done

# Print a newline at the end to move the terminal cursor to the next line
echo