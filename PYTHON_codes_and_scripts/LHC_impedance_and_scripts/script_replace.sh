#!/bin/bash

# change all files of the form *{first_argument}*
# by replacing all occurrences of {second_arguement} in them
# by {third argument}

files=$1
str1=$2
str2=$3

for i in `ls *${files}*`; do

sed "s/${str1}/${str2}/g" $i > ${i}_tmp
mv ${i}_tmp $i

done

