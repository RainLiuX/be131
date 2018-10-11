#!/bin/bash

cmp_file=$1
out_file=comp_${cmp_file}.output
rm -f ${cmp_file}.* ${out_file}
echo "Time for gzip" >> ${out_file}
{ time -p gzip -k ${cmp_file}; } 2>> ${out_file}
echo "Time for bzip2" >> ${out_file}
{ time -p bzip2 -k ${cmp_file}; } 2>> ${out_file}
echo "Time for pbzip2" >> ${out_file}
{ time -p pbzip2 -ck ${cmp_file} > ${cmp_file}.pbz2; } 2>> ${out_file}
echo "Time for Arithmatic compress" >> ${out_file}
{ time -p ArithmeticCompress ${cmp_file} ${cmp_file}.art; } 2>> ${out_file}
echo "Sizes of each file:" >> ${out_file}
ls -l ${cmp_file}* >> ${out_file}
