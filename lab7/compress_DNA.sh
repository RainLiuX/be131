#!/bin/bash

cmp_file="nt_seq"
out_file="compress_01.output"
rm -f ${cmp_file}.* ${out_file}
echo "Time for gzip" >> ${out_file}
{ time gzip -k ${cmp_file}; } 2>> ${out_file}
echo "Time for bzip2" >> ${out_file}
{ time bzip2 -k ${cmp_file}; } 2>> ${out_file}
echo "Time for pbzip2" >> ${out_file}
{ time pbzip2 -ck ${cmp_file} > ${cmp_file}.pbz2; } 2>> ${out_file}
echo "Time for Arithmatic compress" >> ${out_file}
{ time ArithmeticCompress ${cmp_file} ${cmp_file}.art; } 2>> ${out_file}
echo "Zizes of each file:" >> ${out_file}
ls -lh | grep ${cmp_file} >> ${out_file}
