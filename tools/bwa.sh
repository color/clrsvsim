#!/bin/bash

set -e

# Usage: ./bwa.sh <bwa executable location> <reference genome fasta> <forward fastq> <reverse fastq> <bam destination (without the .bam extension)>

bamfile=`mktemp --tmpdir=$(pwd) -t XXXXXX.bam.tmp`
$1 mem -R '@RG\tID:foo\tSM:bar' -M -c 16000 -t 25 $2 $3 $4 | samtools view -uSb /dev/stdin > $bamfile
samtools sort -m 4G $bamfile $5
