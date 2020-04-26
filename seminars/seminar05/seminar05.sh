#!/bin/bash

samtools sort hela.bam -o hela.sorted.bam
samtools index hela.sorted.bam

htseq-count -s yes -a 10 hela.sorted.bam gencode.gtf > hela_count.txt
