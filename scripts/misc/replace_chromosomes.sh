samtools view -H /genomes/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/Assembly/LP2000274-DNA_B11.bam | sed  -e 's/SN:chrM/SN:MT/' -e 's/SN:chr/SN:/'  | samtools reheader - /genomes/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/Assembly/LP2000274-DNA_B11.bam > /genomes/analysis/by_date/2016-11-18/RAREP40001/LP2000274-DNA_B11/coverage/LP2000274-DNA_B11.wochr.bam

