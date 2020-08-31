#BSUB parameters go here, if HPC system is used

#Setting paths and primer
input_data_loc="/workdir/kja68/raw_sequence_files_R1"
perl_code="/workdir/kja68/Trimming_HS_Primers.pl"
primer="GTCGGTAAAACTCGTGCCAGC"

cd $input_data_loc

mkdir Trimmed

date

for f in *.gz ; do gunzip -c "$f" > Trimmed/"${f%.*}"; echo "${f}" | cut -f1 -d'.'; perl $perl_code  < Trimmed/"${f%.*}" $(echo "${f}" | cut -f 1 -d '.') $primer ; done

for f in *TrimmedHS_R1_FASTQ.fastq ; do echo $f.gz; gzip "$f"; mv $f.gz Trimmed/$f.gz; done

for f in *log.txt ; do mv $f Trimmed/$f; done

for f in $input_data_loc/Trimmed/*.fastq ; do rm $f ; done

echo "Trimming completed at: "
date