#!/usr/bin/env bash
# This scripts will apply the ASV pipeline used in Littlefair et al.
# 2019.
# ---------
# Arguments
# ---------_
# 1) file with the list fileset prefix (before the _R1)
# 2) Forward primer
# 3) Reverse primer
# 4) primer name
# 5) prefix of output
# 6) Number of threads/cores
# 7) Minimum amplicon lenght
# 8) Maximum amplicon lenght
# -------------
# Requirements:
# -------------
# 1) cutadapt
# 2) usearch
# 3) vsearch
# 4) pear
# 5) seqkit
# 6) fastqc
# 7) GNU parallel

# some functions
length_stats(){
python - << EOF
from collections import Counter
lengths = {}
name = None
seq = ''
with open("$1") as F:
    for line in F:
        if line.startswith('>'):
            if name is not None:
                lengths[name] = len(seq.strip().replace('\n',''))
            seq = ''
            name = line.strip()
        else:
            seq += line
    if name not in lengths:
        lengths[name] = len(seq.strip().replace('\n',''))
c = Counter(lengths.values())
mc = c.most_common(1)[0][0]
gt = sum([x[1] for x in c.items() if x[0] > mc ])
lt = sum([x[1] for x in c.items() if x[0] < mc ])
line = '\n%d sequences are longer than the mode\n%d sequences are ' \
       'shorter than the mode' % (gt, lt)
sor = sorted(c.items(), key=lambda x: x[0])
with open('${2}_distribution.txt', 'w') as F:
    F.write('\n'.join(['%d:%d' % x for x in sor]))
    F.write(line)
EOF
}

run_QC_single_pair(){
# Argumnets:
# 1) Forward primer
# 2) Reverse primer
# 3) Reverse complement of reverse primer
# 4) Reverse complement of reverse primer
# 5) Path to output directory
# 6) prefix of output
# 7) Prefix of the fileset
# 8) cpus

# run first cutadapt
fileset=`basename ${7}`
sample=`echo ${fileset}|rev|cut -d '.' -f1|rev`
cutadapt -g "${1}" -G "${2}" -a "${3}" -A "${4}" \
-o ${5}/${sample}.1.fastq.gz -p ${5}/${sample}.2.fastq.gz \
-m 130 --match-read-wildcards -q 25 --trim-n -n 2 --untrimmed-output \
${5}/untrimmed.${sample}.fastq --untrimmed-paired-output \
${5}/untrimmed.paired.${sample}.fastq ${7}_R1.fastq.gz \
${7}_R2.fastq.gz > ${5}/${sample}.log1
# Merge reads
pear -f ${5}/${sample}.1.fastq.gz -r ${5}/${sample}.2.fastq.gz \
-o ${5}/${sample}_pear -q 25 -t 100 -s 2 > ${5}/${sample}_pear.log
# check if adapters still there
seqkit -j ${8} locate -d -p "${1}" ${5}/${sample}_pear.assembled.fastq > ${5}/${sample}.fwd
seqkit -j ${8} locate -d -p "${2}" ${5}/${sample}_pear.assembled.fastq > ${5}/${sample}.rev
# cut the adapter if they remain
if [[ $(wc -l <${5}/${sample}.fwd) -ge 2 || $(wc -l <${5}/${sample}.rev) -ge 2 ]]; then
cutadapt -a "$3" -m 100 -n 2 --too-short-output ${5}/${sample}.short.fastq \
-o ${5}/${sample}.3trimmed.fastq ${5}/${sample}_pear.assembled.fastq > ${5}/${sample}.log2
cutadapt -g "${1}" -n 2 -o ${5}/${sample}.5trimmed.fastq \
${5}/${sample}.3trimmed.fastq > ${outdir}/${sample}.log3; else
cp ${5}/${sample}_pear.assembled.fastq ${5}/${sample}.3trimmed.fastq
fi
# dereplicate
usearch -fastx_uniques ${5}/${sample}.3trimmed.fastq -fastaout \
${5}/${sample}.trimmed.derep.fasta -sizeout
# get lenghth distributions
length_stats ${5}/${sample}.trimmed.derep.fasta ${5}/${6}
# get stats for all the steps
seqkit -j ${8} stats ${5}/*${sample}*.fa* > ${5}/${6}.stats
# run a fastqc
fastqc ${5}/${sample}.3trimmed.fastq -o ${5}
}

# Get the command line arguments
file_list=$1
ADAPTER_FWD=$2
ADAPTER_REV=$3
primer_name=$4
prefix=$5
cpus=$6
min_len=$7
max_len=$8
outdir=${prefix}_output${primer_name}

echo `date`
echo "executing $0 $@"
# Process the primers and get reverse complements
Adapter1rc=`echo $ADAPTER_FWD | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
Adapter2rc=`echo $ADAPTER_REV | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
# create output folder
mkdir -p ${outdir}

# run QC in all pairs
while read line
do
    run_QC_single_pair ${ADAPTER_FWD} ${ADAPTER_REV} ${Adapter2rc} \
    ${Adapter1rc} ${outdir} ${prefix} ${line} ${cpus}
done <${file_list}


# Trim the reads to expected lenght
for i in ${outdir}/*3trimmed.fastq
do
    cat ${i} | seqkit -j ${cpus} seq -m ${min_len} \
    -M ${max_len} >${i%%3trimmed.fastq}lengthfilter.fastq
done

# Concatenate all resulting files
cat ${outdir}/*.lengthfilter.fastq > ${outdir}/all.lengthfilter.fastq

# Dereplicate on full set reporting size of identical sequences
vsearch --derep_fulllength  ${outdir}/all.lengthfilter.fastq -sizein -sizeout \
-relabel Uniq -output  ${outdir}/vs_all_lengthfilter.fasta

# Run the unoise3 denoising
usearch -unoise3 ${outdir}/vs_all_lengthfilter.fasta -zotus \
${outdir}/all_lengthfilter_zotus_4.fasta -minsize 4 -tabbedout \
${outdir}/allzotus_lengthfilter_4.txt

# Rename to Zotu
sed -i "s/>O/>Zo/g" ${outdir}/all_lengthfilter_zotus_4.fasta

# reannotate each sample
for i in ${outdir}/*lengthfilter.fastq
do
    usearch -fastx_relabel ${i} -prefix `basename ${i%%.lengthfilter.fastq}`_ \
    -fastqout ${i%%.fastq}_relabel.fastq -keep_annots
done

# Make database with the joined file
usearch -makeudb_usearch ${outdir}/all_lengthfilter_zotus_4.fasta -output \
${outdir}/all_lengthfilter_zotus_4.udb

cat ${outdir}/*lengthfilter_relabel.fastq > ${outdir}/all.lengthfilter_relabel.fastq

# Created the Zotu table with the renamed samples
usearch -otutab ${outdir}/all.lengthfilter_relabel.fastq  -sample_delim _ \
-zotus ${outdir}/all_lengthfilter_zotus_4.udb -otutabout \
${outdir}/all_${primer_name}zotutab_4.txt -mapout \
${outdir}/all_${primer_name}zmap_4.txt

