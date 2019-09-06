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
# 9) unoise minimum abundance parameter
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

#set -e
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
# 9) minimum read lenght

# run first cutadapt
fileset=`basename ${7}`
sample=`echo ${fileset}|rev|cut -d '.' -f1|rev`
echo "Processing  sample ${sample}"
echo -e "\tRemoving primers"
if [[ -s ${5}/${sample}.1.fastq.gz ]]
then
    echo -e "File ${5}/${sample}.1.fastq.gz exist and is not empty... Skipping"
else
    cutadapt -g "${1}" -G "${2}" -a "${3}" -A "${4}" \
    -o ${5}/${sample}.1.fastq.gz -p ${5}/${sample}.2.fastq.gz \
    -m $9 --match-read-wildcards -q 25 --trim-n -n 2 --untrimmed-output \
    ${5}/untrimmed.${sample}.fastq --untrimmed-paired-output \
    ${5}/untrimmed.paired.${sample}.fastq ${7}_R1.fastq.gz \
    ${7}_R2.fastq.gz > ${5}/${sample}.log1
fi
# Merge reads
echo "Running cutadapt"
pear -f ${5}/${sample}.1.fastq.gz -r ${5}/${sample}.2.fastq.gz \
-o ${5}/${sample}_pear -q 25 -t $9 -s 2 > ${5}/${sample}_pear.log
echo "Checking if any adapters remain"
if [[ -s ${5}/${sample}_pear.assembled.fastq ]]
then
    echo -e "\t\tFile ${5}/${sample}_pear.assembled.fastq exist and is not empty... Skipping"
else
    echo -e "\tMerging reads"
    pear -f ${5}/${sample}.1.fastq.gz -r ${5}/${sample}.2.fastq.gz \
    -o ${5}/${sample}_pear -q 25 -t $9 -s 2 > ${5}/${sample}_pear.log
fi
# check if adapters still there
echo -e "\tIdentifying remnant primers"
seqkit -j ${8} locate -d -p "${1}" ${5}/${sample}_pear.assembled.fastq > ${5}/${sample}.fwd
seqkit -j ${8} locate -d -p "${2}" ${5}/${sample}_pear.assembled.fastq > ${5}/${sample}.rev
# cut the adapter if they remain
if [[ $(wc -l <${5}/${sample}.fwd) -ge 2 || $(wc -l <${5}/${sample}.rev) -ge 2 ]]
then
    echo -e "\t\tRemnant primers found, attempting to remove them"
    cutadapt -a "$3" -m $9 -n 2 --too-short-output ${5}/${sample}.short.fastq \
    -o ${5}/${sample}.3trimmed.fastq ${5}/${sample}_pear.assembled.fastq > ${5}/${sample}.log2
    cutadapt -g "${1}" -n 2 -o ${5}/${sample}.5trimmed.fastq \
    ${5}/${sample}.3trimmed.fastq > ${outdir}/${sample}.log3
else
    echo -e "\t\tNo remnant primers found"
    cp ${5}/${sample}_pear.assembled.fastq ${5}/${sample}.3trimmed.fastq
fi
echo "Dereplicating"
usearch -fastx_uniques ${5}/${sample}.3trimmed.fastq -fastaout \
${5}/${sample}.trimmed.derep.fasta -sizeout
echo "Extimating lenghth distributions"
length__stats ${5}/${sample}.trimmed.derep.fasta ${5}/${6}
echp "estimating stats for all the steps"
seqkit -j ${8} stats ${5}/*${sample}*.fa* >> ${5}/${6}.stats
# dereplicate
echo -e "\tRemoving duplicate sequences"
if [[ -s ${5}/${sample}.trimmed.derep.fasta ]]
then
    echo -e "\t\tFile${5}/${sample}.trimmed.derep.fasta exist and is not empty... Skipping"
else
    fs=`du -b ${5}/${sample}.3trimmed.fastq|cut -f1`
    fs=$(( fs / 1000000000 ))
    if [[ "fs" -gt 3 ]]
    then
        p=`python -c "from math import ceil; print(int(ceil($fs/3))"`
        seqkit -j ${8} split2 ${5}/${sample}.3trimmed.fastq -p ${p} -f \
        -O split_trimmed
        for i in split_trimmed/*part_*
        do
            usearch -fastx_uniques ${i} -fastaout split_trimmed/${i}_uniq \
            -sizeout -threads ${8}
        done
        cat split_trimmed/*_uniq > tmp
        usearch -fastx_uniques tmp -fastaout ${5}/${sample}.trimmed.derep.fasta \
        -sizein -sizeout -threads ${8}
    else
        usearch -fastx_uniques ${5}/${sample}.3trimmed.fastq -fastaout \
        ${5}/${sample}.trimmed.derep.fasta -sizeout -threads ${8}
    fi
fi
# get lenghth distributions
length_stats ${5}/${sample}.trimmed.derep.fasta ${5}/${6}
# get stats for all the steps
seqkit -j ${8} stats ${5}/*${sample}*.fa* > ${5}/${6}.stats
# run a fastqc
echo "Running fastqc"
fastqc ${5}/${sample}.3trimmed.fastq -o ${5}
#fastqc ${5}/${sample}.3trimmed.fastq -o ${5}
echo -e "\n\n"}

prog() {
    # Progress bar, courtesy of ilkkachu (stackoverflow)
    local w=80 p=$1;  shift
    printf -v dots "%*s" "$(( $p*$w/100 ))" ""; dots=${dots// /.};
    printf "\r\e[K|%-*s| %3d %% %s" "$w" "$dots" "$p" "$*";
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
min_read_len=$(( (min_len / 2) - 20 ))
usearch_min_size=${9}
echo `date`
echo -e "executing $0 $@\n"
# Process the primers and get reverse complements
Adapter1rc=`echo $ADAPTER_FWD | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
Adapter2rc=`echo $ADAPTER_REV | tr 'ACGTYRSWKMBDHV' 'TGCARYSWMKVHDB' | rev`
# create output folder
mkdir -p ${outdir}

# run QC in all pairs
while read line
do
    run_QC_single_pair ${ADAPTER_FWD} ${ADAPTER_REV} ${Adapter2rc} \
    ${Adapter1rc} ${outdir} ${prefix} ${line} ${cpus} ${min_read_len}
done <${file_list}

# Trim the reads to expected lenght
echo "Trimming all samples between ${min_len} and ${max_len}"
n=`ls ${outdir}/*3trimmed.fastq| wc -l`
perc=0
for i in ${outdir}/*3trimmed.fastq
do
    let perc++
    prog $(( (perc * 100) / n ))
    cat ${i} | seqkit -j ${cpus} seq -m ${min_len} \
    -M ${max_len} >${i%%3trimmed.fastq}lengthfilter.fastq
done

echo -e "\nConcatenating all resulting files"
# Concatenate all resulting files
cat ${outdir}/*.lengthfilter.fastq > ${outdir}/all.lengthfilter.fastq

# Dereplicate on full set reporting size of identical sequences
echo -e "\nDereplicating full set reporting size of identical sequences"
vsearch --derep_fulllength  ${outdir}/all.lengthfilter.fastq -sizein -sizeout \
-relabel Uniq -output ${outdir}/vs_all_lengthfilter.fasta

# execute denoising in full file
echo -e "\nDenoising"
usearch -unoise3 ${outdir}/vs_all_lengthfilter.fasta -zotus \
${outdir}/all_lengthfilter_zotus_${usearch_min_size}.fasta \
-minsize ${usearch_min_size} -tabbedout \
${outdir}/allzotus_lengthfilter_${usearch_min_size}.txt

# Rename to Zotu
echo -e "\n Renaming sequences"
sed -i "s/>O/>Zo/g" ${outdir}/all_lengthfilter_zotus_${usearch_min_size}.fasta

# reannotate each sample
echo -e "\n Re-annotate each sample"
n=`ls ${outdir}/*lengthfilter.fastq| wc -l`
perc=0
for i in ${outdir}/*lengthfilter.fastq
do
    let perc++
    prog $(( (perc * 100) / n ))
    usearch -fastx_relabel ${i} -prefix `basename ${i%%.lengthfilter.fastq}`_ \
    -fastqout ${i%%.fastq}_relabel.fastq -keep_annots
done

# Make database with the joined file
echo -e "\n\nCreating usearch database"
usearch -makeudb_usearch ${outdir}/all_lengthfilter_zotus_${usearch_min_size}.fasta \
-output ${outdir}/all_lengthfilter_zotus_${usearch_min_size}.udb

# Concatenating relabeled sequences
echo -e "\nConcatenating relabeled sequences"
cat ${outdir}/*lengthfilter_relabel.fastq > ${outdir}/all.lengthfilter_relabel.fastq

# get file_size
echo -e "\nGenerating ZOTU/ASV  table"
file_size=`du -b ${outdir}/all.lengthfilter_relabel.fastq|cut -f1`
file_size=$(( file_size / 1000000000 ))
if [[ "$file_size" -gt 3 ]]
    then
        # split file into 3GB chunks
        p=`python -c "from math import ceil; print(int(ceil(${file_size}/3)))"`
        seqkit split2 ${outdir}/all.lengthfilter_relabel.fastq -p ${p} -f \
        -O split
        # execute in smaller files
        for i in split/*part_*
        do
            dir=`dirname ${i}`
            part=`basename ${i}|rev| cut -d'.' -f 2|rev`
            usearch -otutab ${i}  -sample_delim _ \
            -zotus ${outdir}/all_lengthfilter_zotus_${usearch_min_size}.udb \
            -otutabout ${dir}/${part}_zotutab_${usearch_min_size}.txt \
            -mapout ${dir}/${part}zmap_${usearch_min_size}.txt
        done
        # merge them
        # merge the individual otu tables
        files=`echo split/*zotutab_${usearch_min_size}.txt|tr ' ' ','`
        usearch -otutab_merge ${files} -output ${outdir}/allzotus_lengthfilter_${usearch_min_size}.txt
    else
    # execute full
    # Created the Zotu table with the renamed samples
    usearch -otutab ${outdir}/all.lengthfilter_relabel.fastq  -sample_delim _ \
    -zotus ${outdir}/all_lengthfilter_zotus_${usearch_min_size}.udb -otutabout \
    ${outdir}/all_${primer_name}zotutab_${usearch_min_size}.txt -mapout \
    ${outdir}/all_${primer_name}zmap_${usearch_min_size}.txt
fi
