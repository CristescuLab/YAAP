# YAAP
Yet Another Amplicon denoising Pipeline (YAAP), is a pipeline to analyse 
metabarcoding amplicon data. It performs QC, adapter removal, denoising and 
ZOTU table construction

Table of Contents
=================

   * [YAAP](#yaap)
      * [Dependencies](#dependencies)
         * [Installation of dependencies](#installation-of-dependencies)
            * [<strong>The easier (no gurarantees) way</strong>](#the-easier-no-gurarantees-way)
         * [<strong>The not so easy, but easy way</strong>](#the-not-so-easy-but-easy-way)
            * [Linking files to your bin directory (REQUIRES ADMIN PRIVILEGES )](#linking-files-to-your-bin-directory-requires-admin-privileges-)
            * [Exporting the path](#exporting-the-path)
      * [Usage](#usage)
      * [ISSUES](#issues)

## Dependencies
This pipeline makes use of the following programs:
1. cutadapt ([https://cutadapt.readthedocs.io/en/stable/index.html](
https://cutadapt.readthedocs.io/en/stable/index.html))
2. vsearch ([https://github.com/torognes/vsearch](
https://github.com/torognes/vsearch))
3. seqkit ([https://bioinf.shenwei.me/seqkit/](
https://bioinf.shenwei.me/seqkit/))
4. pear ([https://cme.h-its.org/exelixis/web/software/pear/](
https://cme.h-its.org/exelixis/web/software/pear/))
5. usearch ([https://www.drive5.com/usearch/](
https://www.drive5.com/usearch/))

This git contains the binaries (executables) of pear, seqkit, and usearch, 
however, vsearch and cutadapt have to be installed locally.

### Installation of dependencies
All the instructions here are only tested in linux systems. You might have to 
adjust part of it for other operating systems.

#### **The easier (no gurarantees) way**
I recently added a bash script to install all dependencies for you. It might 
not work in all systems, but it might be worth given a try. This script assumes
that you have python3 in your environment ans can be executed as python. If you
are in a Compute canada cluster, you need to load the `scipy-stack/2018b` 
module. So let's assume that you want me to walk you though all the 
installation. Then you need to:
1. Clone this repository (you need to have git installed):
```bash
git clone https://github.com/CristescuLab/YAAP.git
```
2. Get into the folder:
```bash
cd YAAP
```
2.  Because of the licence in usearch I am not allowed to distribute the binary.
So you will need to visit [https://drive5.com/usearch/download.html](
https://www.drive5.com/usearch/) download 
the latest LINUX version from your email, and place it in this folder (the YAAP 
folder)

4. Execute the `install_dependencies.sh` code:
```bash
bash install_dependencies.sh <name of usearch binary>
```
You will need to change `<name of usearch binary>` for the actual name of the 
binary you downloaded from usearch 
5. Test that all dependencies work:
```bash
cutadapt -h
vsearch -h
pear -h
seqkit -h
usearch 
```
If all of the above commands do not give you errors, you are good to go!! If you still have errors, try:
```bash
source ~/.bashrc
```
And then, try again.

### **The not so easy, but easy way**
1. cutadapt: Follow the instructions at [
https://cutadapt.readthedocs.io/en/stable/installation.html](
https://cutadapt.readthedocs.io/en/stable/installation.html). If you are using 
a computer cluster I advice you AGAINST installing cutadapt with conda.
2. vsearch: Follow the instruction at [https://github.com/torognes/vsearch](
https://github.com/torognes/vsearch)

For the rest of the dependencies, you need to point your path to the executable
directory or copy/link the executables to your bin directory.  Alternatively, 
you can install them from scratch in your system.

#### Linking files to your bin directory (REQUIRES ADMIN PRIVILEGES )
If you don't have administrator privileges, go to [the next section](#exporting-the-path). 
If you do, you can just copy (`cp`) or symbolically link (`ln -s`). 
For example, let's say that you cloned YAAP to your home directory and let's 
assume that your username is `user1`. If you are in a linux system you can just
 type:

```
cd /home/user1/YAAP/executables_linux_64/
sudo ln -s * /usr/local/bin
```

Then pear, seqkit, and usearch are available to the pipeline.


#### Exporting the path
If you don't have administrator privileges, you can let the system know to look
 for the executables in the appropriate folder. As before, let's assume that 
 your username is `user1` and that you clone the repository in home. You can 
 export the path by typing:
```
cd /home/user1/YAAP/executables_linux_64/
echo "export PATH=$PATH:$PWD" >> ~/.bashrc
source ~/.bashrc
```

## Usage
YAAP has a single bash script called ASV_pipeline.sh. It requires 8 arguments:
1. File with a list of fileset prefix (one per line)
2. Forward primer sequence
3. Reverse primer sequence
4. Primer name
5. Prefix of the output
6. Number of cpus to use
7. Minimum amplicon length
8. Maximum amplicon length

This pipeline assumes that your files are names fileset_prefix_R1.fastq.gz and 
fileset_prefix_R2.fastq.gz for all the samples. It also assumes that you have 
demultiplexed your samples.

As an example, let's assume that we have two samples named MC1 and MC2. Your 
demultiplexed fastq files are `MI.M03555_0320.001.N701N517.MC1_R1.fastq.gz`, 
`MI.M03555_0320.001.N701N517.MC1_R2.fastq.gz`, 
`MI.M03555_0320.001.N702N517.MC2_R1.fastq.gz`,
`MI.M03555_0320.001.N702N517.MC2_R1.fastq.gz`.
 You will need to create a file that looks like this:
 ```
MI.M03555_0320.001.N701N517.MC1
MI.M03555_0320.001.N702N517.MC2 
 ```
Alternatively, you can also create the with the path where the file sets are. 
Let's assume that you called this file `file_list.txt`. Let's also assumed that
 you are working with the Leray fragments flanked by 
 `GGWACWGGWTGAACWGTWTAYCCYC` as forward, and `TAAACTTCAGGGTGACCAAAAAATCA` as 
 reverse, and in a machine that have 28 cpus. You can run the pipeline as:
```
bash ~/YAAP/ASV_pipeline.sh file_list.txt \
GGWACWGGWTGAACWGTWTAYCCYC TAAACTTCAGGGTGACCAAAAAATCA \
COI test 28 293 333 4
```

This call of the pipeline will focus on reads that contained the primer 
sequences, that were longer than 126 (we focused on reads that are longer than 
(min_len/2) - 20) in either direction (forwards and reverse), and which 
assembly shoud be between 293bp and 333 (around expected length of the Leray 
amplicon). The denoising would filter out reads with lower abundances than 4.
It will create an output folder with the name of the primer (COI in 
the example) and append test to the output files. In the example, the output 
folder will be called `test_outputCOI`.


## ISSUES
If you find any bugs related with the pipeline, please open an issue in the the
github repo, and add some traceback (error) information.
