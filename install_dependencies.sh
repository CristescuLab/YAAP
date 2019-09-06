#!/usr/bin/env bash
cwd=${PWD}
echo "Preparing cutadapt"
python -m pip uninstall cutadapt
python -m pip --user 'cutadapt==1.18'

echo "Preparing vsearch"
git clone https://github.com/torognes/vsearch.git
cd vsearch
./autogen.sh
./configure --prefix=${PWD} --exec-prefix=${PWD}
make
cd bin
chmod +x vsearch
cp vsearch ${cwd}/executables_linux_64/
cd ${cwd}

echo "Setting up seqkit, pear and usearch "
mv $1 executables_linux_64/usearch
cd executables_linux_64
chmod +x *
echo "export PATH=$PATH:$PWD" >> ${HOME}/.bashrc
source ${HOME}/.bashrc
. ${HOME}/.bashrc
cd ${cwd}
