#!/usr/bin/env bash
cwd=${PWD}
echo "Preparing cutadapt"
python -m pip install --user --upgrade cutadapt --user

echo "Preparing vsearch"
git clone https://github.com/torognes/vsearch.git
cd vsearch
./autogen.sh
./configure --prefix=${PWD} --exec-prefix=${PWD}
make
cd bin
echo "export PATH=$PATH:$PWD" >> ~/.bashrc
cd ${cwd}
echo "Setting up seqkit, pear and usearch "
mv $1 executables_linux_64/usearch
cd executables_linux_64
chmod +x usearch
echo "export PATH=$PATH:$PWD" >> ~/.bashrc
source ~/.bashrc
cd ${cwd}
