cwd=${PWD}
echo "Preparing usearch"
mv $1 usearch
chmod +x usearch
echo "export PATH=$PATH:$PWD" >> ~/.bashrc

echo "Preparing cutadapt"
python -m pip install --user --upgrade cutadapt --user

echo "Preparing vsearch"
git clone https://github.com/torognes/vsearch.git
cd vsearch
./autogen.sh
./configure --prefix=${PWD} --exec-prefix=${PWD}
make
echo "export PATH=$PATH:$PWD" >> ~/.bashrc
echo "Setting up seqkit, pear and usearch "
cd executables_linux_64
echo "export PATH=$PATH:$PWD" >> ~/.bashrc
source ~/.bashrc
cd ${cwd}
