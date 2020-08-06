current_path=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

# sra-tools
cd ncbi/ngs
./configure --prefix=./
make -C ngs-sdk
make -C ngs-java
cd ../ncbi-vdb
./configure --prefix=./
make
make install
cd ../sra-tools
./configure --prefix=./
make
make install
# echo "export PATH=${current_path}/ncbi/sra-tools/bin:"'$PATH' >> ~/.bashrc
