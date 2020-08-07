current_path=$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)

# sra-tools
cd ncbi/ngs
./configure --prefix=${current_path}/ncbi/ngs
make -C ngs-sdk
make -C ngs-java
cd ../ncbi-vdb
./configure --prefix=${current_path}/ncbi/ncbi-vdb
make
make install
cd ../sra-tools
./configure --prefix=${current_path}/ncbi/sra-tools
make
make install
echo "export PATH=${current_path}/ncbi/sra-tools/bin:"'$PATH' >> ~/.bashrc
