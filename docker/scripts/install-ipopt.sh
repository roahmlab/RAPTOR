git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
mkdir build
cd build
../configure
make -j$(nproc)
make install