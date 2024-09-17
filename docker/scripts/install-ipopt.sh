git clone https://github.com/coin-or/Ipopt.git
cd Ipopt
mkdir build
cd build
../configure --with-precision=single
make -j$(nproc)
make install