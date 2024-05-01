unzip HSL.zip
cd HSL
mkdir build
cd build
../configure --prefix=/usr/local
make -j$(nproc)
make install