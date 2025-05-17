git clone https://github.com/ompl/ompl.git
cd ompl
git checkout tags/1.7.0
mkdir build
cd build
cmake ..
make -j$(nproc) clean_bindings
make -j$(nproc) update_bindings 
make -j$(nproc) py_ompl
make -j$(nproc)
make install