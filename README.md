QM
----------------
energy force calculation, current only support waterbox, You can set two mode,
use thread or not use thread, just edit `Makefile.am`,you can see such line:

```
AM_CXXFLAGS = -fpic -Wall -Wextra -Wno-unused-parameter -std=c++11 -O2 #-DUSE_THREADS
```
To add thread support, just uncomment the `-DUSE_THREADS`

dependency
----------------
1. `g++` g++ version should support c++11
2. `zlib` current only test for version 1.2.8
3. `boost` current only test for the newest boost, 1.55, you can `yum install boost` to install boost,
or just store in BOOST_PATH, And add include and lib to the .bashrc.
```
export LD_LIBRARY_PATH=$BOOST_PATH/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$BOOST_PATH/lib:$LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$BOOST_PATH/include:$CPLUS_INCLUDE_PATH
```

install
----------------
1. use following method to compile
```
aclocal && automake --add-missing && autoconf
./configure
make
```

Usage
----------------
You can see use `qm -h` to see more infomation

