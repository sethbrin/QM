QM
----------------
energy force calculation, current only support waterbox, You can set two mode,
use thread or not use thread, just edit `Makefile.am`,you can see such line:

```
AM_CXXFLAGS = -fpic -Wall -Wextra -Wno-unused-parameter -std=c++11 -O2 #-DUSE_THREADS
```
To add thread support, just uncomment the `-DUSE_THREADS`

install
----------------

1. `zlib` current only test for version 1.2.8
2. `boost` current only test for the newest boost, 1.55

Usage
----------------
You can see use `qm -h` to see more infomation
