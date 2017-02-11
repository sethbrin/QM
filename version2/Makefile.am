AUTOMAKE_OPTIONS=foreign
AM_CXXFLAGS = -fpic -Wall -Wextra -Wno-unused-parameter -std=c++11 -O2 #-DUSE_THREADS

bin_PROGRAMS=qm
qm_SOURCES= main.cpp \
            grid_structures.cpp \
            read_energy_force_new.cpp \
            common.cpp \
            qm_interpolation.cpp \
            gzstream/gzstream.h

qm_LDADD = -lz -L./gzstream -lgzstream -lboost_filesystem -lboost_system -lboost_thread
