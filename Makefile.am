AUTOMAKE_OPTIONS=foreign
AM_CXXFLAGS = -fpic -Wall -Wextra -std=c++11

bin_PROGRAMS=qm
qm_SOURCES= test.cpp \
            grid_structures.h \
            grid_structures.cpp \
					  read_energy_force_new.h \
						read_energy_force_new.cpp \
						common.h \
						common.cpp \
						gzstream/gzstream.h \
						#trans_rot_coords.h

qm_LDADD = -lz -L./gzstream -lgzstream
