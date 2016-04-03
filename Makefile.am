AUTOMAKE_OPTIONS=foreign
AM_CXXFLAGS = -fpic -Wall -Wextra -Wno-unused-parameter -std=c++11

bin_PROGRAMS=qm
qm_SOURCES= main.cpp \
            grid_structures.cpp \
						read_energy_force_new.cpp \
						common.cpp \
						qm_interpolation.cpp \
						gzstream/gzstream.h \
						likely/TriCubicInterpolator.h\
						#trans_rot_coords.h

qm_LDADD = -lz -L./gzstream -lgzstream
