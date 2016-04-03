#include "./grid_structures.h"
#include "./read_energy_force_new.h"
#include "./common.h"
#include "./qm_interpolation.h"

#include <iostream>
#include <boost/algorithm/string.hpp>

int main()
{
  database::EnergeForceDatabase energe("./Dimer_deltaE_data_mp2_wtr_wtr.txt.gz", "wtr");
  energe.read_file();

  QM::QMInterpolation interpolation("wtr", energe);

  std::cout << interpolation.process("./random/test0002.inp") << std::endl;
  return 0;
}
