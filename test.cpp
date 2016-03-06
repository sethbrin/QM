#include "./grid_structures.h"
#include "./read_energy_force_new.h"
#include "./common.h"
#include <iostream>
#include <boost/algorithm/string.hpp>
//#include "./trans_rot_coords.h"

int main()
{
//  database::WaterStructure ws("saml");

  database::EnergeForceDatabase energe("./Dimer_deltaE_data_mp2_wtr_wtr.txt.gz", "wtr");

  energe.read_file();


  return 0;
}
