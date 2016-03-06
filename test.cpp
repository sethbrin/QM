#include "./grid_structures.h"
#include "./read_energy_force_new.h"
//#include "./trans_rot_coords.h"

int main()
{
  database::WaterStructure ws("saml");

  database::EnergeForceDatabase energe("./gzstream/version.gz", "as");

  energe.read_file();
  return 0;
}
