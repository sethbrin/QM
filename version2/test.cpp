#include <iostream>
#include <getopt.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#ifdef USE_THREADS
#include <queue>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#endif

#include "./qm_interpolation.h"

using std::string;

int main(int argc, char* argv[])
{

  string database_file = "/home/zp/work/waterbox/QM/version2//Dimer_deltaEForceTorque_data_mp2_wtr_wtr.txt.gz";
  string model = "wtr";

  database::EnergeForceDatabase energe(database_file.c_str(), model.c_str());
  energe.read_file();
  QM::QMInterpolation interpolation(model, energe);

  std::vector<std::vector<double>> lhs_point = {
    {-9.7219999999999995, -8.0780000000000012, -10.799999999999999},
     {-9.2196870000000004, -8.5968359999999997, -11.491731},
     {-10.531833000000001, -8.5933530000000005, -10.519677}};
  std::vector<std::vector<double>> rhs_point = {
    {-6.859, -10.388999999999999, -11.500999999999999},
     {-7.3702009999999998, -10.859132000000001, -10.781521999999999},
     {-5.9002359999999996, -10.671218, -11.467464}};
  interpolation.calculate(lhs_point, rhs_point);

  return 0;
}
