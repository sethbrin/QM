#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "read_energy_force_new.h"

using namespace database;

Atom::Atom(std::string& line)
{
  std::vector<std::string> fields;
  boost::trim(line);
  boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

  assert(fields.size() == 5);
  m_name = fields[0];
  m_coord = {boost::lexical_cast<double>(fields[2]),
                         boost::lexical_cast<double>(fields[3]),
                         boost::lexical_cast<double>(fields[4])};
}

void Molecule::addAtom(std::string& line)
{
  m_atoms.push_back(Atom(line));
}


Config::Config():
  m_n1(0), m_n2(0), m_energy(0)
{}

Config::Config(int n1, int n2):
  m_n1(n1), m_n2(n2), m_energy(0)
{}

void Config::setEnergy(double energy)
{
  m_energy = energy;
}
void Config::setForce(std::array<double, 3> ff) {
    m_Fx = ff[0];
    m_Fy = ff[1];
    m_Fz = ff[2];
}

void Config::setTorque(std::array<double, 3> tq) {
    m_Tx = tq[0];
    m_Ty = tq[1];
    m_Tz = tq[2];
}

void Config::setName(std::string name)
{
  m_name = name;
}

double Config::get_prop(std::string name, double w)
{
  if (name == "E") return m_energy;
  else if (name == "Fx") return m_Fx;
  else if (name == "Fy") return m_Fy;
  else if (name == "Fz") return m_Fz;
  else if (name == "Tx") return m_Tx;
  else if (name == "Ty") return m_Ty;
  else if (name == "Tz") return m_Tz;
  else return 0;

  // TODO Here will not reach Now
  //int i = boost::lexical_cast<int>(name.substr(1, 3));
  //int j = boost::lexical_cast<int>(name.substr(4, 3));

  //if (w < 0)
  //{
  //  i = SwapH[i];
  //}
  //if (i < m_n1)
  //{
  //  return m_f
  //}
}

/**
 * EnergeForceDatabase
 *
 */
EnergeForceDatabase::EnergeForceDatabase(const char* filename, const char* fragtype):
  m_filename(filename), m_fragtype(fragtype)
{
  m_pDS = getDataStructure(fragtype);

  for (size_t i = 0; i < m_pDS->m_R_NDX.size(); i++) {
      m_all_config.push_back({});
      for (size_t j = 0; j < static_cast<size_t>(m_pDS->m_nGrid); j++) {
          m_all_config[i].push_back({});
          for (size_t m = 0; m < static_cast<size_t>(m_pDS->m_nConf[i]); m++) {
              m_all_config[i][j].push_back({});
              for (size_t k = 0; k < static_cast<size_t>(m_pDS->m_nNorm[i]); k++) {
                  m_all_config[i][j][m].push_back(NULL);
              }
          }
      }
  }

  m_natoms1 = m_pDS->m_n1;
  m_natoms2 = m_pDS->m_n2;
}

EnergeForceDatabase::~EnergeForceDatabase()
{
    for (size_t i = 0; i < m_pDS->m_R_NDX.size(); i++) {
      for (size_t j = 0; j < static_cast<size_t>(m_pDS->m_nGrid); j++) {
          for (size_t m = 0; m < static_cast<size_t>(m_pDS->m_nConf[i]); m++) {
              for (size_t k = 0; k < static_cast<size_t>(m_pDS->m_nNorm[i]); k++) {
                  delete m_all_config[i][j][m][k];
              }
          }
      }
  }

}

Config*& EnergeForceDatabase::at_all(int dim1, int dim2, int dim3, int dim4)
{
    if (dim1 < 0) {
        dim1 = m_pDS->m_R_NDX.size() - std::abs(dim1);
    }
    if (dim2 < 0) {
        dim2 = m_pDS->m_nGrid - std::abs(dim2);
    }
    if (dim3 < 0) {
        dim3 = m_pDS->m_nConf[dim1] - std::abs(dim3);
    }
    if (dim4 < 0) {
        dim4 = m_pDS->m_nNorm[dim1] - std::abs(dim4);
    }
    return m_all_config[dim1][dim2][dim3][dim4];
}


double EnergeForceDatabase::get_prop(int i, int j, int si, int ni, std::string name, double w, double eheigh=100.0)
{
  if (i < 0) return eheigh;
  else if (i > int(m_pDS->m_R_NDX.size()) - 1) return 0.0;

  Config*& cfg = at_all(i, j, si, ni);
  if (!cfg)
  {
    return eheigh;
  }
  else
  {
    return cfg->get_prop(name, w);
  }
}

// the ind only equal to -1 or zero
void EnergeForceDatabase::read_mole(int ind, igzstream& datafile)
{
  std::string line;
  std::getline(datafile, line);

  if (!boost::starts_with(line, "# No."))
  {
    return;
  }
  assert(ind == -1 || ind == 0);
  std::vector<std::string> fields;
  boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);
  if (fields.back() != boost::lexical_cast<std::string>(ind))
  {
    std::cerr << "Mole "<< ind + 1 <<" not found \n";
    exit(EXIT_FAILURE);
  }

  // skip latter three lines
  std::getline(datafile, line);
  std::getline(datafile, line);
  std::getline(datafile, line);

  boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

  if (ind == -1)
  {
    m_mole1.m_energy = boost::lexical_cast<double>(fields[1]);
  }
  else
  {
    m_mole2.m_energy = boost::lexical_cast<double>(fields[1]);
  }
  // skip latter three lines
  std::getline(datafile, line);
  std::getline(datafile, line);
  std::getline(datafile, line);

  while (true)
  {
    std::getline(datafile, line);
    if (boost::starts_with(line, "# End of No."))
    {
      break;
    }
    if (ind == -1)
    {
      m_mole1.addAtom(line);
    }
    else
    {
      m_mole2.addAtom(line);
    }
  }
  // skip blank line
  std::getline(datafile, line);

}

void EnergeForceDatabase::read_file()
{
  igzstream datafile(m_filename.c_str());
  if ( ! datafile.good()) {
    std::cerr << "ERROR: Opening file '" << m_filename << "' failed.\n";
    exit(EXIT_FAILURE);
  }


  if (!datafile.eof())
  {
    // Read Mole 1 (config No.-1):
    read_mole(-1, datafile);
    // Read Mole 2 (config No.0):
    read_mole(0, datafile);
  }

  std::string line;
  std::vector<std::string> fields;
  // read configs
  while (!datafile.eof())
  {
    std::getline(datafile, line);  // No
    std::getline(datafile, line);  // filename

    if (line == "")
    {
      break;
    }
    std::string name = line.substr(2);
    // remove the last .inp.log from line
    line = boost::trim_right_copy_if(name, boost::is_any_of(".inp.log"));
    boost::split(fields, line , boost::is_any_of("_"), boost::token_compress_on);
    int d, g, c, ic;
    std::vector<double>& R_NDX = m_pDS->m_R_NDX;
    d = std::find(R_NDX.begin(), R_NDX.end(), boost::lexical_cast<double>(fields[1].substr(1))) - R_NDX.begin();
    g = boost::lexical_cast<double>(fields[2].substr(1));
    c = boost::lexical_cast<double>(fields[3].substr(1));

    ic = std::floor(c / m_pDS->m_nConf[d]);
    c = c % m_pDS->m_nConf[d];
    //skip blank
    std::getline(datafile, line);
    std::getline(datafile, line); // energy

    if (line.find("Failed") != std::string::npos)
    {
      // skip three line
      std::getline(datafile, line);
      std::getline(datafile, line);
      std::getline(datafile, line);
      continue;
    }
    Config* cfg = new Config(m_natoms1, m_natoms2);
    cfg->setName(name);
    boost::split(fields, line , boost::is_any_of(" "), boost::token_compress_on);

    cfg->setEnergy(boost::lexical_cast<double>(fields[1]));
    std::getline(datafile, line);

    // Force
    std::getline(datafile, line);
    boost::split(fields, line , boost::is_any_of(" "), boost::token_compress_on);
    cfg->setForce({boost::lexical_cast<double>(fields[1]),
            boost::lexical_cast<double>(fields[2]),
            boost::lexical_cast<double>(fields[3])});

    // Torque
    std::getline(datafile, line);
    boost::split(fields, line , boost::is_any_of(" "), boost::token_compress_on);
    cfg->setTorque({boost::lexical_cast<double>(fields[1]),
            boost::lexical_cast<double>(fields[2]),
            boost::lexical_cast<double>(fields[3])});

    // coords & forces of mole2
    for (int i=0; i<m_natoms2; i++)
    {
      std::getline(datafile, line);
      cfg->m_xmole2.push_back(Atom(line));
    }

    std::getline(datafile, line);
    std::getline(datafile, line);

    m_all_config[d][g][c][ic] = cfg;
  }
  datafile.close();
}


std::pair<double, double> database::weights_for_2_configs(const vector<double>& norm_vec, const vector<database::Atom> config1, const vector<database::Atom> config2, double cut)
{
  std::vector<double> vec1 = get_normal_unit(subtraction(config1[1].m_coord, config1[0].m_coord),
                                             subtraction(config1[2].m_coord, config1[0].m_coord));

  std::vector<double> vec2 = get_normal_unit(subtraction(config2[1].m_coord, config2[0].m_coord),
                                             subtraction(config2[2].m_coord, config2[0].m_coord));

  // treat vec1 as the x-axis, the normal of (vec1,vec2) as the z-axis
  std::vector<double>& vx = vec1;
  std::vector<double> vz = get_normal_unit(vec1, vec2);
  // np.cross
  std::vector<double> vy = {vz[1]*vx[2] - vz[2]*vx[1], vz[2]*vx[0]-vz[0]*vx[2], vz[0]*vx[1]-vz[1]*vx[0]};

  // project the "norm_vec" to the plane of vec1 and vec2, then use symmetry operation.
  double px = dot(norm_vec, vx);
  double py = dot(norm_vec, vy);

  double w1, w2;
  if (std::fabs(px) < cut)
  {
    if (py > 0) w2 = 1.0;
    else w2 = -1.0;
    w1 = 0.0;
  }
  else
  {
    double ang = atan(py / px);

    if (px < 0) ang += M_PI;
    else if (py < 0) ang += 2*M_PI;

    double period = M_PI;

    int nrot = ang / period;
    double rrot = ang - period * nrot;

    double ang2;
    if ((rrot - M_PI/2.0)>0.0)
    {
      ang2 = M_PI - rrot;
    }
    else
    {
      ang2 = rrot;
    }

    w2 = ang2 / (M_PI / 2);
    w1 = 1 - w2;
  }

  return {w1, w2};
}
