#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "read_energy_force_new.h"

using namespace database;

Atom::Atom(std::string line)
{
  std::vector<std::string> fields;
  boost::trim(line);
  boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

  assert(fields.size() == 5);
  m_name = fields[0];
  m_coord = Coord(boost::lexical_cast<double>(fields[2]),
                         boost::lexical_cast<double>(fields[3]),
                         boost::lexical_cast<double>(fields[4]));
}

void Molecule::addAtom(std::string line)
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

void Config::setName(std::string name)
{
  m_name = name;
}

EnergeForceDatabase::EnergeForceDatabase(const char* filename, const char* fragtype):
  m_filename(filename), m_fragtype(fragtype)
{
  m_pDS = getDataStructure(fragtype);

  m_all_config = new Config*[m_pDS->m_R_NDX.size() * m_pDS->m_nGrid * m_pDS->m_nConf * m_pDS->m_nNorm];
  m_wtr_config = new Config*[m_pDS->m_nConf * m_pDS->m_nNorm];

  for (size_t i=0; i<m_pDS->m_R_NDX.size() * m_pDS->m_nGrid * m_pDS->m_nConf * m_pDS->m_nNorm; i++)
  {
    m_all_config[i] = NULL;
  }
  for (int i=0; i<m_pDS->m_nConf * m_pDS->m_nNorm; i++)
  {
    m_wtr_config[i] = NULL;
  }

  m_natoms1 = m_pDS->m_n1;
  m_natoms2 = m_pDS->m_n2;
}

EnergeForceDatabase::~EnergeForceDatabase()
{
  for (size_t i=0; i<m_pDS->m_R_NDX.size() * m_pDS->m_nGrid * m_pDS->m_nConf * m_pDS->m_nNorm; i++)
  {
    delete m_all_config[i];
  }
  delete []m_all_config;

  //for (int i=0; i<m_pDS->m_nConf * m_pDS->m_nNorm; i++)
  //{
  //  delete m_wtr_config[i];
  //}

  delete []m_wtr_config;
}

Config*& EnergeForceDatabase::at_all(int dim1, int dim2, int dim3, int dim4)
{
  int size1 = m_pDS->m_R_NDX.size();
  int size2 = m_pDS->m_nGrid;
  int size3 = m_pDS->m_nConf;
  int size4 = m_pDS->m_nNorm;
  assert(dim1 < size1);
  assert(dim2 < size2);
  assert(dim3 < size3);
  assert(dim4 < size4);

  return m_all_config[dim1 * size2 * size3 * size4 + dim2 * size3 * size4 + dim3 * size4 + dim4];
}

Config*& EnergeForceDatabase::at_wtr(int dim1, int dim2)
{
  int size1 = m_pDS->m_nConf;
  int size2 = m_pDS->m_nNorm;
  assert(dim1 < size1);
  assert(dim2 < size2);

  return m_wtr_config[dim1 * size2 + dim2];
}

void EnergeForceDatabase::set_wtr(int dim1, int dim2, Config* cfg)
{
  int size1 = m_pDS->m_nConf;
  int size2 = m_pDS->m_nNorm;
  assert(dim1 < size1);
  assert(dim2 < size2);

  m_wtr_config[dim1 * size2 + dim2] = cfg;

}
void EnergeForceDatabase::set_all(int dim1, int dim2, int dim3, int dim4, Config* cfg)
{
  int size1 = m_pDS->m_R_NDX.size();
  int size2 = m_pDS->m_nGrid;
  int size3 = m_pDS->m_nConf;
  int size4 = m_pDS->m_nNorm;
  assert(dim1 < size1);
  assert(dim2 < size2);
  assert(dim3 < size3);
  assert(dim4 < size4);

  m_all_config[dim1 * size2 * size3 * size4 + dim2 * size3 * size4 + dim3 * size4 + dim4] = cfg;
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
  //skip blank line
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

    if (c > 25)
    {
      ic = 1;
      c = c - 26;
    }
    else
    {
      ic = 0;
    }

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

    // coords & forces of mole2
    for (int i=0; i<m_natoms2; i++)
    {
      std::getline(datafile, line);
      cfg->m_xmole2.push_back(Atom(line));
    }

    double H_E = H_energy(*cfg);
    cfg->setEnergy(cfg->m_energy - H_E);

    std::getline(datafile, line);
    std::getline(datafile, line);

    set_all(d, g, c, ic, cfg);
    if (at_wtr(c, ic) == NULL)
    {
      set_wtr(c, ic, cfg);
    }
  }
  datafile.close();
}


double EnergeForceDatabase::H_energy(Config& config)
{
  double res = 0;

  if (!m_pDS->m_IsHmm) return res;

  for (size_t i=0; i<m_pDS->m_Hmm.size(); i++)
  {
    double eps = m_pDS->m_params["HC-OW"].dim2[0];
    double sig =  m_pDS->m_params["HC-OW"].dim2[1];

    double r = Coord::distance(m_mole1.m_atoms[i].m_coord, config.m_xmole2[0].m_coord);

    res += (m_pDS->*(m_pDS->m_lj))(r, eps, sig);
    double q0 = m_pDS->m_params["qHC"].dim1;
    double q1 = m_pDS->m_params["qOW"].dim1;
    res += m_pDS->coulomb(q0, q1, r);

    eps = m_pDS->m_params["HC-HW"].dim2[0];
    sig =  m_pDS->m_params["HC-HW"].dim2[1];

    for (int j=0; j<3; j++)
    {
      r = Coord::distance(m_mole1.m_atoms[i].m_coord, config.m_xmole2[j].m_coord);
      res += (m_pDS->*(m_pDS->m_lj))(r, eps, sig);
      res += charge(m_pDS->m_params["qHC"].dim1, m_pDS->m_params["qHW"].dim1, r);
    }
  }

  return res;
}

double EnergeForceDatabase::charge(double q0, double q1, double r)
{
  return q0*q1/r*332.5;
}
