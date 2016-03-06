#include "read_energy_force_new.h"
#include "gzstream/gzstream.h"

using namespace database;

EnergeForceDatabase::EnergeForceDatabase(const char* filename, const char* fragtype):
  m_filename(filename), m_fragtype(fragtype)
{
  m_pDS = getDataStructure(fragtype);

  m_all_config = new Config[m_pDS->m_R_NDX.size() * m_pDS->m_nGrid * m_pDS->m_nConf * m_pDS->m_nNorm];
  m_wtr_config = new Config[m_pDS->m_nConf * m_pDS->m_nNorm];
}

EnergeForceDatabase::~EnergeForceDatabase()
{
  
  delete []m_all_config;

  delete []m_wtr_config;
}

Config& EnergeForceDatabase::at_all(int dim1, int dim2, int dim3, int dim4)
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

Config& EnergeForceDatabase::at_wtr(int dim1, int dim2)
{
  int size1 = m_pDS->m_nConf;
  int size2 = m_pDS->m_nNorm;
  assert(dim1 < size1);
  assert(dim2 < size2);

  return m_wtr_config[dim1 * size2 + dim2];
}

void EnergeForceDatabase::read_file()
{
  igzstream datafile(m_filename.c_str());
  if ( ! datafile.good()) {
    std::cerr << "ERROR: Opening file `" << m_filename << "' failed.\n";
     exit(EXIT_FAILURE);
  }

  const int BUFSIZE = 128;
  char line[BUFSIZE];

  while (!datafile.eof())
  {
    // Read Mole 1 (config No.-1):
    datafile.getline(line, BUFSIZE);
  }
  datafile.close();
}
