#include "./qm_interpolation.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>

using namespace QM;

/**
 * Atom definition
 */
GjfAtom::GjfAtom(std::string line)
{
  std::vector<std::string> fields;
  boost::trim(line);
  boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

  assert(fields.size() == 4);
  m_name = fields[0];
  m_x =  { boost::lexical_cast<double>(fields[1]),
                         boost::lexical_cast<double>(fields[2]),
                         boost::lexical_cast<double>(fields[3])};

}

Atom* GjfAtom::copy()
{
  Atom* newAtom = new GjfAtom(*this);

  return newAtom;
}

InpAtom::InpAtom(std::string line)
{
  std::vector<std::string> fields;
  boost::trim(line);
  boost::split(fields, line, boost::is_any_of(" \t"), boost::token_compress_on);

  assert(fields.size() == 5);
  m_name = fields[0];
  m_x =  { boost::lexical_cast<double>(fields[2]),
                         boost::lexical_cast<double>(fields[3]),
                         boost::lexical_cast<double>(fields[4])};

}

Atom* InpAtom::copy()
{
  Atom* newAtom = new InpAtom(*this);

  return newAtom;
}

// TODO complete it
PdbAtom::PdbAtom(std::string line)
{
}

Atom* PdbAtom::copy()
{
  Atom* newAtom = new PdbAtom(*this);

  return newAtom;
}

/**
 * Coordinates definition
 */
Coordinates::Coordinates(int n1, int n2, std::string& fragtype, std::string& name):
  m_n1(n1), m_n2(n2), m_fragtype(fragtype), m_name(name)
{
  m_pDS = database::getDataStructure(fragtype.c_str());

  m_is_oriented = false;
  m_facendx = {
    {"yx",2},
    {"xy",2},
    {"yz",0},
    {"zy",0},
    {"zx",1},
    {"xz",1},
    {"zarg",5},
    {"zben",6}
  };

  m_symm = {1,1,1};
  m_center = 0;
}

void Coordinates::addAtom(std::string& line, const std::string& ftype)
{
  if (ftype == "gms")
  {
    m_atoms.push_back(boost::shared_ptr<Atom>(new InpAtom(line)));
  }
}

void Coordinates::ReorientToOrigin(double cut=0.0000001)
{
  int cnt_of_atoms = m_atoms.size();
  assert(cnt_of_atoms == 3);
  std::vector<double> dvec = m_pDS->calt_dvec(m_atoms[0]->m_x, m_atoms[1]->m_x, m_atoms[2]->m_x);

  for (int i=0; i<cnt_of_atoms; i++)
  {
    m_atoms[i]->m_x = translate(m_atoms[i]->m_x, dvec);
  }

  // rotate oxygen to x:
  std::pair<std::vector<double>, std::vector<double> > tmp_pair =
    m_pDS->calt_vec1(m_atoms[0]->m_x, m_atoms[1]->m_x, m_atoms[2]->m_x);

  double ang = angle(tmp_pair.first, tmp_pair.second);
  vector<double> ax = get_normal(tmp_pair.first, tmp_pair.second);
  if (ax[0] == 0 && ax[1] == 0 && ax[2] == 0)
  {
    return;
  }

  for (int i=0; i<cnt_of_atoms; i++)
  {
    m_atoms[i]->m_x = rotate(m_atoms[i]->m_x, ax, ang);
  }

  //  rotate C5 to y:
  tmp_pair =
    m_pDS->calt_vec2(m_atoms[0]->m_x, m_atoms[1]->m_x, m_atoms[2]->m_x);

  ang = angle(tmp_pair.first, tmp_pair.second);
  if (std::abs(ang) < cut)
  {
    return;
  }
  else
  {
    if (std::abs(ang - M_PI) < cut)
    {
      ax = {1,0,0};
    }
    else
    {
      ax = get_normal(tmp_pair.first, tmp_pair.second);
    }

    for (int i=0; i<cnt_of_atoms; i++)
    {
      m_atoms[i]->m_x = rotate(m_atoms[i]->m_x, ax, ang);
    }
  }

  m_is_oriented = true;

  spherical_x();
}

//Calculate the coords in spherical coordination system for molecule 2
void Coordinates::spherical_x()
{
   vector<double> x = m_atoms[m_n1]->m_x;
   double r = ::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

   double ang1 = M_PI * 0.5 - ::acos(x[2] / r);
   double ang2;
   if (::abs(x[0] < 0.000001))
   {
     if (x[1] > 0)
     {
       ang2 = M_PI * 0.5;
     }
     else
     {
       ang2 = M_PI * 1.5;
     }
   }
   else
   {
     ang2 = ::atan(x[1] / x[0]);
     if (x[0] < 0)
     {
       ang2 += M_PI;
     }
     else
     {
       ang2 += M_PI * 2;
     }
   }

   m_r = r;
   m_ang1 = ang1;
   m_ang2 = ang2;
   m_center2 = x;
}

/**
 * QMInterpolation definition
 */
const std::vector<std::vector<int> > QMInterpolation::m_aa_ndx = {{1, 2, 3}};
const std::vector<int> QMInterpolation::m_c_ndx = {1};
const std::vector<int> QMInterpolation::m_H_ndx = {2,3,4};
const std::vector<int> QMInterpolation::m_prob_ndx = {4,5,6};

QMInterpolation::QMInterpolation(std::string fftype):
  m_fftype(fftype)
{};

void QMInterpolation::process(std::string filename)
{
  std::ifstream ifs(filename);

  if ( ! ifs.good()) {
    std::cerr << "ERROR: Opening file '" << filename << "' failed.\n";
    exit(EXIT_FAILURE);
  }

  std::string line;
  while (! ifs.eof())
  {
    std::getline(ifs, line);

    if (line == "") {
      std::cerr << "No coord found\n";
      exit(EXIT_FAILURE);
    }

    if (boost::starts_with(line, " $DATA")) {
      // skip later two lines
      std::getline(ifs, line);
      std::getline(ifs, line);
      break;
    }
  }

  std::vector<std::string> lines;
  while (! ifs.eof())
  {
    std::getline(ifs, line);
    boost::trim(line);
    if (line == "" || boost::starts_with(line, " $end") ||
        boost::starts_with(line, " $END"))
    {
      break;
    }

    lines.push_back(line);
  }
  ifs.close();

  for (size_t i=0; i<m_aa_ndx.size(); i++)
  {
    const std::vector<int>& mndx = m_aa_ndx[i];

    char name_str[20];
    sprintf(name_str, "%s%02d", filename.c_str(), mndx[0]);
    std::string name(name_str);
    Coordinates interp = Coordinates(3,3, m_fftype, name);

    for (auto ind : mndx)
    {
      interp.addAtom(lines[ind-1], "gms");
    }

    // Read probe:
    for (auto ind : m_prob_ndx)
    {
      interp.addAtom(lines[ind-1], "gms");
    }
  }


}
