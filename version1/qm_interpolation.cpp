#include "qm_interpolation.h"
#include "TriCubicInterpolator.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <exception>

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
  if (fabs(ang) < cut)
  {
    return;
  }
  else
  {
    if (std::fabs(ang - M_PI) < cut)
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
   vector<double>& x = m_atoms[m_n1]->m_x;
   double r = ::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);

   double ang1 = M_PI * 0.5 - ::acos(x[2] / r);
   double ang2;
   if (fabs(x[0]) < 0.000001)
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
     else if (x[1] < 0)
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
 *According to the coords of the 1st atom in mole2.
 */
void Coordinates::MirrorAll()
{

  vector<double> tmpXZ = {0, 0, -1};
  for (std::string face: m_pDS->m_symface)
  {
    int fndx = m_facendx[face];

    double tempang2;
    bool isMillorXZ = false;
    // Ben
    if (fndx == 6)
    {
      //Rotate to xz (+-)30.0:
      if (m_ang2 < 0)
      {
        tempang2 = m_ang2 + 2.0 * M_PI;
      }
      else
      {
        tempang2 = m_ang2;
      }

      double period = M_PI / 3.0;
      int nrot = tempang2 / period;
      double rrot = tempang2 - period * nrot;

      if (rrot - M_PI / 6.0 > 0.0)
      {
        nrot ++;
        isMillorXZ = true;
      }
      else
      {
        isMillorXZ = false;
      }

      for (size_t i=m_n1; i<m_atoms.size(); i++)
      {
        m_atoms[i]->m_x = rotate(m_atoms[i]->m_x, tmpXZ, nrot * period);
      }
      if (isMillorXZ)
      {
        for (size_t i=m_n1; i<m_atoms.size(); i++)
        {
          m_atoms[i]->m_x[1] *= -1;
        }
      }
      continue;
    }
    else if (fndx == 5)
    {
      //Rotate to xz (+-)60.0:
      if (m_ang2 < 0)
      {
        tempang2 = m_ang2 + 2.0 * M_PI;
      }
      else
      {
        tempang2 = m_ang2;
      }

      double period = M_PI / 3.0 * 2;
      int nrot = tempang2 / period;
      double rrot = tempang2 - period * nrot;

      if (rrot - M_PI / 3.0 > 0.0)
      {
        nrot ++;
        isMillorXZ = true;
      }
      else
      {
        isMillorXZ = false;
      }

      for (size_t i=m_n1; i<m_atoms.size(); i++)
      {
        m_atoms[i]->m_x = rotate(m_atoms[i]->m_x, tmpXZ, nrot * period);
      }
      if (isMillorXZ)
      {
        for (size_t i=m_n1; i<m_atoms.size(); i++)
        {
          m_atoms[i]->m_x[1] *= -1;
        }
      }
      continue;
    }


    // others
    if (m_center2[fndx] < 0.0)
    {
      m_symm[fndx] = -1;

      for (size_t i=m_n1; i<m_atoms.size(); i++)
      {
        m_atoms[i]->m_x[fndx] *= -1;
        m_atoms[i]->m_f[fndx] *= -1;
      }
    }
  }

  spherical_x();
}

void Coordinates::indexing()
{
  if (!m_is_oriented)
  {
    throw std::invalid_argument("Error: indexing beforce reorientation");
  }

  double r = m_r;
  double ang1 = m_ang1;
  double ang2 = m_ang2;

  // ndx of r
  int ir = 10001;
  ir = get_index(r, m_pDS->m_R_NDX);
  if (ir > 10000 || ir < 0)
  {
    m_ir = ir;
    m_ig = 0;
    m_vbis = {0, 0, 0};
    m_vnrm = {0 ,0, 0};
    return;
  }
  int ih = get_index(ang1, m_pDS->m_PHI_angles);

  int ip = -1;
  for (size_t i=1; i<(size_t)m_pDS->m_NTheta[ih]; i++)
  {
    if (ang2 <= m_pDS->m_THETA_angles[ih][i])
    {
      ip = i - 1;
      break;
    }
  }

  if (ip == -1)
  {
    ip = m_pDS->m_NTheta[ih] - 1;
  }

  //ig is the index of all the grids
  int ig = 0;
  for (int i=0; i<ih; i++)
  {
    ig += m_pDS->m_NTheta[i];
  }
  ig += ip;

  //calculate the vectors of bisector and normal of mole2
  vector<double>& a20 = m_atoms[m_n1]->m_x;
  vector<double>& a21 = m_atoms[m_n1 + 1]->m_x;
  vector<double>& a22 = m_atoms[m_n1 + 2]->m_x;

  vector<double> v0 = subtraction(a21, a20);
  vector<double> v1 = subtraction(a22, a20);

  //These two vectors must be unit vector
  vector<double> bisect = get_bisect_unit(v0, v1);
  vector<double> normal = get_normal_unit(v0, v1);

  //relative coords of mole2 in the cubic of the 8 corners
  r -= m_pDS->m_R_NDX[ir];
  ang1 -= m_pDS->m_PHI_angles[ih];
  ang2 -= m_pDS->m_THETA_angles[ih][ip];
  double xr = r / m_pDS->m_DR[ir];

  double tr = ang1 / (m_pDS->m_PHI_angles[ih+1] - m_pDS->m_PHI_angles[ih]);
  double pr;

  if (ip == m_pDS->m_NTheta[ih]-1)
  {
    pr = ang2 / (2 * M_PI + m_pDS->m_THETA_angles[ih][0] - m_pDS->m_THETA_angles[ih][ip]);
  }
  else
  {
    pr = ang2 / (m_pDS->m_THETA_angles[ih][ip+1] - m_pDS->m_THETA_angles[ih][ip]);
  }
  m_ir = ir;
  m_ig = ig;

  // xr: distance
  // tr: Phi angle (between xr to z axis)
  // pr: Theta angle (between xr to xz face)
  m_rel_x = {xr, tr, pr};
  m_vbis = bisect;
  m_vnrm = normal;
}

int Coordinates::get_index(double r, const vector<double>& vec)
{
  int res = 0;
  if (r < vec[0])
  {
    res = -1;
  }
  else
  {
    for (size_t i=1; i<vec.size(); i++)
    {
      if (r <= vec[i])
      {
        res = i - 1;
        break;
      }
    }
  }

  return res;

}


void Coordinates::calt_conf_energy(database::EnergeForceDatabase& allconfig, bool isForce=false, double ehigh=100.0)
{
  double ri = m_ir;

  if (ri > 100.0)
  {
    m_properties = {{"E", 0.0}};
  }
  if (ri < 0.0)
  {
    m_properties = {{"E", ehigh}};
  }

  double gi = m_ig;
  vector<double> relative_x = m_rel_x;
  vector<double> bisv = m_vbis;
  vector<double> nrmv = m_vnrm;

  double wghx, wghy;
  vector<int> grids_sub_ndx = weights_in_subsection(bisv, wghx, wghy);

  std::map<std::string, std::vector<double>> properties = {{"E", {}}};
  std::vector<std::string> propname = {"E"};

  //if (isForce)
  //{
  //  for (size_t i=0; i<m_atoms.size(); i++)
  //  {
  //    for (size_t j=0; j<3; j++)
  //    {
  //      char name[7];
  //      sprintf(name, "f%3d%3d", i, j);
  //      properties[name] = {};
  //      propname.push_back(name);
  //    }
  //  }
  //}
  for (int i=ri; i<ri+2; i++)
  {
    for (int j: m_pDS->m_grid_data[gi])
    {
      std::map<std::string, std::vector<double>> prop = {{"E", {}}};

      for (int ni: grids_sub_ndx)
      {
        std::vector<database::Atom> xconf0 = allconfig.at_all(i, j, ni, 0)->m_xmole2;
        std::vector<database::Atom> xconf1 = allconfig.at_all(i, j, ni, 1)->m_xmole2;

        std::pair<double, double> w = database::weights_for_2_configs(nrmv, xconf0, xconf1);

        for (std::string pp: propname)
        {
          double p0 = allconfig.get_prop(i, j, ni, 0, pp, w.first, ehigh);
          double p1 = allconfig.get_prop(i, j, ni, 1, pp, w.second, ehigh);

          prop[pp].push_back(p1 * fabs(w.second) + p0 * fabs(w.first));
        }
      }

      for (std::string pp: propname)
      {
        vector<double>& tmp = prop[pp];
        double psub = bilinear(tmp[0], tmp[1], tmp[2], tmp[3], wghx, wghy);

        properties[pp].push_back(psub);
      }

    }
  }
  m_properties = {};
  for (std::string pp: propname)
  {
    double pCubic[8] = {
      properties[pp][0], properties[pp][4],
      properties[pp][2], properties[pp][6],
      properties[pp][1], properties[pp][5],
      properties[pp][3], properties[pp][7],
    };

    // occording to the likely.pyx
    likely::TriCubicInterpolator grids8(pCubic, 1.0, 2, 2, 2);
    double stepx=1.0, stepy=1.0,stepz=1.0;
    m_properties[pp] = grids8.test(relative_x[0]/stepx, relative_x[1]/stepy, relative_x[2]/stepz);
  }

  if (isForce)
  {
    // TODO  Add ImMirrorForce function, as in this version it no used
    // so just ignore it
    // ImMirrorForce();
  }

}


double Coordinates::get_interp_energy()
{
  return m_properties["E"];
}

/**
 * QMInterpolation definition
 */
const std::vector<std::vector<int> > QMInterpolation::m_aa_ndx = {{1, 2, 3}};
const std::vector<int> QMInterpolation::m_c_ndx = {1};
const std::vector<int> QMInterpolation::m_H_ndx = {2,3,4};
const std::vector<int> QMInterpolation::m_prob_ndx = {4,5,6};

QMInterpolation::QMInterpolation(std::string fftype, database::EnergeForceDatabase& allconfig):
  m_fftype(fftype), m_allconfig(allconfig)
{};

// return the output string, then write to file
std::string QMInterpolation::process(std::string filename)
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

  double test = 0;
  vector<double> dist = {};
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

    interp.ReorientToOrigin();
    interp.MirrorAll();
    try
    {
       interp.indexing();
    } catch (std::invalid_argument &e)
    {
      fprintf(stderr, "%s\n", name_str);
    }

    interp.calt_conf_energy(m_allconfig);

    test += interp.get_interp_energy();
    dist.push_back(interp.m_r);
  }

  char res[100];
  sprintf(res, "%s %12.7f %.3f", filename.c_str(), test, dist[0]);

  return std::string(res);
}
