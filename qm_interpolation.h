/*
 * The QM interpolation
 *
 * Author: Ping Zeng(zengping10212044@gmail.com)
 */
#ifndef _QM_INTERPOLATION_H_
#define _QM_INTERPOLATION_H_

#include <vector>
#include <string>
#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>

#include "common.h"
#include "grid_structures.h"

namespace QM
{
  class Atom
  {
  public:
    std::vector<double> m_x;
    virtual ~Atom(){}

    virtual Atom* copy() { return NULL;}
  };

  // TODO complete the struct
  class GjfAtom : public Atom
  {
  public:
    explicit GjfAtom(std::string line);
    ~GjfAtom(){}

    Atom* copy();

    std::string m_name;
  };

  // TODO complete the struct
  class PdbAtom : public Atom
  {
  public:
    explicit PdbAtom(std::string line);
    ~PdbAtom(){}

    Atom* copy();

    int m_i_atm;
    std::string m_a_name;
    std::string m_a_res;
    std::string m_a_chn;
    int m_i_res;
  };

  class InpAtom : public Atom
  {
  public:
    explicit InpAtom(std::string line);

    ~InpAtom(){}
    Atom* copy();

    std::string m_name;
  };

  class Coordinates
  {
  public:
    explicit Coordinates(int n1, int n2, std::string& fragtype, std::string& name);

    void addAtom(std::string& line, const std::string& ftype);

    void ReorientToOrigin(double cut);

    int m_n1;
    int m_n2;
    std::string m_fragtype;
    std::string m_name;
    database::DataStructure* m_pDS;
    bool m_is_oriented;
    std::map<const char*, int> m_facendx;
    std::vector<int> m_symm;
    int m_center;
    double m_r;
    double m_ang1;
    double m_ang2;
    std::vector<double> m_center2;
    std::vector<boost::shared_ptr<Atom> > m_atoms;
    std::vector<int> m_operateNdx;
    std::vector<std::vector<double>> m_operations;

  private:
    void spherical_x();
  };

  class QMInterpolation
  {
  public:
    explicit QMInterpolation(std::string fftype);

    void process(std::string filename);

  private:
    std::string m_fftype;

    const static std::vector<std::vector<int> > m_aa_ndx;
    const static std::vector<int> m_c_ndx;
    const static std::vector<int> m_prob_ndx;
    const static std::vector<int> m_H_ndx;
  };
}
#endif
