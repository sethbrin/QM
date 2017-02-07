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
#include <array>

#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>

#include "common.h"
#include "grid_structures.h"
#include "read_energy_force_new.h"

namespace QM
{
  class Atom
  {
  public:
    std::vector<double> m_x;
    std::array<double, 3> m_f;

    Atom(): m_f({0,0,0}) {}
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

    void MirrorAll();

    void indexing();
    void calt_conf_energy(database::EnergeForceDatabase& allconfig, bool isForce, double ehigh);

    double get_interp_energy();

    int m_ir;
    int m_ig;
    std::vector<double> m_vbis;
    std::vector<double> m_vnrm;
    std::vector<double>  m_rel_x;

    std::map<std::string, double> m_properties;
    int m_n1;
    int m_n2;
    std::string m_fragtype;
    std::string m_name;
    database::DataStructure* m_pDS;
    bool m_is_oriented;
    std::map<std::string, int> m_facendx;
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

    // get index of a struct
    int get_index(double r, const vector<double>& vec);
  };

  class QMInterpolation
  {
  public:
    explicit QMInterpolation(std::string fftype, database::EnergeForceDatabase& allconfig);

    std::string process(std::string filename);

  private:
    std::string m_fftype;
    database::EnergeForceDatabase& m_allconfig;

    const static std::vector<std::vector<int> > m_aa_ndx;
    const static std::vector<int> m_c_ndx;
    const static std::vector<int> m_prob_ndx;
    const static std::vector<int> m_H_ndx;
  };
}
#endif
