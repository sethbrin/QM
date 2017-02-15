#include "qm_interpolation.h"
#include "grid_structures.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <fstream>
#include <exception>
#include <set>

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

InpAtom::InpAtom(const std::vector<double>& points, std::string name) {
    m_name = name;
    m_x = points;
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

void Coordinates::addAtom(std::vector<double>& points, std::string name, const std::string& ftype) {
    if (ftype == "gms")
    {
        m_atoms.push_back(boost::shared_ptr<Atom>(new InpAtom(points, name)));
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
    m_operateNdx.push_back(0);
    m_operations.push_back({dvec, 0});

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
    m_operateNdx.push_back(1);
    m_operations.push_back({ax, ang});

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
        m_operateNdx.push_back(2);
        m_operations.push_back({ax, ang});
    }

    m_is_oriented = true;

    spherical_x();
}

void Coordinates::ReorientToOldVec() {
    auto item = m_operations[index(2, m_operateNdx)];
    m_force = rotate(m_force, item.first, item.second * -1);
    m_torque = rotate(m_torque, item.first, item.second * -1);

    item = m_operations[index(1, m_operateNdx)];
    m_force = rotate(m_force, item.first, item.second * -1);
    m_torque = rotate(m_torque, item.first, item.second * -1);
}

//Calculate the coords in spherical coordination system for molecule 2
void Coordinates::spherical_x()
{
    double totalM = 0;
    std::vector<double> x = {0, 0, 0};
    for (int i = m_n1; i < m_atoms.size(); i++) {
        for (int k = 0; k < 3; k++) {
            x[k] += m_atoms[i]->m_x[k] * TMASS[i - m_n1];
        }
        totalM += TMASS[i - m_n1];
    }
    for (int k = 0; k < 3; k++) {
        x[k] /= totalM;
    }
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
    m_origin_com = m_center2;
    for (std::string face: m_pDS->m_symface)
    {
        int fndx = m_facendx[face];

        if (m_center2[fndx] < 0.0)
        {
            m_symm[fndx] = -1;

            for (size_t i=m_n1; i<m_atoms.size(); i++)
            {
                m_atoms[i]->m_x[fndx] *= -1;
            }
        }
    }

    spherical_x();
}

void Coordinates::MirrorBackProperty()
{

    for (std::string face: m_pDS->m_symface)
    {
        int fndx = m_facendx[face];

        if (m_origin_com[fndx] < 0.0)
        {
            m_symm[fndx] = -1;
            m_force[fndx] *= -1;
            for (size_t i=0; i<3; i++) {
                if (i != fndx) {
                    m_torque[i] *= -1;
                }
            }
        }
    }

}


void Coordinates::spherical_orient() {
    //int r = m_orientVec.size();
    // 原来python版本是直接length一个np.array对象，返回1，
    // 觉得有问题，但此处只这样写
    int r = 1;
    const std::vector<double>& x = m_orientVec;

    double ang1 = M_PI * 0.5 - std::acos(x[2]/r);
    double ang2 = 0;

    if (std::abs(x[0]) < 0.000001) {
        if (x[1] > 0) ang2 = M_PI * 0.5;
        else ang2 = M_PI * 1.5;
    } else {
        ang2 = std::atan(x[1]/x[0]);
        if (x[0] < 0) ang2 += M_PI;
        else if (x[1] < 0) ang2 += M_PI * 2;
    }
    m_orient_ang1 = ang1;
    m_orient_ang2 = ang2;
}

void Coordinates::indexing_orient_auto3(int ri)
{
    double ang1 = m_orient_ang1;
    double ang2 = m_orient_ang2;

    database::OrientDataStructure* OrientDS = m_orient_DS[ri];
    // ndx of ang1 (Phi):
    int ih = get_index(ang1, OrientDS->m_PHI_angles);
    int ang1_ndx1 = ih;
    int ang1_ndx2 = ih + 1;
    int ang1_ndx3 = 0;
    if (ang1_ndx1 == OrientDS->m_PHI_angles.size() - 2) {
        ang1_ndx3 = ih - 1;
    } else if (ang1_ndx1 == 0) {
        ang1_ndx3 = ih + 2;
    } else {
        double tmp1 = OrientDS->m_PHI_angles[ih+2] - ang1;
        double tmp2 = ang1 - OrientDS->m_PHI_angles[ih-1];
        if (tmp1 < tmp2) ang1_ndx3 = ih + 2;
        else ang1_ndx3 = ih - 1;
    }

    bool iflinear = false;
    std::set<int> phiList = {ang1_ndx1, ang1_ndx2, ang1_ndx3};
    std::map<int, std::vector<int>> dgrid_sub_ndx = {};
    std::map<int, std::vector<int>> dtheta_ndx = {};

    for (int kk : phiList) {
        dgrid_sub_ndx[kk] = {};
        dtheta_ndx[kk] = {};

        // ndx_of_ang2 (Theta)
        int ip = -1;
        for (size_t i=1; i<(size_t)OrientDS->m_NTheta[kk]; i++)
        {
            if (ang2 <= OrientDS->m_THETA_angles[kk][i])
            {
                ip = i - 1;
                break;
            }
        }
        if (ip == -1)
        {
            ip = OrientDS->m_NTheta[kk] - 1;
        }
        //ig is the index of all the grids
        int ig = 0;
        for (int i=0; i<kk; i++)
        {
            ig += OrientDS->m_NTheta[i];
        }
        ig += ip;

        dgrid_sub_ndx[kk].push_back(ig);
        dtheta_ndx[kk].push_back(ip);

        if (ip == OrientDS->m_NTheta[kk] - 1) {
            if (OrientDS->m_NTheta[kk] == 1) {
                dgrid_sub_ndx[kk].push_back(ig);
                dtheta_ndx[kk].push_back(0);

                if (!iflinear) {
                    dgrid_sub_ndx[kk].push_back(ig);
                    dtheta_ndx[kk].push_back(0);
                }
            } else {
                dgrid_sub_ndx[kk].push_back(ig - OrientDS->m_NTheta[kk]+1);
                dtheta_ndx[kk].push_back(0+OrientDS->m_NTheta[kk]);

                if (!iflinear) {
                    double tmp1 = OrientDS->m_THETA_angles[kk][1] - ang2 + 2 * M_PI;
                    double tmp2 = ang2 - OrientDS->m_THETA_angles[kk][ip-1];

                    if (tmp1 < tmp2) {
                        dgrid_sub_ndx[kk].push_back(ig - OrientDS->m_NTheta[kk]+1+1);
                        dtheta_ndx[kk].push_back(0+OrientDS->m_NTheta[kk]+1);
                    } else {
                        dgrid_sub_ndx[kk].push_back(ig-1);
                        dtheta_ndx[kk].push_back(ip-1);
                    }

                }

            }
        } else {
            dgrid_sub_ndx[kk].push_back(ig+1);
            dtheta_ndx[kk].push_back(ip+1);

            if (!iflinear) {
                double tmp1, tmp2;
                if (ip == OrientDS->m_NTheta[kk] - 2) {
                    tmp1 = 2 * M_PI - ang2;
                } else {
                    tmp1 = OrientDS->m_THETA_angles[kk][ip+2] - ang2;
                }

                if (ip == 0) {
                    tmp2 = ang2 - OrientDS->m_THETA_angles[kk][OrientDS->m_NTheta[kk]-1] +  2*M_PI;
                } else {
                    tmp2 = ang2 - OrientDS->m_THETA_angles[kk][ip-1];
                }

                if (tmp1 < tmp2) {
                    if (ip == OrientDS->m_NTheta[kk] - 2) {
                        dgrid_sub_ndx[kk].push_back(ig + 1- OrientDS->m_NTheta[kk]+1);
                        dtheta_ndx[kk].push_back(0+OrientDS->m_NTheta[kk]);
                    } else {
                        dgrid_sub_ndx[kk].push_back(ig+2);
                        dtheta_ndx[kk].push_back(ip+2);
                    }
                } else {
                    if (ip == 0) {
                        dgrid_sub_ndx[kk].push_back(ig + OrientDS->m_NTheta[kk]-1);
                        dtheta_ndx[kk].push_back(-1);
                    } else {
                        dgrid_sub_ndx[kk].push_back(ig-1);
                        dtheta_ndx[kk].push_back(ip-1);
                    }
                }
            }
        }

    }

    insert_map_element_2_map(dgrid_sub_ndx, ri, m_dgrid_sub_ndx);
    insert_map_element_2_map(dtheta_ndx, ri, m_dtheta_ndx);

}


void Coordinates::indexing_auto3()
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
        m_r_ndxs = {ir};
        m_vbis = {0, 0, 0};
        m_vnrm = {0 ,0, 0};
        return;
    }

    std::vector<int> r_ndxs = {ir, ir+1};

    // find 3 layers which are close the query one
    if (ir == 0) r_ndxs.push_back(ir + 2);
    else if (ir == m_pDS->m_R_NDX.size() - 2) r_ndxs.push_back(ir - 1);
    else {
        if (std::abs(r - m_pDS->m_R_NDX[ir-1]) < m_pDS->m_R_NDX[ir+2] - r) {
            r_ndxs.push_back(ir - 1);
        } else {
            r_ndxs.push_back(ir + 2);
        }
    }

    // ndx of ang1 (Phi):
    int ih = get_index(ang1, m_pDS->m_PHI_angles);
    int ang1_ndx1 = ih;
    int ang1_ndx2 = ih + 1;
    int ang1_ndx3 = 0;
    if (ang1_ndx1 == m_pDS->m_PHI_angles.size() - 2) {
        ang1_ndx3 = ih - 1;
    } else if (ang1_ndx1 == 0) {
        ang1_ndx3 = ih + 2;
    } else {
        double tmp1 = m_pDS->m_PHI_angles[ih+2] - ang1;
        double tmp2 = ang1 - m_pDS->m_PHI_angles[ih-1];
        if (tmp1 < tmp2) ang1_ndx3 = ih + 2;
        else ang1_ndx3 = ih - 1;
    }

    bool iflinear = false;
    std::set<int> phiList = {ang1_ndx1, ang1_ndx2, ang1_ndx3};
    m_dgrid_ndx_layer = {};
    m_dtheta_ndx_layer = {};

    for (int kk : phiList) {
        m_dgrid_ndx_layer[kk] = {};
        m_dtheta_ndx_layer[kk] = {};

        // ndx_of_ang2 (Theta)
        int ip = -1;
        for (size_t i=1; i<(size_t)m_pDS->m_NTheta[kk]; i++)
        {
            if (ang2 <= m_pDS->m_THETA_angles[kk][i])
            {
                ip = i - 1;
                break;
            }
        }
        if (ip == -1)
        {
            ip = m_pDS->m_NTheta[kk] - 1;
        }
        //ig is the index of all the grids
        int ig = 0;
        for (int i=0; i<kk; i++)
        {
            ig += m_pDS->m_NTheta[i];
        }
        ig += ip;

        m_dgrid_ndx_layer[kk].push_back(ig);
        m_dtheta_ndx_layer[kk].push_back(ip);

        if (ip == m_pDS->m_NTheta[kk] - 1) {
            if (m_pDS->m_NTheta[kk] == 1) {
                m_dgrid_ndx_layer[kk].push_back(ig);
                m_dtheta_ndx_layer[kk].push_back(0);

                if (!iflinear) {
                    m_dgrid_ndx_layer[kk].push_back(ig);
                    m_dtheta_ndx_layer[kk].push_back(0);
                }
            } else {
                m_dgrid_ndx_layer[kk].push_back(ig - 1);
                m_dtheta_ndx_layer[kk].push_back(ip-1);

                if (!iflinear) {
                    m_dgrid_ndx_layer[kk].push_back(ig-2);
                    m_dtheta_ndx_layer[kk].push_back(ip-2);
                }

            }
        } else {
            m_dgrid_ndx_layer[kk].push_back(ig+1);
            m_dtheta_ndx_layer[kk].push_back(ip+1);

            if (!iflinear) {
                if (ip == m_pDS->m_NTheta[kk] - 2) {
                    m_dgrid_ndx_layer[kk].push_back(ig-1);
                    m_dtheta_ndx_layer[kk].push_back(ip-1);
                } else if (ip == 0) {
                    m_dgrid_ndx_layer[kk].push_back(ig+2);
                    m_dtheta_ndx_layer[kk].push_back(ip+2);
                } else {
                    double tmp1 = m_pDS->m_THETA_angles[kk][ip+2] - ang2;
                    double tmp2 = ang2 - m_pDS->m_THETA_angles[kk][ip-1];

                    if (tmp1 < tmp2) {
                        m_dgrid_ndx_layer[kk].push_back(ig+2);
                        m_dtheta_ndx_layer[kk].push_back(ip+2);
                    } else {
                        m_dgrid_ndx_layer[kk].push_back(ig-1);
                        m_dtheta_ndx_layer[kk].push_back(ip-1);
                    }
                }
            }
        }

    }


    //calculate the vectors of bisector and normal of mole2
    vector<double>& a20 = m_atoms[m_n1]->m_x;
    vector<double>& a21 = m_atoms[m_n1 + 1]->m_x;
    vector<double>& a22 = m_atoms[m_n1 + 2]->m_x;

    vector<double> v0 = subtraction(a21, a20);
    vector<double> v1 = subtraction(a22, a20);

    //These two vectors must be unit vector
    vector<double> bisect = get_bisect_unit(v0, v1);
    vector<double> normal = get_normal_unit(v0, v1);

    // xr: distance
    // tr: Phi angle (between xr to z axis)
    // pr: Theta angle (between xr to xz face)
    m_r_ndxs = r_ndxs;
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
    for (int ri : m_r_ndxs) {
        if (ri > 100.0)
        {
            m_properties = {{"E", 0.0}};
            return;
        }
        if (ri < 0.0)
        {
            m_properties = {{"E", ehigh}};
            return;
        }
    }

    vector<double> bisv = m_vbis;
    vector<double> nrmv = m_vnrm;
    std::vector<double> grid_ndx_layer = {};
    for (const auto& item : m_dgrid_ndx_layer) {
        grid_ndx_layer.insert(grid_ndx_layer.end(), item.second.begin(), item.second.end());
    }
    m_orientVec = bisv;
    spherical_orient();

    double ang1 = m_orient_ang1;
    double ang2 = m_orient_ang2;
    ang2 = ang2 * R2D + 180;
    ang2 -= (static_cast<int>(ang2) / 360) * 360;
    ang2 = ang2 / R2D;
    m_orient_ang2 = ang2;

    m_orient_DS = {};
    m_dgrid_sub_ndx = {};
    m_dtheta_ndx = {};

    std::map<int, double> orient_pr;
    std::map<int, std::vector<int>> grids_sub_ndx;
    std::map<int, double> wghx1 = {};
    std::map<int, double> wghx2 = {};
    std::map<int, double> wghy = {};
    std::map<int, int> label = {};

    for (int i : m_r_ndxs) {
        double dist = m_pDS->m_R_NDX[i];

        if (dist > 5.5000001) {
            std::vector<int> grids_sub_ndx_tmp;
            double wghx_tmp, wghy_tmp;
            grids_sub_ndx_tmp = weights_in_subsection(bisv, wghx_tmp, wghy_tmp);
            grids_sub_ndx[i] = grids_sub_ndx_tmp;
            wghx1[i] = wghx_tmp / PI4;
            wghx2[i] = wghx_tmp / PI4;
            wghy[i] = wghy_tmp / PI4;
            label[i] = 0;
        } else {
            database::OrientDataStructure* OrientDS = NULL;
            if (dist < 2.5000001) OrientDS = database::getOrientDataStructure2("wtr");
            else if(dist > 2.5000001 && dist < 3.5000001) OrientDS = database::getOrientDataStructure3("wtr");
            else OrientDS = database::getOrientDataStructure2("wtr");

            m_orient_DS[i] = OrientDS;

            indexing_orient_auto3(i);

            if (m_dtheta_ndx[i].size() == 2) {
                // current not support
            } else if (m_dtheta_ndx[i].size() == 3) {
                for (const auto& item : m_dtheta_ndx[i]) {
                    grids_sub_ndx[i].insert(grids_sub_ndx[i].end(),
                            m_dgrid_sub_ndx[i][item.first].begin(),
                            m_dgrid_sub_ndx[i][item.first].end());
                }
                label[i] = 2;
            }
        }
    }

    std::map<std::string, std::vector<double>> properties = {{"E", {}}, {"Fx", {}}, {"Fy", {}}, {"Fz", {}}, {"Tx", {}}, {"Ty", {}}, {"Tz", {}}};
    std::vector<std::string> propname = {"E", "Fx", "Fy", "Fz", "Tx", "Ty", "Tz"};
    std::map<std::string, std::vector<double>> tmp_prop = {{"E", {}}, {"Fx", {}}, {"Fy", {}}, {"Fz", {}}, {"Tx", {}}, {"Ty", {}}, {"Tz", {}}};

    for (int i : m_r_ndxs) {
        for (auto j : grid_ndx_layer) {
            std::map<std::string, std::vector<double>> prop = tmp_prop;
            for (auto ni : grids_sub_ndx[i]) {
                std::vector<std::vector<double>> xvecs;
                for (int ff = 0; ff < m_pDS->m_nNorm[i]; ff++) {
                    std::vector<database::Atom> xconf = allconfig.at_all(i, j, ni, ff)->m_xmole2;
                    xvecs.push_back(norm_prob(xconf, {0, 1, 2}, "wtr"));
                }

                int nvec = xvecs.size();

                if (nvec == 2) {
                    // linear interpolation for normal vectors
                    double w0, w1;
                    int ndx0, ndx1;
                    weights_for_normal_general(nrmv, xvecs, w0, w1, ndx0, ndx1);

                    for (std::string pp: propname) {
                        double p0 = allconfig.get_prop(i, j, ni, ndx0, pp, w0, ehigh);
                        double p1 = allconfig.get_prop(i, j, ni, ndx1, pp, w1, ehigh);

                        double p = p1 * std::abs(w1) + p0 * std::abs(w0) ;
                        prop[pp].push_back(p);
                    }

                } else if (nvec > 2) {
                    auto  res = get_neighors_for_normal(nrmv, xvecs);
                    double angNorm = res.first;
                    int ndx1 = res.second[0];
                    int ndx2 = res.second[1];
                    int ndx3 = res.second[2];


                    double angNorm_1 = ndx1 * M_PI / nvec;
                    double angNorm_2 = ndx2 * M_PI / nvec;
                    double angNorm_3 = ndx3 * M_PI / nvec;

                    for (std::string pp: propname) {
                        if (ndx1 == nvec) ndx1 = 0;
                        if (ndx2 == nvec) ndx2 = 0;
                        if (ndx3 == nvec) ndx3 = 0;
                        double p1 = allconfig.get_prop(i, j, ni, ndx1, pp, 0, ehigh);
                        double p2 = allconfig.get_prop(i, j, ni, ndx2, pp, 0, ehigh);
                        double p3 = allconfig.get_prop(i, j, ni, ndx3, pp, 0, ehigh);

                        double p = lagrange_interp({{angNorm_1, p1}, {angNorm_2, p2}, {angNorm_3, p3}}, angNorm);
                        prop[pp].push_back(p);
                    }

                }
            }

            for (std::string pp: propname) {
                // on the level of orientation, theta and phi
                if (prop[pp].size() == 4) {
                    double psub = bilinear_gen(prop[pp][0], prop[pp][1], prop[pp][2], prop[pp][3], wghx1[i], wghx2[i], wghy[i], label[i]);
                    properties[pp].push_back(psub);
                } else if (prop[pp].size() == 9) {
                    int cn = 0;
                    std::vector<std::pair<double, double>> points_phi;

                    for (auto& item : m_dtheta_ndx[i]) {
                        int kk = item.first;
                        double angPhi = m_orient_DS[i]->m_PHI_angles[kk];

                        std::set<int> tmp_set(m_dtheta_ndx[i][kk].begin(), m_dtheta_ndx[i][kk].end());
                        if (tmp_set.size() == 1) {
                            points_phi.push_back({angPhi, prop[pp][cn]});
                            cn += 3;
                            continue;
                        }

                        std::vector<std::pair<double, double>> points_theta;
                        for (int ip : m_dtheta_ndx[i][kk]) {
                            double angTheta;
                            if (ip >= m_orient_DS[i]->m_NTheta[kk]) {
                                angTheta = 2 * M_PI + m_orient_DS[i]->m_THETA_angles[kk][ip-m_orient_DS[i]->m_NTheta[kk]];
                            } else if (ip < 0) {
                                ip = m_orient_DS[i]->m_NTheta[kk] + ip;
                                angTheta = m_orient_DS[i]->m_THETA_angles[kk][ip] - 2 * M_PI;
                            } else {
                                angTheta = m_orient_DS[i]->m_THETA_angles[kk][ip];
                            }

                            points_theta.push_back({angTheta, prop[pp][cn]});
                            cn += 1;
                        }
                        points_phi.push_back({angPhi, lagrange_interp(points_theta, ang2)});
                    }

                    properties[pp].push_back(lagrange_interp(points_phi, ang1));
                }
            }


        }
    }

    m_properties = {};
    if (m_dtheta_ndx_layer.size() == 2) {

    } else if(m_dtheta_ndx_layer.size() == 3) {
        for (std::string pp: propname) {
            std::vector<double> psub_r;

            for (int m = 0; m < properties[pp].size(); m += 9) {
                int count = 0;
                std::vector<std::pair<double, double>> points_th;

                for (auto& item : m_dtheta_ndx_layer) {
                    int kk = item.first;
                    std::set<int> tmp_set(m_dtheta_ndx_layer[kk].begin(), m_dtheta_ndx_layer[kk].end());
                    if (tmp_set.size() == 1) {
                        points_th.push_back({m_pDS->m_PHI_angles[kk], properties[pp][m+count]});
                        count += 3;
                        continue;
                    }

                    int ip1 = m_dtheta_ndx_layer[kk][0];
                    int ip2 = m_dtheta_ndx_layer[kk][1];
                    int ip3 = m_dtheta_ndx_layer[kk][2];

                    double th1 = m_pDS->m_THETA_angles[kk][ip1];
                    double th2 = m_pDS->m_THETA_angles[kk][ip2];
                    double th3 = m_pDS->m_THETA_angles[kk][ip3];

                    double p = lagrange_interp({{th1, properties[pp][m+count]}, {th2, properties[pp][m+count+1]}, {th3, properties[pp][m+count+2]}}, m_ang2);
                    points_th.push_back({m_pDS->m_PHI_angles[kk], p});
                    count += 3;
                }
                double p = lagrange_interp(points_th, m_ang1);
                psub_r.push_back(p);

            }

            if (psub_r.size() != 3) ERROR("size of psub_r not equal to 3");

            std::vector<std::pair<double, double>> points;
            for (int t = 0; t < m_r_ndxs.size(); t++) {
                points.push_back({m_pDS->m_R_NDX[m_r_ndxs[t]], psub_r[t]});
            }
            m_properties[pp] = lagrange_interp(points, m_r);
        }
    }
}


double Coordinates::get_interp_energy()
{
    return m_properties["E"];
}

void Coordinates::reverse_force_toque() {
    m_force = {m_properties["Fx"],
        m_properties["Fy"], m_properties["Fz"]};
    m_torque = {m_properties["Tx"], m_properties["Ty"], m_properties["Tz"]};

    MirrorBackProperty();
    ReorientToOldVec();
}

std::vector<double> Coordinates::get_interp_force() {
    return m_force;
}

std::vector<double> Coordinates::get_interp_torque() {
    return m_torque;
}


/**
 * QMInterpolation definition
 */
const std::vector<int> QMInterpolation::m_aa_ndx = {1, 2, 3};
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
        if (line == "" || boost::starts_with(line, "$end") ||
                boost::starts_with(line, "$END"))
        {
            break;
        }

        lines.push_back(line);
    }
    ifs.close();

    char name_str[1024];
    sprintf(name_str, "%s%02d", filename.c_str(), m_aa_ndx[0]);
    std::string name(name_str);
    Coordinates interp = Coordinates(3,3, m_fftype, name);

    for (auto ind : m_aa_ndx)
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
        interp.indexing_auto3();
    } catch (std::invalid_argument &e)
    {
        fprintf(stderr, "%s\n", name_str);
    }

    interp.calt_conf_energy(m_allconfig);

    double e_interp = interp.get_interp_energy();
    double dist = interp.m_r;
    interp.reverse_force_toque();


    char res[100];
    //sprintf(res, "%s %12.7f %.3f", filename.c_str(), e_interp, dist);
    sprintf(res, "%s %12.7f %12.7f %12.7f", filename.c_str(), interp.get_interp_force()[0], interp.get_interp_force()[1], interp.get_interp_force()[2]);

    return std::string(res);
}

void QMInterpolation::calculate(const std::map<std::string, std::vector<double>>& lhs,
        const std::map<std::string, std::vector<double>>& rhs) {
    std::string name;
    Coordinates interp = Coordinates(3,3, m_fftype, name);

    for (auto item : lhs)
    {
        interp.addAtom(item.second, item.first, "gms");
    }

    for (auto item : rhs)
    {
        interp.addAtom(item.second, item.first, "gms");
    }

    interp.ReorientToOrigin();
    interp.MirrorAll();
    try
    {
        interp.indexing_auto3();
    } catch (std::invalid_argument &e)
    {
        fprintf(stderr, "indexing_auto3 execption\n");
    }

    interp.calt_conf_energy(m_allconfig);

    m_energy = interp.get_interp_energy();
    interp.reverse_force_toque();
    m_force = interp.get_interp_force();
    m_torque = interp.get_interp_torque();
}
