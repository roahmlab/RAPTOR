#ifndef SPATIAL_HPP
#define SPATIAL_HPP

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace RAPTOR {
namespace Spatial {

// [Xj,S] = jcalc( jtyp, q )
template<typename Scalar>
void jcalc(Eigen::Matrix<Scalar, 6, 6>& Xj, 
           Eigen::Vector<Scalar, 6>& S, 
           const int jtyp, 
           const Scalar& q) {
    Xj.setIdentity();
    S.setZero();

    Scalar c = cos(q);
    Scalar s = sin(q);

    if (jtyp < 0) { // reversed direction
        // c = cos(-q) = cos(q)
        s = -s; // s = sin(-q) = -sin(q)
    }

    switch (jtyp)
    {
        case 0: // fixed joint
            // do nothing, results already set
            break;
        case 1: // revolute X axis 'Rx'
            Xj(1,1) = c;
            Xj(1,2) = s;
            Xj(2,1) = -s;
            Xj(2,2) = c;
            S(0) = 1;
            break;
        case 2: // revolute Y axis 'Ry'
            Xj(0,0) = c;
            Xj(0,2) = -s;
            Xj(2,0) = s;
            Xj(2,2) = c;
            S(1) = 1;
            break;
        case 3: // revolute Z axis 'Rz'
            Xj(0,0) = c;
            Xj(0,1) = s;
            Xj(1,0) = -s;
            Xj(1,1) = c;
            S(2) = 1;
            break;
        case -1: // reversed revolute X axis '-Rx'
            Xj(1,1) = c;
            Xj(1,2) = s;
            Xj(2,1) = -s;
            Xj(2,2) = c;
            S(0) = -1;
            break;
        case -2: // reversed revolute Y axis '-Ry'
            Xj(0,0) = c;
            Xj(0,2) = -s;
            Xj(2,0) = s;
            Xj(2,2) = c;
            S(1) = -1;
            break;
        case -3: // reversed revolute Z axis '-Rz'
            Xj(0,0) = c;
            Xj(0,1) = s;
            Xj(1,0) = -s;
            Xj(1,1) = c;
            S(2) = -1;
            break;
        case 4: // prismatic X axis 'Px'
            Xj(4,2) = q;
            Xj(5,1) = -q;
            S(3) = 1;
            break;
        case 5: // prismatic Y axis 'Py'
            Xj(3,2) = -q;
            Xj(5,0) = q;
            S(4) = 1;
            break;
        case 6: // prismatic Z axis 'Pz'
            Xj(3,1) = q;
            Xj(4,0) = -q;
            S(5) = 1;
            break;
        case -4: // reversed prismatic X axis '-Px'
            Xj(4,2) = -q;
            Xj(5,1) = q;
            S(3) = -1;
            break;
        case -5: // reversed prismatic Y axis '-Py'
            Xj(3,2) = q;
            Xj(5,0) = -q;
            S(4) = -1;
            break;
        case -6: // reversed prismatic Z axis '-Pz'
            Xj(4,2) = q;
            Xj(5,1) = -q;
            S(5) = -1;
            break;
        default:
            throw std::invalid_argument("spatial.hpp: jcalc(): unknown joint type!");
            break;
    }
    
    if (fabs(jtyp) <= 3) {
        Xj.block(3,3,3,3) = Xj.block(0,0,3,3);
    }
}

template<typename Scalar>
void jcalc(Eigen::Matrix<Scalar, 6, 6>& Xj, 
           Eigen::Matrix<Scalar, 6, 6>& dXjdq,
           Eigen::Vector<Scalar, 6>& S, 
           const int jtyp, 
           const Scalar& q) {
    Xj.setIdentity();
    dXjdq.setZero();
    S.setZero();

    Scalar c = cos(q);
    Scalar s = sin(q);
    Scalar dcdq = -s;
    Scalar dsdq = c;

    if (jtyp < 0) { // reversed direction
        // c = cos(-q) = cos(q)
        s = -s; // s = sin(-q) = -sin(q)
        dsdq = -dsdq;
    }

    switch (jtyp)
    {
        case 0: // fixed joint
            // do nothing, results already set
            break;
        case 1: // revolute X axis 'Rx'
            Xj(1,1) = c;
            Xj(1,2) = s;
            Xj(2,1) = -s;
            Xj(2,2) = c;
            dXjdq(1,1) = dcdq;
            dXjdq(1,2) = dsdq;
            dXjdq(2,1) = -dsdq;
            dXjdq(2,2) = dcdq;
            S(0) = 1;
            break;
        case 2: // revolute Y axis 'Ry'
            Xj(0,0) = c;
            Xj(0,2) = -s;
            Xj(2,0) = s;
            Xj(2,2) = c;
            dXjdq(0,0) = dcdq;
            dXjdq(0,2) = -dsdq;
            dXjdq(2,0) = dsdq;
            dXjdq(2,2) = dcdq;
            S(1) = 1;
            break;
        case 3: // revolute Z axis 'Rz'
            Xj(0,0) = c;
            Xj(0,1) = s;
            Xj(1,0) = -s;
            Xj(1,1) = c;
            dXjdq(0,0) = dcdq;
            dXjdq(0,1) = dsdq;
            dXjdq(1,0) = -dsdq;
            dXjdq(1,1) = dcdq;
            S(2) = 1;
            break;
        case -1: // reversed revolute X axis '-Rx'
            Xj(1,1) = c;
            Xj(1,2) = s;
            Xj(2,1) = -s;
            Xj(2,2) = c;
            dXjdq(1,1) = dcdq;
            dXjdq(1,2) = dsdq;
            dXjdq(2,1) = -dsdq;
            dXjdq(2,2) = dcdq;
            S(0) = -1;
            break;
        case -2: // reversed revolute Y axis '-Ry'
            Xj(0,0) = c;
            Xj(0,2) = -s;
            Xj(2,0) = s;
            Xj(2,2) = c;
            dXjdq(0,0) = dcdq;
            dXjdq(0,2) = -dsdq;
            dXjdq(2,0) = dsdq;
            dXjdq(2,2) = dcdq;
            S(1) = -1;
            break;
        case -3: // reversed revolute Z axis '-Rz'
            Xj(0,0) = c;
            Xj(0,1) = s;
            Xj(1,0) = -s;
            Xj(1,1) = c;
            dXjdq(0,0) = dcdq;
            dXjdq(0,1) = dsdq;
            dXjdq(1,0) = -dsdq;
            dXjdq(1,1) = dcdq;
            S(2) = -1;
            break;
        case 4: // prismatic X axis 'Px'
            Xj(4,2) = q;
            Xj(5,1) = -q;
            dXjdq(4,2) = 1;
            dXjdq(5,1) = -1;
            S(3) = 1;
            break;
        case 5: // prismatic Y axis 'Py'
            Xj(3,2) = -q;
            Xj(5,0) = q;
            dXjdq(3,2) = -1;
            dXjdq(5,0) = 1;
            S(4) = 1;
            break;
        case 6: // prismatic Z axis 'Pz'
            Xj(3,1) = q;
            Xj(4,0) = -q;
            dXjdq(3,1) = 1;
            dXjdq(4,0) = -1;                                                                                                                                                                                                                                                                                                                                                                          
            S(5) = 1;
            break;
        case -4: // reversed prismatic X axis '-Px'
            Xj(4,2) = -q;
            Xj(5,1) = q;
            dXjdq(4,2) = -1;
            dXjdq(5,1) = 1;
            S(3) = -1;
            break;
        case -5: // reversed prismatic Y axis '-Py'
            Xj(3,2) = q;
            Xj(5,0) = -q;
            dXjdq(3,2) = 1;
            dXjdq(5,0) = -1;
            S(4) = -1;
            break;
        case -6: // reversed prismatic Z axis '-Pz'
            Xj(4,2) = q;
            Xj(5,1) = -q;
            dXjdq(4,2) = 1;
            dXjdq(5,1) = -1;
            S(5) = -1;
            break;
        default:
            throw std::invalid_argument("spatial.hpp: jcalc(): unknown joint type!");
            break;
    }
    
    if (fabs(jtyp) <= 3) {
        Xj.block(3,3,3,3) = Xj.block(0,0,3,3);
        dXjdq.block(3,3,3,3) = dXjdq.block(0,0,3,3);
    }
}

template<typename Scalar>
void jcalc(Eigen::Matrix<Scalar, 6, 6>& Xj, 
           Eigen::Matrix<Scalar, 6, 6>& dXjdt,
           Eigen::Vector<Scalar, 6>& S, 
           const int jtyp, 
           const Scalar& q,
           const Scalar& q_d) {
    Xj.setIdentity();
    dXjdt.setZero();
    S.setZero();

    Scalar c = cos(q);
    Scalar s = sin(q);
    Scalar dcdt = -s * q_d;
    Scalar dsdt = c * q_d;

    if (jtyp < 0) { // reversed direction
        // c = cos(-q) = cos(q)
        s = -s; // s = sin(-q) = -sin(q)
        dsdt = -dsdt;   
    }

    switch (jtyp)
    {
        case 0: // fixed joint
            // do nothing, results already set
            break;
        case 1: // revolute X axis 'Rx'
            Xj(1,1) = c;
            Xj(1,2) = s;
            Xj(2,1) = -s;
            Xj(2,2) = c;
            dXjdt(1,1) = dcdt;
            dXjdt(1,2) = dsdt;
            dXjdt(2,1) = -dsdt;
            dXjdt(2,2) = dcdt;
            S(0) = 1;
            break;
        case 2: // revolute Y axis 'Ry'
            Xj(0,0) = c;
            Xj(0,2) = -s;
            Xj(2,0) = s;
            Xj(2,2) = c;
            dXjdt(0,0) = dcdt;
            dXjdt(0,2) = -dsdt;
            dXjdt(2,0) = dsdt;
            dXjdt(2,2) = dcdt;
            S(1) = 1;
            break;
        case 3: // revolute Z axis 'Rz'
            Xj(0,0) = c;
            Xj(0,1) = s;
            Xj(1,0) = -s;
            Xj(1,1) = c;
            dXjdt(0,0) = dcdt;
            dXjdt(0,1) = dsdt;
            dXjdt(1,0) = -dsdt;
            dXjdt(1,1) = dcdt;
            S(2) = 1;
            break;
        case -1: // reversed revolute X axis '-Rx'
            Xj(1,1) = c;
            Xj(1,2) = s;
            Xj(2,1) = -s;
            Xj(2,2) = c;
            dXjdt(1,1) = dcdt;
            dXjdt(1,2) = dsdt;
            dXjdt(2,1) = -dsdt;
            dXjdt(2,2) = dcdt;
            S(0) = -1;
            break;
        case -2: // reversed revolute Y axis '-Ry'
            Xj(0,0) = c;
            Xj(0,2) = -s;
            Xj(2,0) = s;
            Xj(2,2) = c;
            dXjdt(0,0) = dcdt;
            dXjdt(0,2) = -dsdt;
            dXjdt(2,0) = dsdt;
            dXjdt(2,2) = dcdt;
            S(1) = -1;
            break;
        case -3: // reversed revolute Z axis '-Rz'
            Xj(0,0) = c;
            Xj(0,1) = s;
            Xj(1,0) = -s;
            Xj(1,1) = c;
            dXjdt(0,0) = dcdt;
            dXjdt(0,1) = dsdt;
            dXjdt(1,0) = -dsdt;
            dXjdt(1,1) = dcdt;
            S(2) = -1;
            break;
        case 4: // prismatic X axis 'Px'
            Xj(4,2) = q;
            Xj(5,1) = -q;
            dXjdt(4,2) = q_d;
            dXjdt(5,1) = -q_d;
            S(3) = 1;
            break;
        case 5: // prismatic Y axis 'Py'
            Xj(3,2) = -q;
            Xj(5,0) = q;
            dXjdt(3,2) = -q_d;
            dXjdt(5,0) = q_d;
            S(4) = 1;
            break;
        case 6: // prismatic Z axis 'Pz'
            Xj(3,1) = q;
            Xj(4,0) = -q;
            dXjdt(3,1) = q_d;
            dXjdt(4,0) = -q_d; 
            S(5) = 1;
            break;
        case -4: // reversed prismatic X axis '-Px'
            Xj(4,2) = -q;
            Xj(5,1) = q;
            dXjdt(4,2) = -q_d;
            dXjdt(5,1) = q_d;
            S(3) = -1;
            break;
        case -5: // reversed prismatic Y axis '-Py'
            Xj(3,2) = q;
            Xj(5,0) = -q;
            dXjdt(3,2) = q_d;
            dXjdt(5,0) = -q_d;
            S(4) = -1;
            break;
        case -6: // reversed prismatic Z axis '-Pz'
            Xj(4,2) = q;
            Xj(5,1) = -q;
            dXjdt(4,2) = q_d;
            dXjdt(5,1) = -q_d;
            S(5) = -1;
            break;
        default:
            throw std::invalid_argument("spatial.hpp: jcalc(): unknown joint type!");
            break;
    }
    
    if (fabs(jtyp) <= 3) {
        Xj.block(3,3,3,3) = Xj.block(0,0,3,3);
        dXjdt.block(3,3,3,3) = dXjdt.block(0,0,3,3);
    }
}

template<typename Scalar>
Eigen::Matrix<Scalar, 6, 6> crm(const Eigen::Vector<Scalar, 6>& v) {
    Eigen::Matrix<Scalar, 6, 6> vcross;
    vcross.setZero();

    // vcross = [  0    -v(3)  v(2)   0     0     0    ;
    //             v(3)  0    -v(1)   0     0     0    ;
    //             -v(2)  v(1)  0      0     0     0    ;
    //             0    -v(6)  v(5)   0    -v(3)  v(2) ;
    //             v(6)  0    -v(4)   v(3)  0    -v(1) ;
    //             -v(5)  v(4)  0     -v(2)  v(1)  0 ];

    vcross(0, 1) = -v(2);
    vcross(0, 2) = v(1);
    vcross(1, 0) = v(2);
    vcross(1, 2) = -v(0);
    vcross(2, 0) = -v(1);
    vcross(2, 1) = v(0);
    vcross.block(3,3,3,3) = vcross.block(0,0,3,3);

    vcross(3, 1) = -v(5);
    vcross(3, 2) = v(4);
    vcross(4, 0) = v(5);
    vcross(4, 2) = -v(3);
    vcross(5, 0) = -v(4);
    vcross(5, 1) = v(3);

    return vcross;
}

template<typename Scalar>
Eigen::Matrix<Scalar, 6, 6> crf(const Eigen::Vector<Scalar, 6>& v) {
    Eigen::Matrix<Scalar, 6, 6> vcross = crm(v);
    return -vcross.transpose();
}

}; // namespace Spatial
}; // namespace RAPTOR

#endif // SPATIAL_HPP
