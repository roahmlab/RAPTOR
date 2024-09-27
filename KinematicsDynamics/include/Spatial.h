#ifndef SPATIAL_HPP
#define SPATIAL_HPP

#include "pinocchio/parsers/urdf.hpp"

#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/kinematics-derivatives.hpp"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace RAPTOR {

typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
 
// [Xj,S] = jcalc( jtyp, q )
void jcalc(Matrix6d& Xj, 
           Vector6d& S, 
           const int jtyp, 
           const double q);

void jcalc(Matrix6d& Xj, 
           Matrix6d& dXjdq,
           Vector6d& S, 
           const int jtyp, 
           const double q);

void jcalc(Matrix6d& Xj, 
           Matrix6d& dXjdt,
           Vector6d& S, 
           const int jtyp, 
           const double q,
           const double q_d);

Matrix6d crm(const Vector6d& v);

Matrix6d crf(const Vector6d& v);

}; // namespace RAPTOR

#endif // SPATIAL_HPP
