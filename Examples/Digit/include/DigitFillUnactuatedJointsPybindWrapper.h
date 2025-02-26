#ifndef KINOVA_PYBIND_WRAPPER_H
#define KINOVA_PYBIND_WRAPPER_H

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include "DigitDynamicsConstraints.h"

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"

namespace RAPTOR {
namespace Digit {

namespace nb = nanobind;

class DigitFillUnactuatedJointsPybindWrapper {
public:
    using Model = pinocchio::Model;
    using Vec3 = Eigen::Vector3d;
    using VecX = Eigen::VectorXd;
    using MatX = Eigen::MatrixXd;

    using nb_1d_double = nb::ndarray<double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;
    using nb_2d_double = nb::ndarray<double, nb::ndim<2>, nb::c_contig, nb::device::cpu>;

    // Constructor
    DigitFillUnactuatedJointsPybindWrapper() = default;

    DigitFillUnactuatedJointsPybindWrapper(const char stanceLeg) {
        const std::string urdf_filename = "Robots/digit-v3/digit-v3-armfixedspecific-floatingbase-springfixed.urdf";
        
        modelPtr_ = std::make_shared<Model>();
        pinocchio::urdf::buildModel(urdf_filename, *modelPtr_);

        Transform stance_foot_T_des(3, -M_PI_2);

        dcPtr_ = std::make_shared<DigitDynamicsConstraints>(modelPtr_, 
                                                            stanceLeg,
                                                            stance_foot_T_des);
    }

    // Destructor
    ~DigitFillUnactuatedJointsPybindWrapper() = default;

    // Class methods
    nb::tuple fill_q(const nb_2d_double& qas) {
        if (qas.shape(0) != NUM_INDEPENDENT_JOINTS) {
            throw std::invalid_argument("Input qa has wrong number of rows");
        }

        qs.resize(modelPtr_->nq, qas.shape(1));
        
        for (int i = 0; i < qas.shape(1); i++) {
            VecX q(modelPtr_->nq);
            VecX qa(NUM_INDEPENDENT_JOINTS);

            for (int j = 0; j < qas.shape(0); j++) {
                qa(j) = qas(j, i);
            }
            dcPtr_->fill_independent_vector(q, qa, true);
            dcPtr_->setupJointPosition(q);

            qs.col(i) = q;
        }

        const size_t shape_ptr[] = {modelPtr_->nq, qas.shape(1)};
        auto result = nb::ndarray<nb::numpy, const double>(qs.data(),
                                                           2,
                                                           shape_ptr,
                                                           nb::handle());

        return nb::make_tuple(result);
    }

    nb::tuple fill_qv(const nb_2d_double& qas, 
                      const nb_2d_double& qa_ds) {
        if (qas.shape(0) != NUM_INDEPENDENT_JOINTS) {
            throw std::invalid_argument("Input qa has wrong number of rows");
        }
        if (qa_ds.shape(0) != NUM_INDEPENDENT_JOINTS) {
            throw std::invalid_argument("Input qa_d has wrong number of rows");
        }
        if (qa_ds.shape(1) != qas.shape(1)) {
            throw std::invalid_argument("Input qa_d has wrong number of columns");
        }

        qs.resize(modelPtr_->nq, qas.shape(1));
        q_ds.resize(modelPtr_->nv, qas.shape(1));

        for (int i = 0; i < qas.shape(1); i++) {
            VecX q(modelPtr_->nq);
            VecX q_d(modelPtr_->nv);
            VecX q_dd(modelPtr_->nv);
            VecX qa(NUM_INDEPENDENT_JOINTS);
            VecX qa_d(NUM_INDEPENDENT_JOINTS);

            for (int j = 0; j < qas.shape(0); j++) {
                qa(j) = qas(j, i);
                qa_d(j) = qa_ds(j, i);
            }
            dcPtr_->fill_independent_vector(q, qa, true);
            dcPtr_->fill_independent_vector(q_d, qa_d, true);
            dcPtr_->setupJointPositionVelocityAcceleration(q, q_d, q_dd);

            qs.col(i) = q;
            q_ds.col(i) = q_d;
        }

        const size_t shape_ptr[] = {modelPtr_->nq, qas.shape(1)};
        auto result1 = nb::ndarray<nb::numpy, const double>(qs.data(),
                                                            2,
                                                            shape_ptr,
                                                            nb::handle());
        auto result2 = nb::ndarray<nb::numpy, const double>(q_ds.data(),
                                                            2,
                                                            shape_ptr,
                                                            nb::handle());
        return nb::make_tuple(result1, result2);
    }

    nb::tuple fill_qva(const nb_2d_double& qas, 
                       const nb_2d_double& qa_ds, 
                       const nb_2d_double& qa_dds) {
        if (qas.shape(0) != NUM_INDEPENDENT_JOINTS) {
            throw std::invalid_argument("Input qa has wrong number of rows");
        }
        if (qa_ds.shape(0) != NUM_INDEPENDENT_JOINTS) {
            throw std::invalid_argument("Input qa_d has wrong number of rows");
        }
        if (qa_dds.shape(0) != NUM_INDEPENDENT_JOINTS) {
            throw std::invalid_argument("Input qa_dd has wrong number of rows");
        }
        if (qa_ds.shape(1) != qas.shape(1)) {
            throw std::invalid_argument("Input qa_d has wrong number of columns");
        }
        if (qa_dds.shape(1) != qas.shape(1)) {
            throw std::invalid_argument("Input qa_dd has wrong number of columns");
        }
        
        qs.resize(modelPtr_->nq, qas.shape(1));
        q_ds.resize(modelPtr_->nv, qas.shape(1));
        q_dds.resize(modelPtr_->nv, qas.shape(1));

        for (int i = 0; i < qas.shape(1); i++) {
            VecX q(modelPtr_->nq);
            VecX q_d(modelPtr_->nv);
            VecX q_dd(modelPtr_->nv);
            VecX qa(NUM_INDEPENDENT_JOINTS);
            VecX qa_d(NUM_INDEPENDENT_JOINTS);
            VecX qa_dd(NUM_INDEPENDENT_JOINTS);

            for (int j = 0; j < qas.shape(0); j++) {
                qa(j) = qas(j, i);
                qa_d(j) = qa_ds(j, i);
                qa_dd(j) = qa_dds(j, i);
            }
            dcPtr_->fill_independent_vector(q, qa, true);
            dcPtr_->fill_independent_vector(q_d, qa_d, true);
            dcPtr_->fill_independent_vector(q_dd, qa_dd, true);
            dcPtr_->setupJointPositionVelocityAcceleration(q, q_d, q_dd);

            qs.col(i) = q;
            q_ds.col(i) = q_d;
            q_dds.col(i) = q_dd;
        }

        const size_t shape_ptr[] = {modelPtr_->nq, qas.shape(1)};
        auto result1 = nb::ndarray<nb::numpy, const double>(qs.data(),
                                                            2,
                                                            shape_ptr,
                                                            nb::handle());
        auto result2 = nb::ndarray<nb::numpy, const double>(q_ds.data(),
                                                            2,
                                                            shape_ptr,
                                                            nb::handle());
        auto result3 = nb::ndarray<nb::numpy, const double>(q_dds.data(),
                                                            2,
                                                            shape_ptr,
                                                            nb::handle());
        return nb::make_tuple(result1, result2, result3);
    }

    // Class members
    std::shared_ptr<Model> modelPtr_;
    std::shared_ptr<DigitDynamicsConstraints> dcPtr_;

    // results
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> qs;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> q_ds;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> q_dds;
};

}; // namespace Kinova
}; // namespace RAPTOR

#endif // KINOVA_PYBIND_WRAPPER_H
