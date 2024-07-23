#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>

namespace RAPTOR {
namespace Utils {

inline double deg2rad(const double deg) {
    return deg * M_PI / 180.0;
}

inline double rad2deg(const double rad) {
    return rad * 180.0 / M_PI;
}

inline Eigen::Vector3d deg2rad(const Eigen::Vector3d& deg) {
    return deg * M_PI / 180.0;
}

inline Eigen::Vector3d rad2deg(const Eigen::Vector3d& rad) {
    return rad * 180.0 / M_PI;
}

inline Eigen::VectorXd deg2rad(const Eigen::VectorXd& deg) {
    return deg * M_PI / 180.0;
}

inline Eigen::VectorXd rad2deg(const Eigen::VectorXd& rad) {
    return rad * 180.0 / M_PI;
}

inline double wrapToPi(const double angle) {
    double res = angle;
    while (res > M_PI) {
        res -= 2.0 * M_PI;
    }
    while (res < -M_PI) {
        res += 2.0 * M_PI;
    }
    return res;
}

inline Eigen::VectorXd wrapToPi(const Eigen::VectorXd& angles) {
    Eigen::VectorXd res = angles;
    for (int i = 0; i < res.size(); i++) {
        res(i) = wrapToPi(res(i));
    }
    return res;
}

inline double sign(double val, double eps = 1e-8) {
    if (val > eps) {
        return 1.0;
    } 
    else if (val < -eps) {
        return -1.0;
    } 
    else {
        return 0.0;
    }
}

inline bool ifTwoVectorEqual(const Eigen::VectorXd& a, 
                             const Eigen::VectorXd& b, 
                             double tol = 1e-10) {
    if (a.size() != b.size()) {
        return false;
    }
    for (int i = 0; i < a.size(); i++) {
        if (fabs(a(i) - b(i)) > tol) {
            return false;
        }
    }
    return true;
}

inline bool ifTwoMatrixEqual(const Eigen::MatrixXd& a, 
                             const Eigen::MatrixXd& b, 
                             double tol = 1e-10) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        return false;
    }
    for (int i = 0; i < a.rows(); i++) {
        for (int j = 0; j < a.cols(); j++) {
            if (fabs(a(i, j) - b(i, j)) > tol) {
                return false;
            }
        }
    }
    return true;
}

inline Eigen::MatrixXd reshape(const Eigen::VectorXd& vec, 
                               int rows, 
                               int cols) {
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), rows, cols);
}

inline Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d res;
    res << 0,    -v(2), v(1),
           v(2),  0,   -v(0),
          -v(1),  v(0), 0;
    return res;
}

inline Eigen::Vector3d skew(const Eigen::Matrix3d& m) {
    Eigen::Vector3d res;
    res << m(2,1) - m(1,2), 
           m(0,2) - m(2,0), 
           m(1,0) - m(0,1);
    return 0.5 * res;
}

inline Eigen::Vector3d unskew(const Eigen::Matrix3d& m) {
    Eigen::Vector3d res;
    res << m(2,1), 
           m(0,2), 
           m(1,0);
    return res;
}

inline Eigen::Matrix<double, 6, 6> plux(const Eigen::Matrix3d& R, 
                                        const Eigen::Vector3d& p) {
    Eigen::Matrix<double, 6, 6> res;
    res << R,            Eigen::MatrixXd::Zero(3, 3),
           -R * skew(p), R;
    return res;
}

inline Eigen::VectorXd initializeEigenVectorFromArray(const double* array, 
                                                      size_t size) {
    Eigen::VectorXd res(size);
    for (int i = 0; i < size; i++) {
        res(i) = array[i];
    }
    return res;
}

inline Eigen::MatrixXd initializeEigenMatrixFromFile(const std::string filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file " + filename);
    }

    std::vector<std::vector<double>> data;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> lineData;
        double value;
        while (iss >> value) {
            lineData.push_back(value);
        }
        data.push_back(lineData);
    }

    Eigen::MatrixXd res(data.size(), data[0].size());
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            res(i, j) = data[i][j];
        }
    }

    return res;
}

inline void writeEigenMatrixToFile(const Eigen::MatrixXd& matrix, 
                                   const std::string filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file " + filename);
    }

    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.cols(); j++) {
            file << matrix(i, j) << " ";
        }
        file << std::endl;
    }

    file.close();
}

inline Eigen::VectorXd uniformlySampleVector(const Eigen::VectorXd& vec, 
                                             int numSamples) {
    Eigen::VectorXd samples(0);

    if (numSamples <= 0 || 
        vec.size() == 0) {
        return samples; 
    }

    if (numSamples >= vec.size()) {
        samples = vec;
        return samples;
    }

    samples.resize(numSamples);

    // Calculate sampling interval
    double interval = static_cast<double>(vec.size()) / numSamples;
    
    for (int i = 0; i < numSamples; i++) {
        int idx = static_cast<int>(i * interval);
        // Ensure the index is within the bounds of the vector
        if (idx >= 0 && idx < vec.size()) {
            samples(i) = vec(idx);
        }
    }

    return samples;
}

inline Eigen::MatrixXd uniformlySampleMatrixInRows(const Eigen::MatrixXd& mat, 
                                                   int numSamples) {
    Eigen::MatrixXd samples(0, 0);

    if (numSamples <= 0 || 
        mat.rows() == 0 || 
        mat.cols() == 0) {
        return samples; 
    }

    if (numSamples >= mat.rows()) {
        samples = mat;
        return samples;
    }

    samples.resize(numSamples, mat.cols());

    // Calculate sampling interval
    double interval = static_cast<double>(mat.rows()) / numSamples;
    
    for (int i = 0; i < numSamples; i++) {
        int idx = static_cast<int>(i * interval);
        if (idx >= 0 && idx < mat.rows()) {
            samples.row(i) = mat.row(idx);
        }
    }

    return samples;
}

inline Eigen::MatrixXd uniformlySampleMatrixInCols(const Eigen::MatrixXd& mat, 
                                                   int numSamples) {
    Eigen::MatrixXd samples(0, 0);

    if (numSamples <= 0 || 
        mat.rows() == 0 || 
        mat.cols() == 0) {
        return samples; 
    }

    if (numSamples >= mat.cols()) {
        samples = mat;
        return samples;
    }

    samples.resize(mat.rows(), numSamples);

    // Calculate sampling interval
    double interval = static_cast<double>(mat.cols()) / numSamples;
    
    for (int i = 0; i < numSamples; i++) {
        int idx = static_cast<int>(i * interval);
        if (idx >= 0 && idx < mat.cols()) {
            samples.col(i) = mat.col(idx);
        }
    }

    return samples;
}

}; // namespace Utils
}; // namespace RAPTOR

#endif // UTILS_H
