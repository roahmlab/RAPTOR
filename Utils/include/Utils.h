#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>

namespace RAPTOR {
namespace Utils {

inline float deg2rad(const float deg) {
    return deg * M_PI / 180.0;
}

inline float rad2deg(const float rad) {
    return rad * 180.0 / M_PI;
}

inline Eigen::Vector3f deg2rad(const Eigen::Vector3f& deg) {
    return deg * M_PI / 180.0;
}

inline Eigen::Vector3f rad2deg(const Eigen::Vector3f& rad) {
    return rad * 180.0 / M_PI;
}

inline Eigen::VectorXf deg2rad(const Eigen::VectorXf& deg) {
    return deg * M_PI / 180.0;
}

inline Eigen::VectorXf rad2deg(const Eigen::VectorXf& rad) {
    return rad * 180.0 / M_PI;
}

inline float wrapToPi(const float angle) {
    float res = angle;
    while (res > M_PI) {
        res -= 2.0 * M_PI;
    }
    while (res < -M_PI) {
        res += 2.0 * M_PI;
    }
    return res;
}

inline Eigen::VectorXf wrapToPi(const Eigen::VectorXf& angles) {
    Eigen::VectorXf res = angles;
    for (int i = 0; i < res.size(); i++) {
        res(i) = wrapToPi(res(i));
    }
    return res;
}

inline float sign(float val, float eps = 1e-8) {
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

inline bool ifTwoVectorEqual(const Eigen::VectorXf& a, 
                             const Eigen::VectorXf& b, 
                             float tol = 1e-10) {
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

inline bool ifTwoMatrixEqual(const Eigen::MatrixXf& a, 
                             const Eigen::MatrixXf& b, 
                             float tol = 1e-10) {
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

inline Eigen::MatrixXf reshape(const Eigen::VectorXf& vec, 
                               int rows, 
                               int cols) {
    return Eigen::Map<const Eigen::MatrixXf>(vec.data(), rows, cols);
}

inline Eigen::Matrix3f skew(const Eigen::Vector3f& v) {
    Eigen::Matrix3f res;
    res << 0,    -v(2), v(1),
           v(2),  0,   -v(0),
          -v(1),  v(0), 0;
    return res;
}

inline Eigen::Vector3f skew(const Eigen::Matrix3f& m) {
    Eigen::Vector3f res;
    res << m(2,1) - m(1,2), 
           m(0,2) - m(2,0), 
           m(1,0) - m(0,1);
    return 0.5 * res;
}

inline Eigen::Vector3f unskew(const Eigen::Matrix3f& m) {
    Eigen::Vector3f res;
    res << m(2,1), 
           m(0,2), 
           m(1,0);
    return res;
}

inline Eigen::Matrix<float, 6, 6> plux(const Eigen::Matrix3f& R, 
                                       const Eigen::Vector3f& p) {
    Eigen::Matrix<float, 6, 6> res;
    res << R,            Eigen::MatrixXf::Zero(3, 3),
           -R * skew(p), R;
    return res;
}

inline Eigen::VectorXf initializeEigenVectorFromArray(const float* array, 
                                                      size_t size) {
    Eigen::VectorXf res(size);
    for (int i = 0; i < size; i++) {
        res(i) = array[i];
    }
    return res;
}

inline Eigen::MatrixXf initializeEigenMatrixFromFile(const std::string filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file " + filename);
    }

    std::vector<std::vector<float>> data;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<float> lineData;
        float value;
        while (iss >> value) {
            lineData.push_back(value);
        }
        data.push_back(lineData);
    }

    Eigen::MatrixXf res(data.size(), data[0].size());
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            res(i, j) = data[i][j];
        }
    }

    return res;
}

inline void writeEigenMatrixToFile(const Eigen::MatrixXf& matrix, 
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

inline Eigen::VectorXf uniformlySampleVector(const Eigen::VectorXf& vec, 
                                             int numSamples) {
    Eigen::VectorXf samples(0);

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
    float interval = static_cast<float>(vec.size()) / numSamples;
    
    for (int i = 0; i < numSamples; i++) {
        int idx = static_cast<int>(i * interval);
        // Ensure the index is within the bounds of the vector
        if (idx >= 0 && idx < vec.size()) {
            samples(i) = vec(idx);
        }
    }

    return samples;
}

inline Eigen::MatrixXf uniformlySampleMatrixInRows(const Eigen::MatrixXf& mat, 
                                                   int numSamples) {
    Eigen::MatrixXf samples(0, 0);

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
    float interval = static_cast<float>(mat.rows()) / numSamples;
    
    for (int i = 0; i < numSamples; i++) {
        int idx = static_cast<int>(i * interval);
        if (idx >= 0 && idx < mat.rows()) {
            samples.row(i) = mat.row(idx);
        }
    }

    return samples;
}

inline Eigen::MatrixXf uniformlySampleMatrixInCols(const Eigen::MatrixXf& mat, 
                                                   int numSamples) {
    Eigen::MatrixXf samples(0, 0);

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
    float interval = static_cast<float>(mat.cols()) / numSamples;
    
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
