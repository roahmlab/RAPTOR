#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

namespace RAPTOR {
namespace Utils {

// Converts degrees to radians
// @param deg: Angle in degrees
// @return: Angle in radians
inline double deg2rad(const double deg) {
    return deg * M_PI / 180.0;
}

// Converts radians to degrees
// @param rad: Angle in radians
// @return: Angle in degrees
inline double rad2deg(const double rad) {
    return rad * 180.0 / M_PI;
}

// Converts a vector of angles in degrees to radians
// @param deg: Vector of angles in degrees
// @return: Vector of angles in radians
inline Eigen::Vector3d deg2rad(const Eigen::Vector3d& deg) {
    return deg * M_PI / 180.0;
}

// Converts a vector of angles in radians to degrees
// @param rad: Vector of angles in radians
// @return: Vector of angles in degrees
inline Eigen::Vector3d rad2deg(const Eigen::Vector3d& rad) {
    return rad * 180.0 / M_PI;
}

// Converts a vector of angles in degrees to radians
// @param deg: Vector of angles in degrees
// @return: Vector of angles in radians
inline Eigen::VectorXd deg2rad(const Eigen::VectorXd& deg) {
    return deg * M_PI / 180.0;
}

// Converts a vector of angles in radians to degrees
// @param rad: Vector of angles in radians
// @return: Vector of angles in degrees
inline Eigen::VectorXd rad2deg(const Eigen::VectorXd& rad) {
    return rad * 180.0 / M_PI;
}

// Wraps an angle to the range [-pi, pi]
// @param angle: Angle in radians
// @return: Wrapped angle in radians
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

// Wraps a vector of angles to the range [-pi, pi]
// @param angles: Vector of angles in radians
// @return: Vector of wrapped angles in radians
inline Eigen::VectorXd wrapToPi(const Eigen::VectorXd& angles) {
    Eigen::VectorXd res = angles;
    for (int i = 0; i < res.size(); i++) {
        res(i) = wrapToPi(res(i));
    }
    return res;
}

// Returns the sign of a value
// @param val: Input value
// @param eps: Tolerance for zero
// @return: 1.0 if val > eps, -1.0 if val < -eps, 0.0 otherwise
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

// Checks if two vectors are equal within a tolerance
// @param a: First vector
// @param b: Second vector
// @param tol: Tolerance for equality
// @return: True if vectors are equal, false otherwise
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

// Checks if two matrices are equal within a tolerance
// @param a: First matrix
// @param b: Second matrix
// @param tol: Tolerance for equality
// @return: True if matrices are equal, false otherwise
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

// Reshapes a vector into a matrix
// @param vec: Input vector
// @param rows: Number of rows in the output matrix
// @param cols: Number of columns in the output matrix
// @return: Reshaped matrix
inline Eigen::MatrixXd reshape(const Eigen::VectorXd& vec, 
                               int rows, 
                               int cols) {
    return Eigen::Map<const Eigen::MatrixXd>(vec.data(), rows, cols);
}

// Computes the skew-symmetric matrix of a vector
// @param v: Input vector
// @return: Skew-symmetric matrix
inline Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d res;
    res << 0,    -v(2), v(1),
           v(2),  0,   -v(0),
          -v(1),  v(0), 0;
    return res;
}

// Computes the vector from a skew-symmetric matrix
// @param m: Input skew-symmetric matrix
// @return: Vector
inline Eigen::Vector3d skew(const Eigen::Matrix3d& m) {
    Eigen::Vector3d res;
    res << m(2,1) - m(1,2), 
           m(0,2) - m(2,0), 
           m(1,0) - m(0,1);
    return 0.5 * res;
}

// Computes the vector from a skew-symmetric matrix
// @param m: Input skew-symmetric matrix
// @return: Vector
inline Eigen::Vector3d unskew(const Eigen::Matrix3d& m) {
    Eigen::Vector3d res;
    res << m(2,1), 
           m(0,2), 
           m(1,0);
    return res;
}

// Computes the plucker transformation matrix
// @param R: Rotation matrix
// @param p: Translation vector
// @return: Plucker transformation matrix
inline Eigen::Matrix<double, 6, 6> plux(const Eigen::Matrix3d& R, 
                                       const Eigen::Vector3d& p) {
    Eigen::Matrix<double, 6, 6> res;
    res << R,            Eigen::MatrixXd::Zero(3, 3),
           -R * skew(p), R;
    return res;
}

// Initializes an Eigen vector from an array
// @param array: Input array
// @param size: Size of the array
// @return: Eigen vector
inline Eigen::VectorXd initializeEigenVectorFromArray(const double* array, 
                                                      size_t size) {
    Eigen::VectorXd res(size);
    for (int i = 0; i < size; i++) {
        res(i) = array[i];
    }
    return res;
}

// Initializes an Eigen matrix from a file
// @param filename: Name of the file
// @return: Eigen matrix
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

    if (data.size() == 0) {
        throw std::runtime_error("Empty file: " + filename);
    }

    Eigen::MatrixXd res(data.size(), data[0].size());
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            res(i, j) = data[i][j];
        }
    }

    return res;
}

// Writes an Eigen matrix to a file
// @param matrix: Input matrix
// @param filename: Name of the file
// @param precision: Precision of the output
inline void writeEigenMatrixToFile(const Eigen::MatrixXd& matrix, 
                                   const std::string filename,
                                   const size_t precision = 6) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file " + filename);
    }

    file << std::setprecision(precision);

    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.cols(); j++) {
            file << matrix(i, j) << " ";
        }
        file << std::endl;
    }

    file.close();
}

// Uniformly samples a vector
// @param vec: Input vector
// @param numSamples: Number of samples
// @return: Sampled vector
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

// Uniformly samples a matrix in rows
// @param mat: Input matrix
// @param samples: Output sampled matrix
// @param numSamples: Number of samples
inline void uniformlySampleMatrixInRows(const Eigen::MatrixXd& mat, 
                                        Eigen::MatrixXd& samples,
                                        int numSamples) {
    if (numSamples <= 0 || 
        mat.rows() == 0 || 
        mat.cols() == 0) {
        throw std::invalid_argument("Invalid input arguments");
    }

    if (numSamples >= mat.rows()) {
        throw std::invalid_argument("Number of samples should be less than the number of rows in the matrix");
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
}

// Uniformly samples a matrix in columns
// @param mat: Input matrix
// @param samples: Output sampled matrix
// @param numSamples: Number of samples
inline void uniformlySampleMatrixInCols(const Eigen::MatrixXd& mat, 
                                        Eigen::MatrixXd& samples,
                                        int numSamples) {
    if (numSamples <= 0 || 
        mat.rows() == 0 || 
        mat.cols() == 0) {
        throw std::invalid_argument("Invalid input arguments");
    }

    if (numSamples >= mat.cols()) {
        throw std::invalid_argument("Number of samples should be less than the number of columns in the matrix");
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
}

}; // namespace Utils
}; // namespace RAPTOR

#endif // UTILS_H
