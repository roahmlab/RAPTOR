#ifndef SYSIDALGFULL_H
#define SYSIDALGFULL_H

#include <vector>
#include <Eigen/Dense>
#include "Utils.h"

class SysIDAlgFull {
public:
    SysIDAlgFull(
        const std::vector<double>& variables,
        const std::vector<double>& AlgOptions,
        const std::vector<int>& na_idx,
        const std::vector<int>& nu_idx,
        const Eigen::VectorXd& lb,
        const Eigen::VectorXd& ub,
        const Eigen::VectorXd& lba,
        const Eigen::VectorXd& uba,
        const Eigen::VectorXd& X0_1,
        const Eigen::VectorXd& T,
        const Eigen::MatrixXd& data,
        const Eigen::MatrixXd& dataFull,
        const Eigen::MatrixXd& W_ip,
        const Eigen::MatrixXd& Ginv,
        const Eigen::MatrixXd& Aid,
        const Eigen::MatrixXd& Ad,
        const Eigen::MatrixXd& Kd
    );

    void runAlgorithm();

    Eigen::VectorXd getX() const;
    Eigen::MatrixXd getWfull() const;
    Eigen::VectorXd getAlphaNew() const;

private:
    // 输入参数
    int n_;
    int m_;
    int b_;
    int d_;
    std::vector<int> na_idx_;
    std::vector<int> nu_idx_;
    Eigen::VectorXd lb_;
    Eigen::VectorXd ub_;
    Eigen::VectorXd lba_;
    Eigen::VectorXd uba_;
    Eigen::VectorXd X0_1_;
    Eigen::VectorXd T_;
    Eigen::MatrixXd data_;
    Eigen::MatrixXd dataFull_;
    Eigen::MatrixXd W_ip_;
    Eigen::MatrixXd Ginv_;
    Eigen::MatrixXd Aid_;
    Eigen::MatrixXd Ad_;
    Eigen::MatrixXd Kd_;

    // 算法选项
    double k_;
    double tol_;
    int MS1_;
    int MS2_;
    bool regroup_;
    bool SearchAlpha_;
    bool includeOffset_;
    bool includeConstraints_;
    int constraintVariant_;

    // 中间变量
    int p_ip_;
    int p_full_;
    int b_full_;
    Eigen::VectorXd X_;
    Eigen::MatrixXd Wfull_;
    Eigen::VectorXd alphanew_;

    // 私有方法
    void initialize();
    void processConstraints();
    void regroupParameters();
    void runFirstOptimization();
    void runSecondOptimization();
    // 其他辅助方法
};

#endif // SYSIDALGFULL_H