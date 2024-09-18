#include "SysIDAlgFull.h"
#include "FirstOptimization.h"
#include "SecondOptimization.h"
#include <iostream>

SysIDAlgFull::SysIDAlgFull(
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
)
    : na_idx_(na_idx), nu_idx_(nu_idx),
      lb_(lb), ub_(ub), lba_(lba), uba_(uba),
      X0_1_(X0_1), T_(T), data_(data), dataFull_(dataFull),
      W_ip_(W_ip), Ginv_(Ginv), Aid_(Aid), Ad_(Ad), Kd_(Kd)
{
    // 初始化变量
    n_ = static_cast<int>(variables[0]);
    m_ = static_cast<int>(variables[1]);
    b_ = static_cast<int>(variables[2]);
    d_ = static_cast<int>(variables[3]);

    // 算法选项
    k_ = AlgOptions[0];
    tol_ = AlgOptions[1];
    MS1_ = static_cast<int>(AlgOptions[2]);
    MS2_ = static_cast<int>(AlgOptions[3]);
    regroup_ = static_cast<bool>(AlgOptions[4]);
    SearchAlpha_ = static_cast<bool>(AlgOptions[7]);
    includeOffset_ = static_cast<bool>(AlgOptions[8]);
    includeConstraints_ = static_cast<bool>(AlgOptions[10]);
    constraintVariant_ = static_cast<int>(AlgOptions[11]);

    initialize();
}

void SysIDAlgFull::initialize() {
    // 初始化参数维度
    p_ip_ = 10 * n_;

    if (includeOffset_) {
        p_full_ = p_ip_ + 4 * n_;
        b_full_ = b_ + 4 * n_;
    } else {
        p_full_ = p_ip_ + 3 * n_;
        b_full_ = b_ + 3 * n_;
    }


    

    // 其他初始化操作
}

void SysIDAlgFull::runAlgorithm() {
    // 处理约束
    processConstraints();

    // 参数重组
    if (regroup_) {
        regroupParameters();
    }

    // 运行第一个优化问题
    runFirstOptimization();

    // 根据需要运行第二个优化问题
    if (SearchAlpha_) {
        runSecondOptimization();
    }

    // 最终结果保存在 X_, Wfull_, alphanew_
}

void SysIDAlgFull::processConstraints() {
    // 根据 includeConstraints_ 和 constraintVariant_ 处理约束
    // 具体实现省略
}

void SysIDAlgFull::regroupParameters() {
    // 参数重组的实现
    // 具体实现省略
}

void SysIDAlgFull::runFirstOptimization() {
    // 构建第一个优化问题的对象
    FirstOptimization firstOpt(/* 传入必要的参数 */);

    // 设置初始点
    firstOpt.setStartingPoint(X0_1_);

    // 运行优化
    firstOpt.solve();

    // 获取结果
    X_ = firstOpt.getSolution();
}

void SysIDAlgFull::runSecondOptimization() {
    // 构建第二个优化问题的对象
    SecondOptimization secondOpt(/* 传入必要的参数 */);

    // 设置初始点
    Eigen::VectorXd X0_2; // 根据需要初始化
    secondOpt.setStartingPoint(X0_2);

    // 运行优化
    secondOpt.solve();

    // 获取结果
    alphanew_ = secondOpt.getAlphaNew();
    // 更新 X_
    // 具体实现省略
}

Eigen::VectorXd SysIDAlgFull::getX() const {
    return X_;
}

Eigen::MatrixXd SysIDAlgFull::getWfull() const {
    return Wfull_;
}

Eigen::VectorXd SysIDAlgFull::getAlphaNew() const {
    return alphanew_;
}
