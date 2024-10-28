#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include "Utils.h"  
using namespace RAPTOR;

int main() {

    std::string basePath = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/";
    const int numFiles = 4;  
    std::vector<Eigen::MatrixXd> posDataList;
    std::vector<Eigen::MatrixXd> velDataList;
    std::vector<Eigen::MatrixXd> accDataList;
    std::vector<Eigen::MatrixXd> torqueDataList;


    for (int i = 1; i <= numFiles; ++i) {

        std::string posFile = basePath + "q_downsampled_"+ std::to_string(i) + ".csv";
        std::string velFile = basePath + "q_d_downsampled_"+ std::to_string(i) + ".csv";
        std::string accFile = basePath + "q_dd_downsampled_"+ std::to_string(i) + ".csv";
        std::string torqueFile = basePath + "tau_downsampled_"+ std::to_string(i) + ".csv";


        Eigen::MatrixXd posData = Utils::initializeEigenMatrixFromFile(posFile);
        Eigen::MatrixXd velData = Utils::initializeEigenMatrixFromFile(velFile);
        Eigen::MatrixXd accData = Utils::initializeEigenMatrixFromFile(accFile);
        Eigen::MatrixXd torqueData = Utils::initializeEigenMatrixFromFile(torqueFile);

        posDataList.push_back(posData);  
        velDataList.push_back(velData);
        accDataList.push_back(accData);
        torqueDataList.push_back(torqueData);

    }

    
    for (size_t i = 1; i < posDataList.size(); ++i) {
        if (posDataList[i].cols() != posDataList[0].cols() ||
            velDataList[i].cols() != velDataList[0].cols() ||
            accDataList[i].cols() != accDataList[0].cols() ||
            torqueDataList[i].cols() != torqueDataList[0].cols()) 
        {
            std::cerr << "Error: The number of cols in file " << i + 1 << " does not match the first file." << std::endl;
            return 1;
        }
    }

    int totalRows = 0;
    for (const auto& posData : posDataList) {
        totalRows += posData.rows();
    }

    Eigen::MatrixXd combinedPosData(totalRows, posDataList[0].cols());
    Eigen::MatrixXd combinedAccData(totalRows, posDataList[0].cols());
    Eigen::MatrixXd combinedVelData(totalRows, posDataList[0].cols());
    Eigen::MatrixXd combinedTorqueData(totalRows, posDataList[0].cols());

    int colOffset_1 = 0;
    for (const auto& posData : posDataList) {
        combinedPosData.block(colOffset_1, 0, posData.rows(), posData.cols()) = posData;
        colOffset_1 += posData.rows();
    }

    int colOffset_2 = 0;
    for (const auto& velData : velDataList) {
        combinedVelData.block(colOffset_2, 0,  velData.rows(), velData.cols()) = velData;
        colOffset_2 += velData.rows();
    }
    int colOffset_3 = 0;
    for (const auto& accData : accDataList) {
        combinedAccData.block(colOffset_3, 0,  accData.rows(), accData.cols()) = accData;
        colOffset_3 += accData.rows();
    }

    int colOffset_4 = 0;
    for (const auto& torqueData : torqueDataList) {
        combinedTorqueData.block(colOffset_4, 0, torqueData.rows(), torqueData.cols()) = torqueData;
        colOffset_4 += torqueData.rows();
    }

    const std::string posoutputFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_downsampled_6.csv";
    Utils::writeEigenMatrixToFile(combinedPosData, posoutputFile);

    const std::string veloutputFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_d_downsampled_6.csv";
    Utils::writeEigenMatrixToFile(combinedVelData, veloutputFile);

    const std::string accoutputFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/q_dd_downsampled_6.csv";
    Utils::writeEigenMatrixToFile(combinedAccData, accoutputFile);

    const std::string torqueoutputFile = "../Examples/Kinova/SystemIdentification/IterativeSystemIdentification/full_params_data/tau_downsampled_6.csv";
    Utils::writeEigenMatrixToFile(combinedTorqueData, torqueoutputFile);
 

    return 0;
}
