#pragma once
#include "data.h"
#include "so3.h"
#include "se3.h"



void calib()
{
	Eigen::Matrix4d ql, qr;
    Eigen::Matrix4d R1, R2;


    int frameSize = 10;

    Eigen::MatrixXd A(frameSize * 4, 4); //构建公式（7）中的A矩阵


    double qw = Eigen::Quaterniond(R1).w();
    Eigen::Vector3d qv = Eigen::Quaterniond(R1).vec();
    ql.block<3, 3>(0, 0) = qw * Eigen::Matrix3d::Identity() + skewSymmetric(qv);
    ql.block<3, 1>(0, 3) = qv;
    ql.block<1, 3>(3, 0) = -qv.transpose();

    Eigen::Quaterniond(R1).angularDistance

    qw = Eigen::Quaterniond(R2).w();
    qv = Eigen::Quaterniond(R2).vec();
    qr.block<3, 3>(0, 0) = qw * Eigen::Matrix3d::Identity() - skewSymmetric(qv);
    qr.block<3, 1>(0, 3) = qv;
    qr.block<1, 3>(3, 0) = -qv.transpose();

   // A.block<4, 4>((i - 1) * 4, 0) = huber * (ql - qr);

    


}