#ifndef SOPHUS_LOCAL_PARAMETERIZATION_SE3_HPP
#define SOPHUS_LOCAL_PARAMETERIZATION_SE3_HPP

#include <ceres/local_parameterization.h>
#include <ceres/internal/eigen.h>
#include "se3.h"

class SE3Parameterization :public ceres::LocalParameterization
{
public:
    SE3Parameterization() {}
    virtual ~SE3Parameterization() {}
    virtual bool Plus(const double* x,
        const double* delta,
        double* x_plus_delta)const
    {
        Sophus::Vector6d lie(x);
        Sophus::Vector6d lie_delta(delta);
        Sophus::SE3 T = Sophus::SE3::exp(lie);
        Sophus::SE3 T_delta = Sophus::SE3::exp(lie_delta);

        Eigen::Map<Sophus::Vector6d> x_plus(x_plus_delta);
        x_plus = (T_delta * T).log();
        return true;

    }
    //流形到其切平面的雅克比矩阵
    virtual bool ComputeJacobian(const double* x,
        double* jacobian) const
    {
        ceres::MatrixRef(jacobian, 6, 6) = ceres::Matrix::Identity(6, 6);//ceres::MatrixRef()函数也类似于Eigen::Map,通过模板参数传入C++指针将c++数据映射为Ceres::Matrix
        return true;
    }
    //定义流形和切平面维度：在本问题中是李代数到李代数故都为6
    virtual int GlobalSize()const { return 6; }
    virtual int LocalSize()const { return 6; }
};
#endif // SOPHUS_LOCAL_PARAMETERIZATION_SE3_HPP
