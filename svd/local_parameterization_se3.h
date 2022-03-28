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
    //���ε�����ƽ����ſ˱Ⱦ���
    virtual bool ComputeJacobian(const double* x,
        double* jacobian) const
    {
        ceres::MatrixRef(jacobian, 6, 6) = ceres::Matrix::Identity(6, 6);//ceres::MatrixRef()����Ҳ������Eigen::Map,ͨ��ģ���������C++ָ�뽫c++����ӳ��ΪCeres::Matrix
        return true;
    }
    //�������κ���ƽ��ά�ȣ��ڱ����������������������ʶ�Ϊ6
    virtual int GlobalSize()const { return 6; }
    virtual int LocalSize()const { return 6; }
};
#endif // SOPHUS_LOCAL_PARAMETERIZATION_SE3_HPP
