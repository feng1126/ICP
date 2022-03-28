// svd.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <vector>
#include "../eigen/Eigen/Eigen"
#include "../eigen/Eigen/Geometry"
#include "ceres/ceres.h"

#include "ceres/ceres.h"
#include "so3.h"
#include "se3.h"
#include "data.h"
#include "lm.h"
#include "local_parameterization_se3.h"


inline Eigen::Quaterniond YPR2Quaterniond(float roll, float pitch, float yaw)
{
    Eigen::Vector3d eulerAngle(yaw, pitch, roll);
    Eigen::AngleAxisd rollAngle(Eigen::AngleAxisd(eulerAngle(2), Eigen::Vector3d::UnitX()));
    Eigen::AngleAxisd pitchAngle(Eigen::AngleAxisd(eulerAngle(1), Eigen::Vector3d::UnitY()));
    Eigen::AngleAxisd yawAngle(Eigen::AngleAxisd(eulerAngle(0), Eigen::Vector3d::UnitZ()));
    Eigen::Quaterniond quaternion;
    quaternion = yawAngle * pitchAngle * rollAngle;
    return quaternion;
}




class external_cali
{
public:
    external_cali(Eigen::Vector3d point1, Eigen::Vector3d point2)
    {
        pd1 = point1;
        pd2 = point2;
    }

    template <typename T>
    bool operator()(const T* _q, const T* _t, T* residuals) const
    {
        Eigen::Quaternion<T> q_incre{ _q[3], _q[0], _q[1], _q[2] };

        Eigen::Matrix<T, 3, 1> p_l{ T(pd1[0]) , T(pd1[1]) ,T(pd1[2]) };
        Eigen::Matrix<T, 3, 1> t_incre(_t[0], _t[1], _t[2]);
        Eigen::Matrix<T, 3, 1> p_c = q_incre.toRotationMatrix() * p_l + t_incre;

        residuals[0] = T(pd2[0]) - p_c[0];
        residuals[1] = T(pd2[1]) - p_c[1];
        residuals[2] = T(pd2[2]) - p_c[2];

        return true;
    }

    static ceres::CostFunction* Create(Eigen::Vector3d point1, Eigen::Vector3d point2)
    {
        return (new ceres::AutoDiffCostFunction<external_cali, 3, 4, 3>(new external_cali(point1, point2)));
    }
private:
    Eigen::Vector3d  pd1,pd2;

};

class external_cali2 : public ceres::SizedCostFunction<3, 6>
{
public:
    external_cali2(Eigen::Vector3d point1, Eigen::Vector3d point2)
    {
        pd1 = point1;
        pd2 = point2;
    }

    virtual bool Evaluate(double const* const* parameters,
        double* residuals,
        double** jacobians) const
    {

        Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(parameters[0]);
        Sophus::SE3 T = Sophus::SE3::exp(lie); //指数映射到se3

        auto rp1 = T * pd1;


        residuals[0] = pd2[0] - rp1[0];
        residuals[1] = pd2[1] - rp1[1];
        residuals[2] = pd2[2] - rp1[2];

        Eigen::Matrix<double, 3,6> jaci = Eigen::Matrix<double, 3, 6>::Zero();
        jaci.block<3, 3>(0, 0) = -Eigen::Matrix3d::Identity();
        jaci.block<3, 3>(0, 3) = Sophus::SO3::hat(rp1);

        int k = 0;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                if (jacobians) {
                    if (jacobians[0])
                        jacobians[0][k] = jaci(i, j);//0-18
                }
                k++;
            }
        }

        return true;
    
    }


private:
    Eigen::Vector3d  pd1, pd2;

};


int main()
{

    Eigen::Matrix<double, 3, 6> Jac = Eigen::Matrix<double, 3, 6>::Zero();


    Eigen::Matrix< double, 6, 1 > Vector6d;

   // Vector6d.setZero();

    //std::cout << Vector6d[5] << std::endl;

    //std::cout << Jac << "\n";

    std::cout << "Hello World!\n";
    std::vector<Eigen::Vector3d > points1, points2;

    Eigen::Quaterniond q = YPR2Quaterniond(2.52, 2.23, 1.58);

    Eigen::Vector3d trans = Eigen::Vector3d{ 12,2.3, 2 };
    int num = 1000;

    for (int i = 0; i < num; i++)
    {
       Eigen::Vector3d temp(1 + 2.1 * i, 1.2 + 3.8 * i, 1 + 1.5 * i);
       points1.push_back(temp);
       Eigen::Vector3d temp2 = q * temp + trans;
       points2.push_back(temp2);
    }

    std::string  pathE1 = "D:\\vs_project\\dataSync\\dataSync\\E1.txt";
    std::string  pathMPU = "D:\\vs_project\\dataSync\\dataSync\\MPU.txt";
    std::vector<data > data1 = readData(pathE1);
    std::vector<data > data2 = readData(pathMPU);

    points1.clear(); points2.clear();
    for (int i = 0; i < data1.size(); i++)
    {
        points1.push_back(data1[i].xyz);

    }

    for (int i = 0; i < data2.size(); i++)
    {
        points2.push_back(data2[i].xyz);

    }


        Eigen::Vector3d srcMean1, srcMean2;
        srcMean1.setZero();
        srcMean2.setZero();
        for (int i = 0; i < num; i++)
        {

            srcMean1  = srcMean1 +  points1[i];
            srcMean2  = srcMean2 +  points2[i];
        }

        srcMean1 = srcMean1 / num;
        srcMean2 = srcMean2 / num;

        std::vector<Eigen::Vector3d > centerpoints1, centerpoints2;
        for (int i = 0; i < num; i++)
        {

            Eigen::Vector3d temp1 =  points1[i] - srcMean1;
            Eigen::Vector3d temp2 = points2[i] - srcMean2;
            centerpoints1.push_back(temp1);
            centerpoints2.push_back(temp2);
        }

        Eigen::Matrix3d  w;
        w.setZero();

        for (int i = 0; i < num; i++)
        {
            w += centerpoints1[i] * centerpoints2[i].transpose();

        }

       Eigen::JacobiSVD<Eigen::Matrix3d> svd(w, Eigen::ComputeFullU | Eigen::ComputeFullV);

        const Eigen::Matrix3d U = svd.matrixU();
        const Eigen::Matrix3d V = svd.matrixV();

        Eigen::Matrix3d Vt = V;
       Eigen::Matrix3d R  = V * U.transpose();

       //std::cout << "R.determinant " << R.determinant() << std::endl;

       if (R.determinant() < 0) 
       {
           Vt.block<1, 3>(2, 0) *= -1;
           R = Vt * U.transpose();
       }

        std::cout << "R 0" << q.toRotationMatrix() << std::endl << std::endl;
        Eigen::Vector3d t = srcMean2- R * srcMean1;

        std::cout <<" 1 " << std::endl;
        std::cout << "R " << R << std::endl << std::endl;
        std::cout << "t " <<  t  << std::endl << std::endl;

        ceres::LocalParameterization* q_parameterization = new ceres::EigenQuaternionParameterization();
        ceres::Problem problem;
        Eigen::Quaterniond q1;
        q1 = R;

        double ext[4];
        ext[0] = q1.x();
        ext[1] = q1.y();
        ext[2] = q1.z();
        ext[3] = q1.w();
        problem.AddParameterBlock(ext, 4, q_parameterization);


        double tran[3];

        tran[0] = t[0];
        tran[1] = t[1];
        tran[2] = t[2];

        for (int i=0; i< num; i++ )
        {
            ceres::CostFunction* cost_function;
            cost_function = external_cali::Create(points1[i], points2[i]);
            problem.AddResidualBlock(cost_function, NULL, ext, tran);
        }

        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.minimizer_progress_to_stdout = false;
        options.max_num_iterations = 100;
        options.min_trust_region_radius = 1e-128;

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << std::endl << std::endl;

        

        Eigen::Quaterniond qext;
        qext.w() = ext[3];
        qext.x() = ext[0];
        qext.y() = ext[1];
        qext.z() = ext[2];


        std::cout << "2 " << std::endl;
        std::cout << "R " << qext.toRotationMatrix() << std::endl;
        std::cout << "tran " << tran[0] << " " << tran[1] << " " << tran[2] << std::endl;

        {
            double pose[6];
            Sophus::SE3 T = Sophus::SE3(R, t);
            Sophus::Vector6d T3 = T.log();


            for (int j = 0; j < 6; ++j) 
            {
                pose[j] = T3(j);
            }

            ceres::Problem problem;

            ceres::LocalParameterization* local_param = new SE3Parameterization();
            problem.AddParameterBlock(pose, 6, local_param);

            for (int i = 0; i < num; i++)
            {
                ceres::CostFunction* cost_function = new external_cali2(points1[i], points2[i]);
                problem.AddResidualBlock(cost_function, NULL,
                    pose);
            }

            ceres::Solver::Options options;
            options.linear_solver_type = ceres::DENSE_QR;
            options.minimizer_progress_to_stdout = false;
            options.max_num_iterations = 100;
            options.min_trust_region_radius = 1e-128;

            ceres::Solver::Summary summary;
            ceres::Solve(options, &problem, &summary);
            std::cout << summary.BriefReport() << std::endl;

            Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(pose);
            T = Sophus::SE3::exp(lie); //指数映射到se3


            std::cout << "3 " << std::endl;
            std::cout << "R " << T.rotation_matrix() << std::endl;
            std::cout << "tran " << T.translation().transpose() << std::endl;

        }


        {
            double pose[6];
            Sophus::SE3 T = Sophus::SE3(R, t);
            Sophus::Vector6d T3 = T.log();

            for (int j = 0; j < 6; ++j)
            {
                pose[j] = T3[j];

                std::cout << "pose[j] " << j << " " << pose[j] << std::endl;
            }

            LM(points1, points2, pose);
            GN(points1, points2, pose);
        }
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
