#pragma once
#include "data.h"
#include "so3.h"
#include "se3.h"






void GD(const std::vector<Eigen::Vector3d >& xSet, const std::vector<Eigen::Vector3d>& ySet, double const* parameters)
{
	int iterations = 1700;
	double cost, lastCost = DBL_MAX;

	double aph = 0.000001;

	Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(parameters);
	Sophus::SE3 T = Sophus::SE3::exp(lie);
	for (int iter = 0; iter < iterations; iter++)
	{
		cost = 0;
		Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
		Eigen::Matrix< double, 6, 1 > g = Eigen::Matrix< double, 6, 1 >::Zero();

		for (int i = 0; i < xSet.size(); i++)
		{

			Eigen::Matrix<double, 3, 6> jaci = Eigen::Matrix<double, 3, 6>::Zero();
			jaci.block<3, 3>(0, 0) = -Eigen::Matrix3d::Identity();
			jaci.block<3, 3>(0, 3) = Sophus::SO3::hat(T * xSet[i]);


			

			Eigen::Vector3d dt =  - jaci.block<3, 3>(0, 0)  * T.log().head(3);


			Eigen::Vector3d dr = -jaci.block<3, 3>(0, 3) * T.log().tail(3);

			//std::cout << "dx " << dt <<   std::endl;
			//std::cout << "dr " << dr << std::endl;
			
			

			//std::cout << "jaci.block<3, 3>(0, 0) " << jaci.block<3, 3>(0, 0) << std::endl;

			//std::cout << "T.log().head(3) " << T.log().head(3) << std::endl;

			Sophus::Vector6d dx;

			dx.head(3) = dt;
			dx.tail(3) = dr;
			dx = -aph * dx;
			//std::cout << "dx " << dx << std::endl;

			T  = Sophus::SE3::exp(dx) * T;
			Eigen::Vector3d ei = ySet[i] - T * xSet[i];

			cost += (ei[0] * ei[0] + ei[1] * ei[1] + ei[2] * ei[2]);
	
		}

		std::cout << lastCost << " " << cost << std::endl;
		if (cost > lastCost)
		{

			//std::cout << "pig " << std::endl;
			//aph = aph * 0.1;
		}
		lastCost = cost;

		std::cout << "iter " << iter << " cost " << cost << std::endl;
	}

	std::cout << " cost=" << std::cout.precision(12) << cost << std::endl;
	std::cout << "4 " << std::endl;
	std::cout << "R " << T.rotation_matrix() << std::endl;
	std::cout << "tran " << T.translation().transpose() << std::endl;

}



void LM(const std::vector<Eigen::Vector3d >& xSet, const std::vector<Eigen::Vector3d>& ySet, double const* parameters)
{
	int iterations = 100000;
	double cost, costnew,lastCost = 0;

	Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(parameters);
	Sophus::SE3 T = Sophus::SE3::exp(lie);
	double tao = 1e-10;
	double u = 0;
	double v = 2;


	cost = 0; costnew = 0;
	Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
	Eigen::Matrix< double, 6, 1 > g = Eigen::Matrix< double, 6, 1 >::Zero();

	for (int i = 0; i < xSet.size(); i++)
	{

		Eigen::Vector3d ei = ySet[i] - T * xSet[i];

		Eigen::Matrix<double, 3, 6> jaci = Eigen::Matrix<double, 3, 6>::Zero();
		jaci.block<3, 3>(0, 0) = -Eigen::Matrix3d::Identity();
		jaci.block<3, 3>(0, 3) = Sophus::SO3::hat(T * xSet[i]);

		H += jaci.transpose() * jaci;
		g += -jaci.transpose() * ei;
		cost += (ei[0] * ei[0] + ei[1] * ei[1] + ei[2] * ei[2]);

	}

	Eigen::Matrix<double, 6, 6>  I = Eigen::Matrix<double, 6, 6>::Identity();
	u = tao * H.maxCoeff();


	int iter = 0;
	while(1)
	{
		cost = 0; costnew = 0;
		Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
		Eigen::Matrix< double, 6, 1 > g = Eigen::Matrix< double, 6, 1 >::Zero();

		for (int i = 0; i < xSet.size(); i++)
		{

			Eigen::Vector3d ei = ySet[i] - T * xSet[i];

			Eigen::Matrix<double, 3, 6> jaci = Eigen::Matrix<double, 3, 6>::Zero();
			jaci.block<3, 3>(0, 0) = -Eigen::Matrix3d::Identity();
			jaci.block<3, 3>(0, 3) = Sophus::SO3::hat(T * xSet[i]);

			H += jaci.transpose() * jaci;
			g += -jaci.transpose() * ei;
			cost += (ei[0] * ei[0] + ei[1] * ei[1] + ei[2] * ei[2]);

		}

		Eigen::Matrix<double, 6, 6>  I = Eigen::Matrix<double, 6, 6>::Identity();

		H = H + u * I ;
		Sophus::Vector6d dx;
		dx = H.ldlt().solve(g);

		Sophus::SE3 tempT = Sophus::SE3::exp(dx) * T;

		for (int i = 0; i < xSet.size(); i++)
		{

			Eigen::Vector3d ei = ySet[i] - tempT * xSet[i];
			costnew += (ei[0] * ei[0] + ei[1] * ei[1] + ei[2] * ei[2]);
		}

		double theta = 2 * (cost - costnew) / ((dx.transpose() * (u * dx + g)) + 1e-3);
		if (theta > 0)
		{
			T = tempT;
			u = u * std::max(0.3333333, (1 - (2 * theta - 1) * (2 * theta - 1) * (2 * theta - 1)));
			v = 2;
		}
		else
		{
			//不满足,扩大范围,更接近最速下降
			u = u * v;
			v = 2 * v;
		}
		if (cost < 2000) break;
		std::cout << "iter " << iter++ << " cost " << cost << std::endl;
		lastCost = cost;
	}

	std::cout << " cost=" << std::cout.precision(12) << cost << std::endl;
	std::cout << "4 " << std::endl;
	std::cout << "R " << T.rotation_matrix() << std::endl;
	std::cout << "tran " << T.translation().transpose() << std::endl;

}


void GN(const std::vector<Eigen::Vector3d >& xSet, const std::vector<Eigen::Vector3d>& ySet, double const* parameters)
{
	int iterations = 10000;
	double cost, lastCost = 0;

	Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(parameters);
	Sophus::SE3 T = Sophus::SE3::exp(lie); 
	double tao = 1e-10;
	for (int iter = 0; iter < iterations; iter++)
	{
		cost = 0;
		Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
		Eigen::Matrix< double, 6, 1 > g = Eigen::Matrix< double, 6, 1 >::Zero();

		for (int i = 0; i < xSet.size(); i++)
		{

			Eigen::Vector3d ei = ySet[i] - T * xSet[i];
			Eigen::Matrix<double, 3, 6> jaci = Eigen::Matrix<double, 3, 6>::Zero();
			jaci.block<3, 3>(0, 0) = -Eigen::Matrix3d::Identity();
			jaci.block<3, 3>(0, 3) = Sophus::SO3::hat(T * xSet[i]);

			H += jaci.transpose() * jaci;
			g += -jaci.transpose() * ei ;
			cost += (ei[0] * ei[0] + ei[1] * ei[1] + ei[2] * ei[2]);

		}

		Eigen::Matrix<double, 6, 6>  I = Eigen::Matrix<double, 6, 6>::Identity();
	

		Sophus::Vector6d dx;
		dx = H.ldlt().solve(g);
		
		//if (dx.norm() < 1e-15)
		//{
		//	break;
		//}

		T = Sophus::SE3::exp(dx) * T;

		//std::cout << "iter " << iter << " cost " << cost << std::endl;
		lastCost = cost;
	}

	std::cout << " cost=" << std::cout.precision(12) << cost << std::endl;
	std::cout << "4 " << std::endl;
	std::cout << "R " << T.rotation_matrix() << std::endl;
	std::cout << "tran " << T.translation().transpose() << std::endl;

}