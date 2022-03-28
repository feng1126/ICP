#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../eigen/Eigen/Eigen"

struct data
{
	double timeStamp;
	Eigen::Vector3d xyz;

};

std::vector<std::string> splitString(const std::string& strs)
{
	std::string temp;
	std::vector<std::string> splitOut;
	splitOut.clear();
	for (int i = 0; i < strs.size(); ++i)
	{
		if (strs[i] == ',')
		{
			splitOut.push_back(temp);
			temp = "";
		}
		else
		{
			temp = temp + strs[i];
		}
		if (i == strs.size() - 1)
		{
			splitOut.push_back(temp);
		}
	}
	return splitOut;
}

std::vector<data > readData(std::string path)
{
	std::string strs;
	std::ifstream in_fp;
	in_fp.open(path.c_str());
	std::vector<data > datas;
	datas.clear();
	while (std::getline(in_fp, strs))
	{
		data temp;
		std::vector<std::string> logData = splitString(strs);
		temp.timeStamp  = stod(logData[0]);
		temp.xyz = Eigen::Vector3d(stod(logData[1]),stod(logData[2]),stod(logData[3]));
		datas.push_back(temp);
	}

	return datas;
}

template <typename Derived>
static Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymmetric(const Eigen::MatrixBase<Derived>& q)
{
	Eigen::Matrix<typename Derived::Scalar, 3, 3> ans;
	ans << typename Derived::Scalar(0), -q(2), q(1),
		q(2), typename Derived::Scalar(0), -q(0),
		-q(1), q(0), typename Derived::Scalar(0);
	return ans;
}