#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

int cal()
{
	Eigen::Matrix3d mat_3d(1,2,3);	
	for (int i = 0;i < 3;i++)
		for (int j = 0;j < 3;j++)
			cout << mat_3d(i, j) << (j == 2 ? '\n' : ' ');

	return 0;
}