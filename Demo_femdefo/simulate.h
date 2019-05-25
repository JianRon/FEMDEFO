#ifndef SIMULATE_H
#define SIMULATE_H

#include <mesh_io.h>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>
using namespace std;

// stop while any difference less than EPS
const double EPS = 1e-6;

vector<Eigen::Matrix3d> B;
vector<double> W;
bool firstSet = true;
// reset model
void reset(meshModel& model)
{

	// precomputation
	// --------------
	if (firstSet)
	{
		for (auto& v : model.mTetrahedras)
		{
			Eigen::Matrix3d D;
			getMatrixD(D, v.vertices, model);
			B.push_back(D.inverse());
			W.push_back((double)(1.0 / 6.0) * D.determinant());
		}
		firstSet = false;
	}
	model = meshModel(model.direction);
}

// make vertices radom distribute in a cube
void random(meshModel& model)
{

	/*default_random_engine generator;
	uniform_real_distribution<double> distribution(-0.01,0.01);
	for (auto& u : model.mMesh.vertices)
		u.position.x = distribution(generator),
		u.position.y = distribution(generator),
		u.position.z = distribution(generator);*/
	/*
	for (auto& u : model.mMesh.vertices)
		u.position.x /= 10,
		u.position.y /= 10,
		u.position.z /= 10;
	*/
}

void getMatrixD(Eigen::Matrix3d& d,const vector<unsigned int>& u,const meshModel& model)
{
	const auto& Tet = model.mMesh.vertices;
	glm::vec3 i(Tet[u[0]].position);
	glm::vec3 j(Tet[u[1]].position);
	glm::vec3 k(Tet[u[2]].position);
	glm::vec3 l(Tet[u[3]].position);

	d(0, 0) = (double)(i.x - l.x);d(0, 1) = (double)(j.x - l.x);d(0, 2) = (double)(k.x - l.x);
	d(1, 0) = (double)(i.y - l.y);d(1, 1) = (double)(j.y - l.y);d(1, 2) = (double)(k.y - l.y);
	d(2, 0) = (double)(i.z - l.z);d(2, 1) = (double)(j.z - l.z);d(2, 2) = (double)(k.z - l.z);

}
void getMatrixP(Eigen::Matrix3d& P, Eigen::Matrix3d F)
{
	double miu = 0.1;
	double lambda = 0.1;
	
	
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d sigma = svd.singularValues();
	sigma *= sigma;
	P = miu * (F - miu * F.inverse().transpose());
}

void move(meshModel& model, float deltatime)
{
	glm::vec3 acceletation;
	float powDeltatime = deltatime * deltatime;

	for (auto& u : model.mMesh.vertices)
	{
		double invMass = 1.0 / u.mass;
		glm::mat4 temp(1.0f);
		temp = glm::scale(temp, glm::vec3(invMass, invMass, invMass));
		glm::vec4 tvec;
		tvec = temp * glm::vec4(u.force,1.0f);
		u.force = tvec;

		// s = 1/2 * a * t^2 + v * t;
		glm::mat4 Time(1.0f),powTime(1.0f);
		Time = glm::scale(Time, glm::vec3(deltatime, deltatime, deltatime));
		powTime = glm::scale(powTime, glm::vec3(powDeltatime / 2, powDeltatime / 2, powDeltatime / 2));
		glm::vec4 Temp;
		Temp = Time * glm::vec4(u.velocity, 1.0f) + powTime * glm::vec4(acceletation, 1.0f);
		u.position.x += Temp.x;
		u.position.y += Temp.y;
		u.position.z += Temp.z;
		u.velocity = acceletation * deltatime;
	}
}

void simulate(meshModel& model,float deltaTime)
{
	// reset force;
	for (auto& u : model.mMesh.vertices)
		u.force = glm::vec3(0.0f);

	// compute force
	// -------------
	unsigned int index = 0;
	for (auto& v : model.mTetrahedras)
	{
		Eigen::Matrix3d D;
		getMatrixD(D, v.vertices, model);

		Eigen::Matrix3d F;
		F = D * B[index];
		
		Eigen::Matrix3d P;
		getMatrixP(P, F);

		Eigen::Matrix3d H;
		H = -1 * W[index] * P *B[index].transpose();
		index++;

		auto& vertex = model.mMesh.vertices;
		const auto& u = v.vertices;
		
		glm::vec3 f1(H(0,0),H(1,0),H(2,0));
		glm::vec3 f2(H(0,1),H(1,1),H(2,1));
		glm::vec3 f3(H(0,2),H(1,2),H(2,2));
		
		vertex[u[0]].force += f1;
		vertex[u[1]].force += f2;
		vertex[u[2]].force += f3;
		vertex[u[3]].force -= f1 + f2 + f3;	
	}
	move(model, deltaTime);
}
#endif // !SIMULATE_H
