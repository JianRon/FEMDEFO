#ifndef SIMULATE_H
#define SIMULATE_H

#include <mesh_io.h>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>

using namespace std;

vector<Eigen::Matrix3d> B;
vector<double> W;
bool FirstSet = true;
int cnt;

// stop while any difference less than EPS
const double DAMPING = 0; // 
const double SPEED_MUL = 6;
const double VOL_MUL = 0.5; //¡¡0.2 error
const double EPS = 1e-6;
const double miu = 0.5;
const double lambda = 1;

void check(meshModel& model)
{
	cout << " cnt -> " << cnt << endl;
	for (auto& u : model.mMesh.vertices)
	{

		cout << "force (" << u.force.x << ","
			<< u.force.y << ","
			<< u.force.z << ")" << endl;
		cout << "postion (" << u.position.x << ","
			<< u.position.y << ","
			<< u.position.z << ")" << endl;
		cout << "velocity (" << u.velocity.x << ","
			<< u.velocity.y << ","
			<< u.velocity.z << ")" << endl;
	}
	cout << endl;
}

void getMatrixD(Eigen::Matrix3d& d, const vector<unsigned int>& u, const meshModel& model)
{
	const auto& Tet = model.mMesh.vertices;
	glm::vec3 i(Tet[u[0]].position);
	glm::vec3 j(Tet[u[1]].position);
	glm::vec3 k(Tet[u[2]].position);
	glm::vec3 l(Tet[u[3]].position);
	//cout << "? " << i.x << " " << i.y << " " << i.z << endl;
	d(0, 0) = (double)(i.x - l.x);d(0, 1) = (double)(j.x - l.x);d(0, 2) = (double)(k.x - l.x);
	d(1, 0) = (double)(i.y - l.y);d(1, 1) = (double)(j.y - l.y);d(1, 2) = (double)(k.y - l.y);
	d(2, 0) = (double)(i.z - l.z);d(2, 1) = (double)(j.z - l.z);d(2, 2) = (double)(k.z - l.z);
}
// reset model
void reset(meshModel& model)
{
	model = meshModel(model.direction);
	FirstSet = true;
}

// make vertices radom distribute in a cube
void random(meshModel& model)
{
	// reset force;
	cnt++;
	if (FirstSet)
	{
		FirstSet = false;
		// precomputation
		// --------------
		for (auto& v : model.mTetrahedras)
		{
			Eigen::Matrix3d D;
			getMatrixD(D, v.vertices, model);
			B.push_back(D.inverse());
			W.push_back((double)(1.0 / 6.0) * fabs(D.determinant()));
		}
	}
	
	for (auto& u : model.mMesh.vertices)
		u.position *= VOL_MUL;
	Sleep(1000);
	/*default_random_engine generater;
	uniform_real_distribution<double> distribution(-1,1);
	for (auto& u : model.mMesh.vertices)
		u.position = glm::vec3(distribution(generater),distribution(generater),distribution(generater));*/
	// check(model);
}


void getMatrixP(Eigen::Matrix3d& P, Eigen::Matrix3d F)
{

	/*for (int i = 0;i < 3;i++)
		for (int j = 0;j < 3;j++)
			cout << F(i, j) << (j == 2 ? '\n' : ' ');*/
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix<double,3,1> sigma = svd.singularValues();
	
	double I3 = 1;
	for (int i = 0;i < 3;i++) I3 *= sigma(i, 0) * sigma(i, 0);
	
	//cout << I3 << endl;
	Eigen::Matrix3d ITF = F.inverse().transpose();
	P = miu * (F - ITF) + lambda * log(I3) / 2.0 * ITF;
}

void move(meshModel& model, float deltatime)
{
	glm::vec3 acceletation;
	deltatime *= SPEED_MUL;
	float powDeltatime = deltatime * deltatime;

	for (auto& u : model.mMesh.vertices)
	{
		
		glm::vec3 realForce = u.force;
		realForce.x -= u.velocity.x * DAMPING;
		realForce.y -= u.velocity.y * DAMPING;
		realForce.z -= u.velocity.z * DAMPING;

		acceletation.x = realForce.x / u.mass;
		acceletation.y = realForce.y / u.mass;
		acceletation.z = realForce.z / u.mass;
		/*double invMass = 1.0 / u.mass;
		glm::mat4 temp(1.0f);
		temp = glm::scale(temp, glm::vec3(invMass, invMass, invMass));
		acceletation = temp * glm::vec4(u.force,1.0f);*/

		// s = 1/2 * a * t^2 + v * t;
		glm::mat4 Time(1.0f),powTime(1.0f);
		Time = glm::scale(Time, glm::vec3(deltatime, deltatime, deltatime));
		powTime = glm::scale(powTime, glm::vec3(powDeltatime / 2, powDeltatime / 2, powDeltatime / 2));
		glm::vec4 Temp;
		Temp = Time * glm::vec4(u.velocity, 1.0f) + powTime * glm::vec4(acceletation, 1.0f);
		
		u.position.x += Temp.x;
		u.position.y += Temp.y;
		u.position.z += Temp.z;

		u.velocity += acceletation * deltatime;
	}
}


void simulate(meshModel& model,float deltaTime = 0.1)
{

	cout << "SIMULATING..." << endl;
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
		// cout << P * B[index].transpose() << endl;
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
	cout << "SIMULATED OVER" << endl;
	// check(model);
}
#endif // !SIMULATE_H
