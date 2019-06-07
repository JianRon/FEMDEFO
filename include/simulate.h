#ifndef SIMULATE_H
#define SIMULATE_H

#include <mesh_io.h>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <random>

using namespace std;

typedef Eigen::Matrix<double, 4, 3> Matrix43;
vector<Eigen::Matrix3d> B;
vector<double> W;
bool FirstSet = true;

int cnt;

// stop while any difference less than EPS
const double DAMPING = 0.0; // 
const double SPEED_MUL = 0.1;
const double VOL_MUL = 0.5; //¡¡0.2 error
const double EPS = 1e-6;
const double MIU = 1;
const double LAMBDA = 1;
const double GAMA = 1;
const double deltaT = 0.01;

void assium(Eigen::Vector3d& a, glm::vec3 b)
{
	a = Eigen::Vector3d(b.x, b.y, b.z);
}
void assium(glm::vec3& a, Eigen::Vector3d b)
{
	a = glm::vec3(b(0), b(1), b(2));
}

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
void getMatrixDeltaD(Eigen::Matrix3d& d, const vector<unsigned int>& u, const meshModel& model)
{
	const auto& Tet = model.mMesh.vertices;
	glm::vec3 i(Tet[u[0]].Deltax);
	glm::vec3 j(Tet[u[1]].Deltax);
	glm::vec3 k(Tet[u[2]].Deltax);
	glm::vec3 l(Tet[u[3]].Deltax);
	//cout << "? " << i.x << " " << i.y << " " << i.z << endl;
	d(0, 0) = (double)(i.x - l.x);d(0, 1) = (double)(j.x - l.x);d(0, 2) = (double)(k.x - l.x);
	d(1, 0) = (double)(i.y - l.y);d(1, 1) = (double)(j.y - l.y);d(1, 2) = (double)(k.y - l.y);
	d(2, 0) = (double)(i.z - l.z);d(2, 1) = (double)(j.z - l.z);d(2, 2) = (double)(k.z - l.z);
}

void getMatrixD(Eigen::Matrix3d& d, const Matrix43& x)
{
	Eigen::Vector3d i(x(0, 0), x(0, 1), x(0, 2));
	Eigen::Vector3d j(x(1, 0), x(1, 1), x(1, 2));
	Eigen::Vector3d k(x(2, 0), x(2, 1), x(2, 2));
	Eigen::Vector3d l(x(3, 0), x(3, 1), x(3, 2));

	d(0, 0) = (double)(i(0) - l(0));d(0, 1) = (double)(i(1) - l(1));d(0, 2) = (double)(i(2) - l(2));
	d(1, 0) = (double)(j(0) - l(0));d(1, 1) = (double)(j(1) - l(1));d(1, 2) = (double)(j(2) - l(2));
	d(2, 0) = (double)(k(0) - l(0));d(2, 1) = (double)(k(1) - l(1));d(2, 2) = (double)(k(2) - l(2));

}
// reset model
void reset(meshModel& model)
{
	model = meshModel(model.direction);
	FirstSet = true;
}



void getMatrixP(Eigen::Matrix3d& P, Eigen::Matrix3d F)
{

	/*for (int i = 0;i < 3;i++)
		for (int j = 0;j < 3;j++)
			cout << F(i, j) << (j == 2 ? '\n' : ' ');*/
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix<double, 3, 1> sigma = svd.singularValues();

	double I3 = 1;
	for (int i = 0;i < 3;i++) I3 *= sigma(i, 0) * sigma(i, 0);

	//cout << I3 << endl;
	Eigen::Matrix3d ITF = F.inverse().transpose();
	P = MIU * (F - ITF) + LAMBDA * log(I3) / 2.0 * ITF;
}

void getDeltaP(Eigen::Matrix3d& DeltaP, Eigen::Matrix3d F, Eigen::Matrix3d DeltaF, int flag = 1)
{
	if (flag == 1) {
		// Stress differentials for St. Venant-Kirchhoff materials
		// -------------------------------------------------------

		Eigen::Matrix3d E, identity, DeltaE;
		identity = Eigen::Matrix3d::Identity(3, 3);

		DeltaE = 0.5 * (DeltaF.transpose() * F + F.transpose() * DeltaF);

		E = 0.5 * (F.transpose()*F - identity);

		DeltaP = DeltaF * (2 * MIU * E + LAMBDA * E.trace() * identity) +
			F * (2 * MIU * DeltaE + LAMBDA * DeltaE.trace() * identity);
	}
	else
	{
		// Stress differentials for Neohookean materials
		// ---------------------------------------------
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix<double, 3, 1> sigma = svd.singularValues();

		double I3 = 1;
		for (int i = 0;i < 3;i++) I3 *= sigma(i, 0) * sigma(i, 0);
		double J = sqrt(I3);
		//cout << I3 << endl;
		Eigen::Matrix3d ITF = F.inverse().transpose();
		DeltaP = MIU * DeltaF + (MIU - LAMBDA * log(J)) * ITF*DeltaF.transpose()*ITF +
			LAMBDA * (F.inverse() * DeltaF).trace() * ITF;
	}
}

vector<glm::vec3> computerForceDifferenials(meshModel model) {
	vector<glm::vec3> res;
	for (auto& u : model.mMesh.vertices)
	{
		u.force  = glm::vec3(0.0f);
		u.Deltaf = glm::vec3(0.0f);
	}
	unsigned int index = 0;
	for (auto& v : model.mTetrahedras)
	{
		Eigen::Matrix3d D;
		getMatrixD(D, v.vertices, model);

		Eigen::Matrix3d DeltaD;
		getMatrixDeltaD(DeltaD, v.vertices, model);

		Eigen::Matrix3d F;
		F = D * B[index];

		Eigen::Matrix3d DeltaF;
		DeltaF = DeltaD * B[index];

		Eigen::Matrix3d DeltaP;
		getDeltaP(DeltaP, F, DeltaF);

		Eigen::Matrix3d DeltaH;
		DeltaH = -1 * W[index] * DeltaP * B[index].transpose();

		glm::vec3 f1(DeltaH(0, 0), DeltaH(1, 0), DeltaH(2, 0));
		glm::vec3 f2(DeltaH(0, 1), DeltaH(1, 1), DeltaH(2, 1));
		glm::vec3 f3(DeltaH(0, 2), DeltaH(1, 2), DeltaH(2, 2));

		auto& vertex = model.mMesh.vertices;
		const auto& u = v.vertices;
		vertex[u[0]].Deltaf += f1;
		vertex[u[1]].Deltaf += f2;
		vertex[u[2]].Deltaf += f3;
		vertex[u[3]].Deltaf -= f1 + f2 + f3;
		index++;
	}
	for (auto& u : model.mMesh.vertices)
		res.push_back(u.Deltaf);
	return res;
}
vector<glm::vec3> computerElasticForces(meshModel model) {
	vector<glm::vec3> res;
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
		H = -1 * W[index] * P * B[index].transpose();
		index++;

		auto& vertex = model.mMesh.vertices;
		const auto& u = v.vertices;

		glm::vec3 f1(H(0, 0), H(1, 0), H(2, 0));
		glm::vec3 f2(H(0, 1), H(1, 1), H(2, 1));
		glm::vec3 f3(H(0, 2), H(1, 2), H(2, 2));

		vertex[u[0]].force += f1;
		vertex[u[1]].force += f2;
		vertex[u[2]].force += f3;
		vertex[u[3]].force -= f1 + f2 + f3;
	}
	for (auto& u : model.mMesh.vertices)
		res.push_back(u.force);
	return res;
}

bool ok(meshModel& model)
{
	for (auto& u : model.mMesh.vertices)
	{
		if (u.Deltax.x > EPS || u.Deltax.y > EPS || u.Deltax.z > EPS) return true;
	}
	return false;
}
vector<glm::vec3> ConjugateGradient(meshModel model, const meshModel& old)
{
	vector<glm::vec3> res, fe, fd;
	res.resize(model.numVertices, glm::vec3(0));
	vector<Eigen::Vector3d> R, P, KW, AP;
	meshModel temp = model;
	fe = computerElasticForces(temp);
	for (auto& u : temp.mMesh.vertices)
	{
		u.Deltax.x = GAMA * u.velocity.x;
		u.Deltax.y = GAMA * u.velocity.y;
		u.Deltax.z = GAMA * u.velocity.z;
	}
	fd = computerForceDifferenials(temp);

	for (unsigned int i = 0;i != model.mMesh.vertices.size();i++) {
		glm::vec3 b = old.mMesh.vertices[i].position - model.mMesh.vertices[i].position;
		b.x *= MASS / deltaT;
		Eigen::Vector3d vb;
		assium(vb, b);
		R.push_back(vb);
		P.push_back(vb);
	}

	int T = 5;
	while (T--)
	{
		unsigned int index = 0;
		for (auto&u : model.mMesh.vertices)
			assium(u.Deltax, -P[index]);

		vector<glm::vec3> KW = computerForceDifferenials(model);
		for (auto&u : model.mMesh.vertices)
		{
			Eigen::Vector3d tmp;
			glm::vec3       tmpv;
			assium(tmp, KW[index]);
			double ta, c, d;
			c = 1 + GAMA / deltaT;
			d = 1 / deltaT / deltaT;
			ta = (R[index].transpose() * R[index]);

			ta /= (P[index].transpose() * (c * tmp + d * P[index]));

			assium(tmpv, ta * P[index]);
			res[index] = res[index] + tmpv;
			Eigen::Vector3d oldR;
			oldR = R[index];
			R[index] = R[index] - ta * tmp;

			double beta = (R[index].transpose() * R[index]);
			beta /= (oldR.transpose() * oldR);
			P[index] = R[index] + beta * P[index];

			index++;
		}
	}
	return res;
}

void simulate(meshModel& model, float deltaTime = 0.1, int flag = 1)
{

	cout << "SIMULATING..." << endl;

	meshModel now = model,old = model;

	int T = 5;
	while (T--)
	{
		vector<glm::vec3> dx = ConjugateGradient(now, old);
		for (unsigned int i = 0;i != dx.size();i++) {
			now.mMesh.vertices[i].position = old.mMesh.vertices[i].position + dx[i];
		}
		old = now;
	}
	check(now);
	model = now;
	cout << "SIMULATED OVER" << endl;

	// check(model);
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
		// simulate(model, (float)deltaT, 1);
		vector<glm::vec3> F = computerElasticForces(model);
		unsigned int  index = 0;
		for (auto& u : model.mMesh.vertices) {
			Eigen::Vector3d temp;
			assium(temp, F[index]);
			assium(u.velocity, deltaT * temp);
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

#endif // !SIMULATE_H
