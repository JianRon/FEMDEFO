#ifndef ME	SHIO_H
#define MESHIO

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <shader.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include <Eigen/Core>
using namespace std;

const double MASS = 1;

struct Vertex
{
	glm::vec3 position;
	glm::vec3 force;
	glm::vec3 velocity;
	glm::vec3 normal;
	double mass;

	Vertex()
	{
		mass = MASS;
		velocity = glm::vec3(0.0f);
		force = glm::vec3(0);
		normal = glm::vec3(0);
	}
};

struct Mesh
{
public:
	vector<unsigned int> indices;
	vector<Vertex> vertices;
	// vector<Eigen::Vector3d> normals;
	unsigned int VAO;
	bool firstFlag;
	Mesh() { firstFlag = true; }
	~Mesh() 
	{
		indices.clear();
		vertices.clear();
	}
	void Draw(Shader shader)
	{
		if (firstFlag)
		{
			setupMesh();
			setTexture();
			firstFlag = false;
		}
		setupVBO(VBO);
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);

		// set texture to default
		glActiveTexture(GL_TEXTURE0);
	}
	void getNormal()
	{
		for (auto& u : vertices) u.normal = glm::vec3(0);

		for (int i = 0;i * 3 < (unsigned int)indices.size();i++)
		{
			glm::vec3 a = vertices[indices[i * 3 + 0]].position;
			glm::vec3 b = vertices[indices[i * 3 + 1]].position;
			glm::vec3 c = vertices[indices[i * 3 + 2]].position;
			
			glm::vec3 ab = b - a;
			glm::vec3 ac = c - a;

			Eigen::Vector3d AB,AC,NR;
			
			AB << ab.x,ab.y,ab.z;
			AC << ac.x,ac.y,ac.z;
			NR = AB.cross(AC);
			NR.normalize();
			
			// normals.push_back(NR);
			for (int j = 0;j < 3;j++)
				vertices[indices[i * 3 + j]].normal.x += NR(0),
				vertices[indices[i * 3 + j]].normal.y += NR(1),
				vertices[indices[i * 3 + j]].normal.z += NR(2);
		}
	}
private:
	unsigned int VBO, EBO;
	void setTexture()
	{
		// TODO
	}
	void setupMesh()
	{
		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &EBO);
		setupVAO(VAO);
		setupVBO(VBO);
		setupEBO(EBO);

		// vertex position
		// ---------------
		bindAttrib(0, 3, sizeof(Vertex), (void*)0);
		bindAttrib(1, 3, sizeof(Vertex), (void*)offsetof(Vertex, normal));
		// take off
		// --------
		glBindVertexArray(0);
	}
	// bind attrib of vectices
	void bindAttrib(int location, GLint size, GLsizei sizeType,const void* offset)
	{
		glEnableVertexAttribArray(location);
		glVertexAttribPointer(location, size, GL_FLOAT, GL_FALSE, sizeType, offset);
	}
	// setup VAO
	void setupVAO(unsigned int& VAO)
	{
		glBindVertexArray(VAO);
	}
	// setup VBO
	void setupVBO(unsigned int& VBO)
	{
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), &vertices[0], GL_DYNAMIC_DRAW);

	}
	// setup EBO
	void setupEBO(unsigned int& EBO)
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_DYNAMIC_DRAW);
	}
};

struct Tetrahedra
{
public:
	vector<unsigned int> vertices;
	Tetrahedra() {}
	~Tetrahedra() { vertices.clear(); }
};

struct meshModel
{
public:
	string direction;
	vector<Tetrahedra> mTetrahedras;
	Mesh mMesh;
	unsigned int numTetrahedra;
	unsigned int numMesh;
	unsigned int numVertices;
	
	// Initializes(load) model
	~meshModel()
	{
		mTetrahedras.clear();
	}
	meshModel(const string& path)
	{
		direction = path;
		loadModel(path);
	}
	void Draw(Shader shader)
	{
		mMesh.Draw(shader);
	}
private:
	void loadModel(const string& path)
	{
		
		fstream inFile(path);
		if (!inFile)
		{
			cerr << "ERROR::MESHIO::LOAD FAIL" << endl;
		}
		string typeName;
		while (inFile.peek() != EOF)
		{
			inFile >> typeName;
			int ignore;

			if (typeName == "Vertices")
			{
				inFile >> numVertices;
				mMesh.vertices.resize(numVertices);
				for (auto& u : mMesh.vertices)
					inFile >> u.position.x >> u.position.y >> u.position.z >> ignore;
			}
			else if (typeName == "Tetrahedra")
			{
				inFile >> numTetrahedra;
				mTetrahedras.resize(numTetrahedra);
				for (auto& u : mTetrahedras)
				{
					u.vertices.resize(4);
					for (auto& v : u.vertices)
					{
						inFile >> v;
						v--;
					}
					inFile >> ignore;
				}
			}
			else if (typeName == "Triangles")
			{
				inFile >> numMesh;
				mMesh.indices.resize(numMesh * 3);
				int cnt = 0;
				for (auto& u : mMesh.indices)
				{
					cnt++;
					int ig;
					
					if (cnt % 4 == 0)
					{
						cnt++;
						inFile >> ig;
						// cout << ig << endl;
					}
						inFile >> u;
					u--;
				}
			}
		}
		mMesh.getNormal();
		
	}
};

#endif // !MESHIO_H
