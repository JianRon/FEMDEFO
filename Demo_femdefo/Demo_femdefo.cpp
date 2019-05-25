#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <shader.h>
#include <camera.h>
//#include <mesh.h>
//#include <model.h>
//#include <mesh_io.h>
#include <simulate.h>
#include <iostream>
const int SCR_WIDTH = 800;
const int SCR_HEIGHT = 600;
void framebuffSizeCallback(GLFWwindow* window, int width,int height);
void mouseCallback(GLFWwindow* window,double xpos,double ypos);
void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

float lastFrame = 0;
float deltaTime = 0;
double lastX = SCR_WIDTH / 2.0;
double lastY = SCR_HEIGHT / 2.0;

bool firstMouse = true;

Camera camera(glm::vec3(0.0f, 0.0f, 10.0f));

meshModel ourModel("reference\\ditto.mesh");

int main()
{
	//for (auto& u : ourModel.mMesh.vertices)
	//	u.position.x *= 50,
	//	u.position.y *= 50,
	//	u.position.z *= 50;
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE,GLFW_OPENGL_CORE_PROFILE);

	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH,SCR_HEIGHT, "FEMDEFO",NULL,NULL);
	
	// config glfw window
	// ------------------
	glfwMakeContextCurrent(window);
	glfwSetFramebufferSizeCallback(window, framebuffSizeCallback);
	glfwSetCursorPosCallback(window, mouseCallback);
	glfwSetScrollCallback(window, scrollCallback);

	// capture mouse
	// -------------
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	
	// load glad
	// ---------
	gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
	
	// configure globle opengl state
	// -----------------------------
	glEnable(GL_DEPTH_TEST);

	// create shader
	// -------------
	Shader ourShader("shader\\shaders.vs", "shader\\shaders.fs");
	
	// load model
	// ----------
	
	
	/*for (int i = 0;i < ourModel.numVertices;i++)
	{
		
		std::cout << ourModel.mMesh.vertices[i].position.x << " "
			<< ourModel.mMesh.vertices[i].position.y << " "
			<< ourModel.mMesh.vertices[i].position.z << " "
			<< std::endl;
	}
	for (int i = 0;i < ourModel.numMesh;i++)
	{
		std::cout << ourModel.mMesh.indices[i * 3] << " "
			<< ourModel.mMesh.indices[i * 3 + 1] << " "
			<< ourModel.mMesh.indices[i * 3 + 2] << " "
			<< std::endl;
	}*/
	//Model ourModel("reference\\nanosuit.obj");
	// cout << " numMesh " << ourModel.numMesh << " numVertices " << ourModel.numVertices << " numTet" << ourModel.numTetrahedra;
	//while (true);
	while (!glfwWindowShouldClose(window))
	{
		float currentFrame = glfwGetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;

		glm::mat4 projection;
		glm::mat4 view;
		glm::mat4 model;

		processInput(window);

		glClearColor(0.1f, 0.3f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		ourShader.use();
		ourShader.setVec3("objectColor", 0.2f, 0.2f, 0.2f);
		ourShader.setVec3("lightColor", 1.0f, 1.0f, 1.0f);
		ourShader.setVec3("lightPos", -2.0f, 2.0f, 3.0f);
		ourShader.setVec3("viewPos", camera.Position);

		projection = glm::perspective(glm::radians(camera.Zoom), (float)(SCR_WIDTH) / (SCR_HEIGHT), 0.1f, 100.0f);
		view = camera.GetViewMatrix();
		ourShader.setMat4("projection", projection);
		ourShader.setMat4("view", view);

		model = glm::mat4(1.0f);
		model = glm::translate(model, glm::vec3(0.0f, -1.75f, 0.0f)); // translate it down so it's at the center of the scene
		model = glm::scale(model, glm::vec3(1.0f, 1.0f, 1.0f));	// it's a bit too big for our scene, so scale it down
		ourShader.setMat4("model", model);
		ourModel.Draw(ourShader);
		
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	glfwTerminate();
	return 0;
}
void processInput(GLFWwindow *window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
		glfwSetWindowShouldClose(window, true);
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		camera.ProcessKeyboard(FORWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		camera.ProcessKeyboard(BACKWARD, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		camera.ProcessKeyboard(LEFT, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		camera.ProcessKeyboard(RIGHT, deltaTime);
	if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
		reset(ourModel);
	if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS)
		random(ourModel);
	if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
		simulate(ourModel);
	if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
		check(ourModel);
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		for (int i = 0;i < 10;i++)
		{
			cout << i + 1 << endl;
			simulate(ourModel, deltaTime);
		}
}

void framebuffSizeCallback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void mouseCallback(GLFWwindow* window, double xpos, double ypos)
{
	if(firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}
	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos;

	lastX = xpos;
	lastY = ypos;

	camera.ProcessMouseMovement(xoffset, yoffset);
}
void scrollCallback(GLFWwindow* window, double xoffset, double yoffset)
{
	camera.ProcessMouseScroll(yoffset);
}