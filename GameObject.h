#ifndef _GAME_OBJECT_H
#define _GAME_OBJECT_H

/*
File Name : Gameobject.h
Copyright © 2018
Original authors : Sanketh Bhat
Written under the supervision of David I.Schwartz, Ph.D., and
supported by a professional development seed grant from the B.Thomas
Golisano College of Computing & Information Sciences
(https ://www.rit.edu/gccis) at the Rochester Institute of Technology.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.

Description:

*/

#include "Model.h"
#include "GLIncludes.h";



using namespace std;

class GameObject
{
public:


	const glm::vec3 GRAVITY = glm::vec3(0, -0.00000000000001, 0); // No gravity here

	glm::vec3 position;
	glm::vec3 velocity;
	glm::vec3 acceleration;

	//float inverseMass;			//I tend to use inverse mass as opposed to mass itself. It saves lots of divides when forces are involved.
	//float restitution;			//How elastic this object is (1.0f is perfectly elastic, 0.0f is perfectly inelastic)
	float dynamicFriction;		//Dynamic coefficient of Friction
	float staticFriction;       //Static coefficient of Friction

	glm::vec3 netForce;			//Forces over time
	glm::vec3 netImpulse;		//Instantaneous forces

	//Rotational Variables
	glm::vec3 angP;
	glm::vec3 angA;
	glm::vec3 torque;
	glm::vec3 angV;


	float inertia;
	Model *mesh;
	float mass;
	float radius = 0.02f;


	glm::mat4 translation;
	glm::mat4 rotation;
	glm::mat4 scale;
	glm::mat4 transformation; //Model matrix


	glm::quat currentRot;
	glm::vec3 totalForce;
	glm::mat4 MVP;
	int integType = 1;
	float width;
	float height;


	void update( float dt, glm::mat4 PV);

	void addForces();

	void calcTransform();

	void eulerIntegration(float dt);

	void setModel(Model *m);

	void rotate(glm::vec3 angle);

	void setRotation(glm::vec3 angle);


	Model* model()
	{
		return mesh;
	}
	
	void render(GLuint uniMVP)
	{
		glUniformMatrix4fv(uniMVP, 1, GL_FALSE, glm::value_ptr(MVP));
		mesh->Draw();
	}

	glm::vec3 Position()
	{
		return position;
	}

	glm::vec3 Position(glm::vec3 p)
	{
		position = p;
		return position;
	}

	glm::vec3 Velocity()
	{
		return velocity;
	}

	glm::vec3 Velocity(glm::vec3 v)
	{
		velocity = v;
		return velocity;
	}

	glm::vec3 Acceleration()
	{
		return acceleration;

	}

	glm::vec3 Acceleration(glm::vec3 a)
	{
		acceleration = a;
		return acceleration;
	}


	

	std::vector<glm::vec3> Vertices();

	GameObject(Model* m);
	~GameObject();
};



#endif _GAME_OBJECT_H