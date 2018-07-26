

#include "Model.h"
#include "GameObject.h"

using namespace std;



void GameObject::update(float dt, glm::mat4 PV)
{	
	//Summing all forces acting on the object
	addForces();
	//Calculating position based off kinematics equations
	eulerIntegration(dt);
	//translating object
	translation = glm::translate(glm::vec3(position));
	//Setting rotation
	setRotation(angP);
	//Calculating model matrix and then updating the objects MVP
	calcTransform();
	MVP = PV * transformation;
	//Zeroing force and acceleration
	totalForce = acceleration = torque = angA = glm::vec3();



	
	
}

void GameObject::addForces()
{
	//totalForce += GRAVITY * mass;

	if (position.y < -1)
	{
		position.y = -1;

		if (velocity.y < 0)
			totalForce += velocity * -1.0f; // Bounce force
	}



	acceleration = (totalForce / mass); //Force=m * a so a = F/m
	angA = torque / inertia; //Torque= I*ang.a so ang.a =T/I
}

//Calculates Transformation matrix -> T * R * S
void GameObject::calcTransform()
{
	transformation = translation * rotation * scale;

}

//Updates velocity based off pure Euler integration
void GameObject::eulerIntegration(float dt)
{
	//LINEAR

	//Kinematic equation x= v*dt+ 0.5(a *dt^2)
	glm::vec3 posFristdegree = velocity * dt;
	position += posFristdegree + 0.5f*(acceleration * dt*dt);

	

	//Updating velocity
	velocity += acceleration * dt;

	//ROTATIONAL

	//Same equation as the linear one, except the variables are switched to their angular counterparts
	angP += angV * dt + 0.5f*(angA*dt*dt);

	angV += angA * dt;//Updating velocity 

}
// + netImpulse * netForce

void GameObject::setModel(Model *m)
{
	mesh = m;
}


//Return global vertices
std::vector<glm::vec3> GameObject::Vertices()
{

	std::vector<glm::vec3> vertices;

	for (int i = 0; i < mesh->NumVertices(); i++)
	{
		if (mesh->Vertices()[i].position != glm::vec3(0, 0, 0))
		{
			glm::vec4 vertex = transformation * glm::vec4(mesh->Vertices()[i].position.x, mesh->Vertices()[i].position.y, 0, 1);
			vertices.push_back(glm::vec3(vertex.x,vertex.y,0));
		}
	}
	return vertices;

}

//Rotates by x,y and z radians
void GameObject::rotate(glm::vec3 angle)
{
	glm::quat qDelta = glm::quat(angle);

	currentRot *= qDelta;

	rotation = glm::toMat4(currentRot);

	calcTransform();
}

// Sets rotation by x, y and z radians
void GameObject::setRotation(glm::vec3 angle)
{
	currentRot = glm::quat(angle);

	rotation = glm::toMat4(currentRot);

	calcTransform();
}

//Initalizing values
GameObject::GameObject(Model* m)
{
	velocity = totalForce = acceleration = position = glm::vec3();

	translation = rotation = scale = transformation = glm::mat4();

	angA = angP = angV = glm::vec3();

	netForce = glm::vec3(0.0f);
	netImpulse = glm::vec3(0.0f);

	dynamicFriction = 1.0f;
	staticFriction = 1.0f;
	mass = 1.0f;
	mesh = m;
	radius = 0.2f;

	//inverseMass = 1.0f;
	//restitution = 1.0f;

	
	

}







GameObject::~GameObject()
{

}
