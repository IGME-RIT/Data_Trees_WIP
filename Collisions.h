#ifndef _COLLISIONS_H
#define _COLLISIONS_H

#include "GameObject.h"

using namespace std;

bool CheckCollision(GameObject* c1, GameObject* c2)
{
	float dist = glm::length(c1->Position() - c2->Position());
	return (0.03f) > dist;
}



void resolveCollisionimpulse(GameObject* A, GameObject* B)
{
	//This function resolves the collision between two particles.
	//There are two ways to resolve collision, using impulse based systems
	//or using the momentum and energy equations to calculate the final velocities.

	glm::vec3 n = (B->Velocity() - A->Velocity());
	//glm::vec3 n = B->Position() - A->Position();

	if (glm::length(n) > FLT_EPSILON)
	{
		n = glm::normalize(n);

		//This section of commented code decouples the two particles
		// We have commented out this code, as it makes the demo more jittery
		//A->Position() += 1 / 2.0f * n;
		//B->Position() -= 1 / 2.0f * n;

		glm::vec3 An = glm::dot(A->Velocity(), n) * n;
		glm::vec3 Bn = glm::dot(B->Velocity(), n) * n;

		float nMagA = fabsf(glm::dot(A->Velocity(), n));
		float nMagB = fabsf(glm::dot(B->Velocity(), n));

		// Apn and Bpn store those component of velocity which will not be affected in this collision
		glm::vec3 Apn = A->Velocity() - An;
		glm::vec3 Bpn = B->Velocity() - Bn;

		float relVel = nMagA - nMagB;

		float num;
		num = -2.0f * relVel;
		float denom;
		denom = A->mass + B->mass; //1.0f/A.mass + 1.0f/B.mass;

								   // Now that we have the velocities in that one dimension, solve them and add the components to the respective
								   // resultant to get the final velocities.
		glm::vec3 u1, u2;
		u1 = An;
		u2 = Bn;

		float j = num / denom;
		//A->Velocity() += j * n / A->mass;
		//B->Velocity() -= j * n / B->mass;

		A->Velocity(((((A->mass - B->mass)*u1) + 2.0f * B->mass*u2) / denom) + Apn);
		B->Velocity(((2.0f * A->mass*u1 + ((B->mass - A->mass)*u2)) / denom) + Bpn);
	}

}


#endif _COLLISIONS_H