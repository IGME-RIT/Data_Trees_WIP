#ifndef _GRID_UP_H
#define _GRID_UP_H

#include "all_Fixed_Data.h"


float randposition()
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::minstd_rand gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<float> dis(-1.0f, 1.0f); //https://en.cppreference.com/w/cpp/numeric/random
	float posrand = dis(gen);
	return posrand;

	/*
	if std::uniform_real_distribution is not working right use this 
	return  static_cast <float> ((0.999995f - 0.001f) * rand()) / static_cast <float> (RAND_MAX + 0.001f);
	float w = randposition() * 2 - 1;
	float h = randposition() * 2 - 1;
	//float x = (randposition() * 2 - 1)*0.005f;
	//float y = (randposition() * 2 - 1)*0.005f;
	*/
}

void grid()
{

	for (int j = 0; j < bodies.size(); j++)
	{

		float w = randposition();
		float h = randposition();
		bodies[j]->Position(glm::vec3(w, h, 0));




		point_body_x[j] = bodies[j]->Position(glm::vec3(w, h, 0)).x;
		point_body_y[j] = bodies[j]->Position(glm::vec3(w, h, 0)).y;


		float x = randposition()*0.01f;
		float y = randposition()*0.01f;
		bodies[j]->Velocity(glm::vec3(0, 0, 0));



		
		//rounding the x positions of the body
		if (point_body_x[j] >= -1.0f  && point_body_x[j] <= -0.75f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else if (point_body_x[j] >= -0.75f  && point_body_x[j] <= -0.5f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else if (point_body_x[j] >= -0.5f  && point_body_x[j] <= -0.25f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else if (point_body_x[j] >= -0.25f  && point_body_x[j] <= 0.00000001f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else if (point_body_x[j] >= 0.00000001f  && point_body_x[j] <= 0.25f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else if (point_body_x[j] >= 0.25f  && point_body_x[j] <= 0.5f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else if (point_body_x[j] >= 0.5f  && point_body_x[j] <= 0.75f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else if (point_body_x[j] >= 0.75f  && point_body_x[j] <= 1.0f)
		{
			point_tree_x[j] = roundf(point_body_x[j] * 10) / 10;
		}
		else
		{
			point_tree_x[j] = 0;
		}

		//rounding the y positions of the body
		if (point_body_y[j] >= -1.0f  && point_body_y[j] <= -0.75f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else if (point_body_y[j] >= -0.75f  && point_body_y[j] <= -0.5f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else if (point_body_y[j] >= -0.5f  && point_body_y[j] <= -0.25f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else if (point_body_y[j] >= -0.25f  && point_body_y[j] <= 0.00000001f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else if (point_body_y[j] >= 0.00000001f  && point_body_y[j] <= 0.25f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else if (point_body_y[j] >= 0.25f  && point_body_y[j] <= 0.5f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else if (point_body_y[j] >= 0.5f  && point_body_y[j] <= 0.75f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else if (point_body_y[j] >= 0.75f  && point_body_y[j] <= 1.0f)
		{
			point_tree_y[j] = roundf(point_body_y[j] * 10) / 10;
		}
		else
		{
			point_tree_y[j] = 0;
		}
		


		
		//for every body positions
		//std::cout << " (" << point_body_x[j] << ", " << point_body_y[j] << ") " << "\n"; 
		//for the body that is position 0 in the array																				 
		//std::cout << " (" << point_body_x[0] << ", " << point_body_y[0] << ") " << "\n"; 
		//for every body positions that will be going into the data tree
		//std::cout << " (" << point_tree_x[j] << ", " << point_tree_y[j] << ") " << "\n"; 
		
		




	}
	
	//std::cout << " ("  << point_body_x[0] << ", " << point_body_y[0] << ") " << "\n";
}


#endif _GRID_UP_H