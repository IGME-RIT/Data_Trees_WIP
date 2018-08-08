/*
File Name: main.cpp
Copyright © 2018
Original authors: Sanketh Bhat, Benjamin Evans
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Description:
space partition data structures

References:
CS 3250: Computational Geometry Spring 2017 by Laura Toma at Bowdoin College
Bella Tumaneng - https://github.com/itumanen/kd_tree
https://avdongre.wordpress.com/2011/06/14/kd-tree-in-c/
https://www.geeksforgeeks.org/k-dimensional-tree/
*/


#include "GLIncludes.h"
#include "GLRender.h"
#include "GameObject.h"
#include "Grid_Up.h"
#include "Model.h"
#include "Collisions.h"
#include "all_Fixed_Data.h"
#include "Rec_lines.h"

#include "Data_Tree_Maker.h"




using namespace std;

#pragma region util_functions

void setup()
{
	for (int i = 0; i < number_of_objects; i++)
	{
		bodies.push_back(new GameObject(mesh));
	}
}

void update()
{


	for each (GameObject* body in bodies)
	{

		tempVelocity = body->Velocity(); //real time change )
		tempPosition = body->Position();
		
		for each (GameObject* body2 in bodies)
		{
			if (CheckCollision(body, body2))
			{
				resolveCollisionimpulse(body, body2);
			}

		}
		
		if (tempPosition.x < min_x)
		{
			body->Velocity(glm::vec3(min_x_2 * tempVelocity.x, tempVelocity.y, tempVelocity.z));
		}

		if (tempPosition.y < min_y)
		{
			body->Velocity(glm::vec3(tempVelocity.x, min_y_2 * tempVelocity.y, tempVelocity.z));
		}
		if (tempPosition.x > max_x)
		{
			body->Velocity(glm::vec3(min_x_2 * tempVelocity.x, tempVelocity.y, tempVelocity.z));
		}
		if (tempPosition.y > max_y)
		{
			body->Velocity(glm::vec3(tempVelocity.x, min_y_2 * tempVelocity.y, tempVelocity.z));
		}
		
		body->update(timestep, PV);

	}
}




// This function is used to handle key inputs.
// It is a callback funciton. i.e. glfw takes the pointer to this function (via function pointer) and calls this function every time a key is pressed in the during event polling.
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{

	//Sets the current window to a close flag when ESC is pressed
	if (key == GLFW_KEY_ESCAPE && ((action == GLFW_PRESS) || action == GLFW_REPEAT))
	{
		glfwSetWindowShouldClose(window, 1);
	}

}
#pragma endregion

//Makes a 2D polygon based on the vertex number arguement
Model* setupModel(int n, glm::vec4 color, float size)
{
	std::vector<GLuint> indices;
	std::vector<VertexFormat> vertices;
	VertexFormat center(glm::vec3(0.0f, 0.0f, 0.0f), color);
	
	//Only add indices if you drawing a polygon with more than 3 points.
	if (n > 3)
	{
		//Indices are added in threes to form tris
		for (int i = 1; i < n+1; i++)
		{
			indices.push_back(0); //Start at the center
			if(i==n) //If we are at the last last vertex, go back to the first  non-center vertex and add it
			{
				indices.push_back(i);
				indices.push_back(1);

			}
			else
			{	//Adds current vertex and the next one
				indices.push_back(i);
				indices.push_back(i + 1);

			}

		}

		//Only 3+ point polygons need a center vertex
		vertices.push_back(center);
	}	

	float theta = 360.0f / n;
	

	//vertex generation
	for (int i = 0; i < n; i++)
	{
		//The point at angle theta  are fed into the buffer.
		vertices.push_back(VertexFormat(glm::vec3(size * cos(glm::radians(i*theta)), size * sin(glm::radians(i*theta)), 0.0f), color));
		
	}
	 return new Model(vertices.size(), vertices.data(), indices.size(), indices.data());

}

void main()
{
	// Initializes most things needed before the main loop
	init();

	// Sends the funtion as a funtion pointer along with the window to which it should be applied to.

	glfwSetKeyCallback(window, key_callback);


	//Creating model

	mesh = setupModel(20, glm::vec4(1.0f, 0.0f, 0.0f, 1.0f), 0.02f);

	//Sets up bodies in the scene

	setup();
	grid();
	coords_tree();
	//_2norm_dist();




	//printing KDtree
	KD_tree* tree;
	tree = new KD_tree(Points_Tree);
	tree->printInfo();
	tree->printTree();
	
	//making Quadtree
	std::vector<std::shared_ptr<QT_data_point>> QT_pointlist;
	quadtree_maker qt(QT_Data_Box{ -1.5f, -1.5f, 1.5f, 1.5f });
	for (int i = 0; i < number_of_objects; i++)
	{
		QT_pointlist.push_back(std::shared_ptr<QT_data_point>(new QT_data_point{ Points_Tree[i].xpos, Points_Tree[i].ypos }));
	}
	for (auto & p : QT_pointlist)
	{
		qt.place(p.get());
	}
	



	//testing
	for (int i = 0; i < number_of_objects; i++)
	{
		std::cout << " (" << Point_Coord[i].x1 << ", " << Point_Coord[i].y1 << ") " << "\n";
	}
	

	// Enter the main loop.
	while (!glfwWindowShouldClose(window))
	{
		
		// Call to checkTime() which will determine how to go about updating via a set physics timestep as well as calculating FPS.
		update();

		drawKDTree(*(tree->get_Root()));
		// Call the render function(s).
		renderScene();
		renderLines();

		for each (GameObject* body in bodies)
		{

			body->render(uniMVP);
		}

		// Swaps the back buffer to the front buffer
		// Remember, you're rendering to the back buffer, then once rendering is complete, you're moving the back buffer to the front so it can be displayed.
		glfwSwapBuffers(window);

		// Checks to see if any events are pending and then processes them.
		glfwPollEvents();

	

	}



	delete mesh;
	
	for each (GameObject* body in bodies)
	{
		delete(body);
	}

	//Cleans shaders and the program and frees up GLFW memory
	cleanup();

	return;
}

















/*
void gridSetup()
{
for (int i = 0; i < bodies.size(); i++)
{
float w = randposition() * 2 - 1;
float h = randposition() * 2 - 1;
bodies[i]->Position(glm::vec3(w, h, 0));


//float x = (randposition() * 2 - 1)*0.005f;
//float y = (randposition() * 2 - 1)*0.005f;
//bodies[i]->Velocity(glm::vec3(x, y, 0));

}
}
void K_DTreePos(GameObject* A)
{






}

*/


/*
glm::vec3 GetMousePos()
{

	double x, y;
	glfwGetCursorPos(window, &x, &y);

	//normalized device coordinates
	float NDCX = (x / 400.0f) - 1.0f;
	float NDCY = -(y / 400.0f) + 1.0f;

	return glm::vec3(NDCX, NDCY, 0.0f);


}
*/



//for each (GameObject* body in bodies)
//{
//glm::vec3 temp = glm::vec3(0.005, 0.005, 0); // change intial velocity for every object
//body->Velocity(temp);
//}




/*
float CheckDynamicCollision(GameObject* c1, GameObject* c2, const glm::vec3 &mvmt, float tStart, float tEnd)
{
	//Get the midpoint time of the interval
	float tMid = ((0.5f* (tEnd - tStart)) + tStart);
	//Determine the position at time tMid
	glm::vec3 pos = tMid * mvmt + c1->Position();

	float mvmtLength = (glm::length((tEnd - tStart) * mvmt));
	float a = (glm::distance(pos, c2->position));
	//Determine the radius of the bounding circle needed to encapsulate this movement
	float rad = mvmtLength + 0.0002f;

	//Create the bounding circle


	bool colCheck = ((rad + 0.00002) >= glm::distance(pos, c2->position));
	//check the bounding circle colliding with c2
	if (colCheck)
	{
		//Set the result as the furthest point along the current interval.
		//This ensures that when the function returns t will represent the position along the movement vector
		//Which will cause you to be entering collision.
		//
		//If you had wanted to be just before the collision, you could set it to tStart.
		float result = tMid;

		//Note: This if statement determines the accuracy of the algorithm!
		//The smaller you make the right side, the more accurate the algorithm will become!
		//However this is at a big cost of speed.
		if (mvmtLength > 0.00002f)
		{

			//Divide the interval in 2 and check the first half
			result = CheckDynamicCollision(c1, c2, mvmt, tStart, tMid);
			//If the result is false, check the second half
			if (result < 0.0f)
			{
				result = CheckDynamicCollision(c1, c2, mvmt, tMid, tEnd);
			}
		}
		//Return the result
		return result;
	}
	else
	{
		return -1.0f;
	}

}



void continuousCollisionUpdate(float timestep, GameObject* c1, GameObject* c2)
{

	//Determine relative velocity of circle 1 from circle 2
	glm::vec3 relV = ((c1->Velocity()) - (c2->Velocity()));

	float t = (CheckDynamicCollision(c1, c2, relV * timestep, 0.1f, 9.9f));
	if (t <= 0.0f)
	{
		//Reposition at point of intersection
		c1->Position() += c1->Velocity() * timestep * t;
		c2->Position() += c2->Velocity() * timestep * t;

	}
	else
	{
		//Move the circlebody
		c1->Position() += c1->Velocity() * timestep;
		c2->Position() += c2->Velocity() * timestep;

	}


	std::sort(domain.begin(), domain.end(), [](float a, float b) {
	return a > b;
	});
	for (auto a : domain) {
	std::cout << a << " ";
	}
	std::cout << '\n';



}

*/


//to get the Vertices for collision
//bodies[i]->Vertices();

/*
void gettingPoints()
{


//bodies[0];
//std::vector<glm::vec3>  a = (bodies[1]->Vertices());






for (auto i = bodies.begin(); i != bodies.end(); i++)
{


//std::vector<glm::vec3>  a = (bodies[i]->Vertices());
}
}

*/
