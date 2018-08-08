
#include "GLIncludes.h"
#include "all_Fixed_Data.h"

//Similar code to renderBody(), however here we use GL_LINES instead of GL_TRIANGLES
void renderLines()
{
	glBufferData(GL_ARRAY_BUFFER, sizeof(VertexFormat) * lineBuffer.numberOfVertices, &lines[0], GL_STATIC_DRAW);


	glUniformMatrix4fv(uniMVP, 1, GL_FALSE, glm::value_ptr(MVP));
	glBindBuffer(GL_ARRAY_BUFFER, lineBuffer.vbo);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(VertexFormat), (void*)0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(VertexFormat), (void*)16);

	glDrawArrays(GL_LINES, 0, lineBuffer.numberOfVertices);



}

void drawRect()
{
	//Clearing buffer so we dont overflow it
	glDeleteBuffers(1, &lineBuffer.vbo);
	lines.clear();
	//rectangle box for Quadtree
	lines.push_back(VertexFormat(glm::vec3(-1, 1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(1, 1, 0), glm::vec4(1, 0, 0, 1)));

	lines.push_back(VertexFormat(glm::vec3(1, 1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(1, -1, 0), glm::vec4(1, 0, 0, 1)));

	lines.push_back(VertexFormat(glm::vec3(1, -1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(-1, -1, 0), glm::vec4(1, 0, 0, 1)));

	lines.push_back(VertexFormat(glm::vec3(-1, -1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(-1, 1, 0), glm::vec4(1, 0, 0, 1)));




}

struct bounds {

	float xMax = 1;
	float xMin = -1;
	float yMax = 1;
	float yMin = -1;
} base;


void drawNode(KD_Node node, bounds space)
{
	float x = node.Get_point().xpos;
	float y = node.Get_point().ypos;


	bounds leftSpace = space;
	bounds rightSpace = space;

	if (node.Get_type() == 100)
	{

		lines.push_back(VertexFormat(glm::vec3(x, space.yMax, 0), glm::vec4(1, 0, 0, 1)));
		lines.push_back(VertexFormat(glm::vec3(x, space.yMin, 0), glm::vec4(1, 0, 0, 1)));

		leftSpace.xMax = x;
		rightSpace.xMin = x;
	}
	else if (node.Get_type() == 101)
	{

		lines.push_back(VertexFormat(glm::vec3(space.xMin, y, 0), glm::vec4(1, 0, 0, 1)));
		lines.push_back(VertexFormat(glm::vec3(space.xMax, y, 0), glm::vec4(1, 0, 0, 1)));

		leftSpace.yMax = y;
		rightSpace.yMin = y;
	}




	if (!node.isleaf())
	{
		if (node.Get_leftTree() != nullptr)
			drawNode(*node.Get_leftTree(), leftSpace);

		if (node.Get_rightTree() != nullptr)
			drawNode(*node.Get_rightTree(), rightSpace);
	}

}


void drawKDTree(KD_Node root)
{
	drawRect();


	drawNode(root, base);
	lineBuffer.initBuffer(lines.size(), &lines[0]);
}
