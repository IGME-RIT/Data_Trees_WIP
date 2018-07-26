
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

	//rectangle box for Quadtree
	lines.push_back(VertexFormat(glm::vec3(-1.1, 1.1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(1.1, 1.1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(1.1, 1.1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(1.1, -1.1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(1.1, -1.1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(-1.1, -1.1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(-1.1, -1.1, 0), glm::vec4(1, 0, 0, 1)));
	lines.push_back(VertexFormat(glm::vec3(-1.1, 1.1, 0), glm::vec4(1, 0, 0, 1)));


	lineBuffer.initBuffer(lines.size(), &lines[0]);

}
