#ifndef _KNN_H
#define _KNN_H

#include "all_Fixed_Data.h"

float dx;
float dy;
float tempx;
float tempy;
/*
float Dist_calculate(const _Point_xy& P, const _Point_xy& Q)
{
	float dx = P.x - Q.x;
	float dy = P.y - Q.y;
	//return sqrt((dx*dx) + (dy*dy));
	return roundf(sqrt((dx*dx) + (dy*dy)) * 10) / 10;
}

void _2norm_dist()
{
	
	for (int di = 0; di < Point_Coord.size(); di++)
	{ 
		for (int dj = 0; dj<Point_Coord.size(); dj++)
		{ 
			

			dist_matrix[di][dj] = Dist_calculate(Point_Coord[di], Point_Coord[dj]);


			//std::cout << dist_matrix[di][dj] << ' ';
		}
		//std::cout << std::endl;
	}
	
	//dist_matrix;

}
*/





#endif _KNN_H