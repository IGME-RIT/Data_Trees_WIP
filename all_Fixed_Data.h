#ifndef _ALL_FIXED_DATA_H
#define _ALL_FIXED_DATA_H


#include "GameObject.h"
#include "Model.h"
#include "GLIncludes.h"

#define number_of_objects 20 // Dont have this number be greater than 1000
#define EPSILON 0.00001

using namespace std;

#pragma region Coord_maker

struct Point_xy
{
	float x;
	float y;

	Point_xy(float u = 0, float v = 0) : x(u), y(v)
	{}
};
struct _Point_xy
{
	float x1;
	float y1;

	_Point_xy(float u = 0.0f, float v = 0.0f) : x1(u), y1(v)
	{}
};

typedef struct Tree_Point
{
	float xpos, ypos;

}Tree_Point;


#pragma endregion

#pragma region program specific Data members

float timestep = 0.016f;
float overlap;
bool collision = false;
int lastState;
const int k = 2;
float min_x = -1.0f;
float min_y = -1.0f;
float max_x = 1.0f;
float max_y = 1.0f;

float min_x_2 = -1.0f;
float min_y_2 = -1.0f;

glm::vec3 pointofcollision;
glm::vec3 minimumTranslationVector;
glm::vec3 mousePos;
glm::vec3 tempVelocity;
glm::vec3 tempPosition;
glm::vec3 tempPosition2;
glm::vec3 tempVelocity2;

float point_body_x[number_of_objects];
float point_body_y[number_of_objects];

float point_tree_x[number_of_objects];
float point_tree_y[number_of_objects];


// vector of scene bodies
std::vector<VertexFormat> lines;
std::vector<GameObject*> bodies;

//vectors for data-trees
std::vector<float> domain;
std::vector<float> range;
std::vector<Tree_Point> Points_Tree;
std::vector<_Point_xy> Point_Coord;





//norm distance vector
float dist_matrix[number_of_objects][number_of_objects];


VertexBuffer lineBuffer;
Model* mesh;
Tree_Point point_temp;


#pragma endregion

#pragma region KD_Tree and Node maker

enum cutType {
	VERTICAL = 100,
	HORIZONTAL,
	LEAF,
	INITIAL_CUT = VERTICAL
};

// Node in the k-d tree
class KD_Node
{
public:
	//constructor and destructor
	KD_Node();
	KD_Node(Tree_Point t_point);
	~KD_Node();

	bool isleaf() { return (type == LEAF); }
	bool isRoot() { return root; }

	//set and get functions
	Tree_Point Get_point() { return this->p_tree; }
	KD_Node* Get_leftTree() { return this->left_Tree; };
	KD_Node* Get_rightTree() { return this->right_Tree; };
	int Get_type() { return this->type; }
	int Get_depth() { return this->depth; }


	
	void set_Type(int T) { this->type = T; }
	void set_Root(bool bool_r) { this->root = bool_r; }
	void set_Depth(int D) { this->depth = D; }
	void set_Point(Tree_Point t_point) { this->p_tree = t_point; }
	void set_Left(KD_Node* L) { this->left_Tree = L; }
	void set_Right(KD_Node* R) { this->right_Tree = R; }
	

	
	// print functions
	void print_Info();
	void print_Type();
	void print_Point(Tree_Point p) { this->p1 = p; };
	
private:
	int type;
	bool root;
	int depth;
	Tree_Point p_tree;
	KD_Node* left_Tree;
	KD_Node* right_Tree;
	Tree_Point p1;
	

};

class KD_tree
{
public:
	//constructor and destructor
	KD_tree();
	KD_tree(std::vector<Tree_Point> all_points);
	~KD_tree();

	//helper functions
	void initial_tree(std::vector<Tree_Point> all_points);
	KD_Node* KD_Maker(std::vector<Tree_Point> X_point, std::vector<Tree_Point> Y_point, int depth);
	void remaker_KD(std::vector<Tree_Point> all_points);
	int compute_KD_Height();
	void levelOrderPts();
	void add_Level(KD_Node* node, int depth);
	void deallocate_tree(KD_Node* node);


	// check for even points
	bool Even_pts() { return (num_kd_nodes % 2 == 0); }

	//get functions
	int get_treeHeight() { return this->height_Tree; }
	int get_num_KDNodes() { return this->num_kd_nodes; }
	KD_Node * get_Root() { return this->root; }
	std::vector<KD_Node*> get_Points() { return this->ordered_pts_set; }

	//set functions
	void set_TreeRoot(KD_Node * r) { this->root = r; }
	void set_TreeHeight(int h) { this->height_Tree = h; }
	void init_Height() { this->height_Tree = 1; }
	void set_numKDNodes(int num) { this->num_kd_nodes = num; }
	void set_Pts(std::vector<Tree_Point> all_points) { this->pts_set = all_points; }

	//print functions
	void printTree();
	void printInfo();
	void printNumNodes();


private:
	KD_Node* root;
	int height_Tree;
	int num_kd_nodes;
	std::vector<Tree_Point> pts_set;
	std::vector<KD_Node*> ordered_pts_set;



	struct sort_X_coords
	{
		bool operator() (const Tree_Point &P, const Tree_Point &Q)
		{
			if (fabs(P.xpos - Q.xpos) < EPSILON)
			{
				return (P.ypos < Q.ypos);
			}
			else
			{
				return (P.xpos < Q.xpos);
			}

		}

	}sort_X_coords;

	struct _sort_Y_coords
	{
		bool operator() (const Tree_Point &P, const Tree_Point &Q)
		{
			if (fabs(P.ypos - Q.ypos) < EPSILON)
			{
				return (P.xpos < Q.xpos);
			}
			else
			{
				return (P.ypos < Q.ypos);
			}

		}

	}sort_Y_coords;
};

#pragma endregion

#pragma region Quad_Tree and Node maker
/*
struct QTNode
{
	_Point_xy pos;
	int data;

	QTNode(_Point_xy _pos, int _data)
	{
		pos = _pos;
		data = _data;
	}
	QTNode()
	{
		data = 0;
	}
};
*/



class QT_PointGetter
{
public:
	virtual float get_X() const = 0;
	virtual float get_Y() const = 0;
	virtual ~QT_PointGetter() = default;

	virtual bool equals(const QT_PointGetter & point) const;

};

class QT_data_point : public QT_PointGetter
{
private:
	float X;
	float Y;

public:
	//void set_X(const float& x) { X = x; }
	//void set_Y(const float& y) { Y = y; }

	QT_data_point::QT_data_point() : X(0.0f), Y(0.0f) {}
	QT_data_point::QT_data_point(const float& x, const float& y) : X(x), Y(y) {}

	float QT_data_point::get_X() const { return X; }
	float QT_data_point::get_Y() const { return Y; }
};



struct QT_Data_Box
{
private:
	float mX1;
	float mY1;
	float mX2;
	float mY2;

public:
	void set_X1(const float &x1) { mX1 = x1; }
	void set_Y1(const float &y1) { mY1 = y1; }
	void set_X2(const float &x2) { mX2 = x2; }
	void set_Y2(const float &y2) { mY2 = y2; }
	virtual float get_X1()const;
	virtual float get_Y1()const;
	virtual float get_X2()const;
	virtual float get_Y2()const;

	bool intersects_p(const QT_PointGetter& point) const;

	bool intersects_b(const QT_Data_Box &box) const;

	QT_Data_Box::QT_Data_Box() : mX1(0.0f), mY1(0.0f), mX2(0.0f), mY2(0.0f) {}

	QT_Data_Box::QT_Data_Box(const float& X1, const float& Y1, const float& X2, const float& Y2) :
		mX1(X1), mY1(Y1), mX2(X2), mY2(Y2) {}

}; 

typedef std::vector<QT_PointGetter*> qtpoint_vec;

class quadtree_maker
{
private:
	quadtree_maker* node_NW;
	quadtree_maker* node_NE;
	quadtree_maker* node_SW;
	quadtree_maker* node_SE;

	QT_Data_Box t_bounds;

	qtpoint_vec pointData;

	void split_tree();

	void join();

public:
	quadtree_maker(const QT_Data_Box& bounds);
	~quadtree_maker();

	quadtree_maker(const quadtree_maker& qtpoint);
	quadtree_maker& operator=(const quadtree_maker& qtpoint);

	bool place(QT_PointGetter* point);
	bool remove(QT_PointGetter* point);

	void query_b(const QT_Data_Box& bounds, qtpoint_vec& points) const;
	void query_p(const QT_data_point& point, qtpoint_vec& points) const;

	void pull_b(const QT_Data_Box& bounds, qtpoint_vec& points);
	void pull_p(const QT_data_point& point, qtpoint_vec& points);

	void clear();

	



};














#pragma endregion


#endif _ALL_FIXED_DATA_H











































/*




/*
struct AABB
{
_Point_xy min;
_Point_xy max;

AABB(const _Point_xy &minVal, const _Point_xy &maxVal)
{
min = minVal;
max = maxVal;
}
AABB()
{
min = _Point_xy(0.0f);
max = _Point_xy(0.0f);
}
};

*/


/*
struct AABB
{
Point_xy center;
Point_xy halfDIM;

AABB(Point_xy center = Point_xy(), Point_xy halfDIM = Point_xy()): center(center), halfDIM(halfDIM) {};
};

struct node_QUAD
{
Point_xy pos;
int idx;

node_QUAD(Point_xy pos = Point_xy(), int data = 0) : pos(pos), idx(data) {};
};


struct tree_Coord {

tree_Coord(float x, float y) :_x(x), _y(y)
{
}
/*
// Node in the k-d tree
struct K_D_Node {
Tree_Point coord;
K_D_Node* leftChild;
K_D_Node* rightChild;
};

float _x, _y;
};


struct POINT_XY {

float X, Y;

POINT_XY(float paramx, float paramy) : X(paramx), Y(paramy) {}

};

//vector<std::pair <float, float>> pair_array[number_of_objects];
//pair_array[a].push_back(std::make_pair(_x, _y));


std::vector<float> _point;



struct POINT_XY_picker {

POINT_XY XY_;

POINT_XY_picker(POINT_XY pointxy) : XY_(pointxy) {}

};




//For Testing Output look
/*
for (int b = 0; b < Point_XY[a].size(); b++)
{
std::cout << Point_XY[a][b] << " ";
}
cout << endl;

//std::vector<float> Point_at;
//std::vector<std::vector<float>> Point_XY;
//Point_at = { _x,_y };
//Point_XY.push_back(Point_at);

*/





/*
/*
struct POINT_XY {

float X, Y;


POINT_XY(float paramx, float paramy) : X(paramx), Y(paramy) {}
std::vector<POINT_XY> XY;
};

void coords_tree()
{
//restoring the x and y coords into an array
for (int a = 0; a < number_of_objects; a++)
{
XY.push_back(POINT_XY(point_tree_x[a] , point_tree_y[a]));

}
//std::array<float, 2> Point_at;
//std::vector<std::array<float, 2>> Point_XY;
//Point_at = { _x,_y };
//Point_XY.push_back(Point_at);
}
/*
//restoring the x and y coords into an array
for (int a = 0; a < number_of_objects; a++)
{
XY.push_back(POINT_XY(point_tree_x[a], point_tree_y[a]));

}
*/

/*
int s = XY.size();
for (int i = 0; i<s; i++)
{
// Accessing structure members using their
// names.
cout << XY[i].X << ", " << XY[i].Y
<< ", " << endl;
}


float Points[number_of_objects][2];

for (int i = 0; i < number_of_objects; ++i)
{
float _x = point_tree_x[i];
float _y = point_tree_y[i];
Points[i][0] = _x;
Points[i][1] = _y;
//for testing
//std::cout << " (" << Points[i][0] << ", " << Points[i][1] << ") " << "\t" << " (" << point_tree_x[i] << ", " << point_tree_y[i] << ") " << "\n";


}


*/
//cout << XY[0].X << ", " << XY[0].Y << endl;

/*
enum Quadrant {
NW, NE, SW, SE
};


class bound_BOX
{
public:
bound_BOX();
bound_BOX( const float &a, const float &b);

float get_width() { return this->width; }
float get_height() { return this->height; }


void set_Height(float h) { this->height = h; }
void set_Width(float w) { this->width = w; }


private:

float width;
float height;


};




class QUAD_Node
{
public:
//constructor
QUAD_Node();
QUAD_Node(const _Point_xy &tp, const bound_BOX &t_box);

_Point_xy Get_point() { return this->p_Quad; }
bound_BOX Get_box() { return this->t_bounds; }

void set_Point(_Point_xy tp) { this->p_Quad = tp; }
void set_Box(bound_BOX t_box) { this->t_bounds = t_box; }





private:
bound_BOX t_bounds;
_Point_xy p_Quad;


};



class Quadtree_
{
public:
Quadtree_();
Quadtree_(const QUAD_Node &Qnode, const int t_level);

QUAD_Node get_node() { return this->tempnode; }

void set_node(QUAD_Node Qnode) { this->tempnode = Qnode; }

private:
QUAD_Node tempnode;
//std::vector<int> node_idx;
int t_level;
std::unique_ptr<Quadtree_> t_subnodes[4];

};




*/
