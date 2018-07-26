
#ifndef _DATA_TREE_NODE_H
#define _DATA_TREE_NODE_H



#include "all_Fixed_Data.h"
#include "GLIncludes.h"

using namespace std;

#pragma region coords_tree


void coords_tree()
{
	//restoring the x and y coords into an array
	for (int a = 0; a < number_of_objects; a++)
	{
		float _x = point_tree_x[a];
		float _y = point_tree_y[a];

		domain.push_back(_x);
		range.push_back(_y);

		point_temp.xpos = _x;
		point_temp.ypos = _y;
		//point_temp.id = a;

		Points_Tree.push_back(point_temp);

		Point_Coord.push_back(_Point_xy(_x, _y));

	}
	//PointMatrix;
	//Points_Tree;
	//Point_Coord;
}



#pragma endregion

#pragma region Kd tree


KD_Node::KD_Node()
{
	Tree_Point empty_point;
	empty_point.xpos = 0.0f;
	empty_point.ypos = 0.0f;
	set_Point(empty_point);
	set_Type(LEAF);
	set_Root(false);
	set_Depth(0);
	set_Left(NULL);
	set_Right(NULL);
}

KD_Node::KD_Node(Tree_Point t_point)
{
	set_Point(t_point);
	set_Type(LEAF);
	set_Root(false);
	set_Depth(1);
	set_Left(NULL);
	set_Right(NULL);
}



void KD_Node::print_Info()
{
	print_Point(Get_point());
	print_Type();
	if (!isRoot()) {
		printf("Depth: %d\n", Get_depth());
	}

}



void KD_Node::print_Type()
{
	string type;
	if (isRoot()) return; // tree will state if it's printing the root
	if (Get_type() == VERTICAL) {
		type = "VERTICAL";
	}
	else if (Get_type() == HORIZONTAL) {
		type = "HORIZONTAL";
	}
	else if (Get_type() == LEAF) {
		type = "LEAF";
	}
	printf("Type: %s\n", type.c_str());
}


KD_Node::~KD_Node() {

}

KD_tree::KD_tree()
{
	set_TreeRoot(NULL);
	set_TreeHeight(0);
	set_numKDNodes(0);

}

KD_tree::KD_tree(std::vector<Tree_Point> all_points)
{
	sort(all_points.begin(), all_points.end(), sort_X_coords);
	set_numKDNodes(all_points.size());
	set_Pts(all_points);
	init_Height();

	initial_tree(all_points);
}

void KD_tree::remaker_KD(std::vector<Tree_Point> all_points)
{
	pts_set.clear();
	ordered_pts_set.clear();
	deallocate_tree(get_Root());
	set_TreeRoot(NULL);

	

	sort(all_points.begin(), all_points.end(), sort_X_coords);
	set_numKDNodes(all_points.size());
	set_Pts(all_points);
	init_Height();

	initial_tree(all_points);
}

void KD_tree::initial_tree(std::vector<Tree_Point> all_points)
{

	if (num_kd_nodes == 1) //if there is only one coord in all_points
	{
		KD_Node * root = new KD_Node(pts_set[0]);
		root->set_Depth(get_treeHeight());
		root->set_Root(true);
		set_TreeRoot(root);


	}
	else if (num_kd_nodes == 2) //if there is two coords in all_points
	{
		KD_Node * root = new KD_Node(pts_set[1]);
		root->set_Depth(get_treeHeight()+1);
		root->set_Type(INITIAL_CUT);
		root->set_Root(true);
		set_TreeRoot(root);

		KD_Node * left_Tree = new KD_Node(pts_set[0]);
		root->set_Depth(get_treeHeight() + 1);
		root->set_Root(true);
		set_TreeRoot(root);

		set_TreeHeight(compute_KD_Height());

	}
	else
	{
		// Initialize root node
		int median_idx = num_kd_nodes / 2;

		KD_Node * root = new KD_Node(all_points[median_idx]);
		root->set_Type(INITIAL_CUT);
		root->set_Depth(get_treeHeight());
		root->set_Root(true);
		set_TreeRoot(root);


		// making vectors for sorted list of positions
		std::vector<Tree_Point> x_left_arr, x_right_arr, y_left_arr, y_right_arr;

		// initialize by copying data from pts_set
		// for loops done separately because left/right arrays could
		// be of different sizes 
		// if the number of points is even, the right array is one 
		// index smaller than the left
		// Example:
		// 0 1 2 3 4, median index 2 --> left [0,1] right [3,4]
		// 0 1 2 3, median 2 --> left [0,1] right [3]

		for (int i = 0; i < median_idx; i++)
		{
			x_left_arr.push_back(pts_set[i]);
			y_left_arr.push_back(pts_set[i]);
		}
		for (int i = median_idx + 1; i < num_kd_nodes; i++) 
		{
			x_right_arr.push_back(pts_set[i]);
			y_right_arr.push_back(pts_set[i]);
		}

		sort(x_left_arr.begin(), x_left_arr.end(), sort_X_coords);
		sort(x_right_arr.begin(), x_right_arr.end(), sort_X_coords);
		sort(y_left_arr.begin(), y_left_arr.end(), sort_Y_coords);
		sort(y_right_arr.begin(), y_right_arr.end(), sort_Y_coords);

		root->set_Left(KD_Maker(x_left_arr, y_left_arr, root->Get_depth() + 1));
		root->set_Right(KD_Maker(x_right_arr, y_right_arr, root->Get_depth() + 1));

		set_TreeHeight(compute_KD_Height());
	}

	pts_set.clear();
	levelOrderPts();
}

KD_Node* KD_tree::KD_Maker(std::vector<Tree_Point> X_point, std::vector<Tree_Point> Y_point, int depth)
{
	int num = X_point.size();
	int middle_position = num / 2;

	if (num == 1)
	{
		KD_Node* node = new KD_Node(X_point[0]);
		node->set_Depth(depth);
		node->set_Left(NULL);
		node->set_Right(NULL);

		set_TreeHeight(depth);
		return node;
	}

	else if (depth % 2 != 0 && num > 1) //height is odd at new node
	{
		KD_Node* odd_node = new KD_Node(X_point[middle_position]);
		odd_node->set_Type(VERTICAL);
		odd_node->set_Depth(depth);

		// making vectors for sorted list of odd node positions
		std::vector<Tree_Point> x_left_odd, x_right_odd, y_left_odd, y_right_odd;

		for (int i = 0; i < middle_position; i++)
		{
			x_left_odd.push_back(X_point[i]);
			y_left_odd.push_back(X_point[i]);
		}
		for (int i = middle_position + 1; i < num; i++)
		{
			x_right_odd.push_back(X_point[i]);
			y_right_odd.push_back(X_point[i]);
		}

		//sort y-coords
		sort(y_left_odd.begin(), y_left_odd.end(), sort_Y_coords);
		sort(y_right_odd.begin(), y_right_odd.end(), sort_Y_coords);

		odd_node->set_Left(KD_Maker(x_left_odd, y_left_odd, depth + 1));
		odd_node->set_Right(KD_Maker(x_right_odd, y_right_odd, depth + 1));

		set_TreeHeight(depth);
		return odd_node;
	}
	else if (depth % 2 == 0 && num > 1) //height is even at new node
	{
		KD_Node * even_node = new KD_Node(Y_point[middle_position]);
		even_node->set_Type(HORIZONTAL);
		even_node->set_Depth(depth);

		// making vectors for sorted list of even node positions
		std::vector<Tree_Point> x_left_even, x_right_even, y_left_even, y_right_even;
		x_right_even;
		for (int i = 0; i < middle_position; i++)
		{
			x_left_even.push_back(Y_point[i]);
			y_left_even.push_back(Y_point[i]);
		}
		for (int i = middle_position + 1; i < num; i++)
		{
			x_right_even.push_back(Y_point[i]);
			y_right_even.push_back(Y_point[i]);
		}

		//sort x-coords
		sort(x_left_even.begin(), x_left_even.end(), sort_X_coords);
		sort(x_right_even.begin(), x_right_even.end(), sort_X_coords);

		even_node->set_Left(KD_Maker(x_left_even, y_left_even, depth + 1));
		even_node->set_Right(KD_Maker(x_right_even, y_right_even, depth + 1));

		set_TreeHeight(depth);
		return even_node;
	}

	return NULL;
}

int KD_tree::compute_KD_Height()
{
	KD_Node* temp_node = root;
	int height = 1;

	while (!temp_node->isleaf())
	{
	height++;
	temp_node = temp_node->Get_leftTree();
	}
	return height;
}

void KD_tree::levelOrderPts()
{
	for (int depth = 1; depth < height_Tree + 1; depth++)
	{
		add_Level(get_Root(), depth);
	}
}

void KD_tree::add_Level(KD_Node* tree_node, int depth)
{
	if (!tree_node)
	{
		return;
	}
	if (tree_node->Get_depth() == depth)
	{
		ordered_pts_set.push_back(tree_node);
return;
	}
	add_Level(tree_node->Get_leftTree(), depth);
	add_Level(tree_node->Get_rightTree(), depth);

}

void KD_tree::printInfo()
{

	printf("TREE INFO\n");
	printNumNodes();
	printf("Height: %d\n", get_treeHeight());
	//get_Points().at(0);
	printf("Root: ");

}

void KD_tree::printTree()
{
	if (!get_Root())
	{
		printf("Root is NULL; tree is empty.\n");
	}
	printf("\n\nPRINTING TREE:\n");
	for (int i = 0; i < ordered_pts_set.size(); i++)
	{
		ordered_pts_set[i]->print_Info();
		printf("\n");
	}

}



void KD_tree::printNumNodes() {
	printf("Number of nodes: %d\n", get_num_KDNodes());
}


void KD_tree::deallocate_tree(KD_Node* tree_node)
{
	if (tree_node->isleaf() || !tree_node)
	{
		return;
	}
	deallocate_tree(tree_node->Get_leftTree());
	deallocate_tree(tree_node->Get_rightTree());
	delete(tree_node);
}

KD_tree::~KD_tree()
{
	deallocate_tree(get_Root());
}

#pragma endregion


#pragma region QuadTree
//testing stuff



bool operator==(const QT_PointGetter & p1, const QT_PointGetter & p2)
{
	return (p1.get_X() == p2.get_X() && p1.get_Y() == p2.get_Y());
}
bool operator!=(const QT_PointGetter & p1, const QT_PointGetter & p2) 
{
	return !(p1 == p2);
}

bool QT_PointGetter::equals(const QT_PointGetter & point) const
{
	return !(get_X() != point.get_X() || get_Y() != point.get_Y());
}

float QT_Data_Box::get_X1() const { return mX1; }
float QT_Data_Box::get_Y1() const { return mY1; }
float QT_Data_Box::get_X2() const { return mX2; }
float QT_Data_Box::get_Y2() const { return mY2; }

bool QT_Data_Box::intersects_p(const QT_PointGetter &point) const
{
	return !((point.get_X() < get_X1()) || (point.get_Y() < get_Y1()) || 
			(point.get_X() > get_X2()) || (point.get_Y() > get_Y2()));
}

bool QT_Data_Box::intersects_b(const QT_Data_Box &box) const
{
	return !((get_X1() > box.get_X2()) || (get_Y1() > box.get_Y2()) || 
			(get_X2() < box.get_X1()) || (get_Y2() < box.get_Y1()));
}

quadtree_maker::quadtree_maker(const QT_Data_Box& bounds)
	: node_NW{ nullptr }, node_NE{ nullptr },
node_SW{ nullptr }, node_SE{ nullptr }, t_bounds{ bounds }
{}

quadtree_maker::~quadtree_maker()
{
	if (node_NW != nullptr) 
	{
		delete node_NW;
		delete node_NE;
		delete node_SW;
		delete node_SE;

	}
}

quadtree_maker::quadtree_maker(const quadtree_maker& qtpoint) : quadtree_maker(qtpoint.t_bounds)
{
	this->operator=(qtpoint);
}

quadtree_maker & quadtree_maker::operator=(const quadtree_maker& qtpoint)
{
	t_bounds = qtpoint.t_bounds;
	if (node_NW != nullptr)
	{
		delete node_NW;
		delete node_NE;
		delete node_SW;
		delete node_SE;

		node_NW = nullptr;
		node_NE = nullptr;
		node_SW = nullptr;
		node_SE = nullptr;
	}
	if (qtpoint.node_NW != nullptr)
	{
		split_tree();

		*node_NW = *qtpoint.node_NW;
		*node_NE = *qtpoint.node_NE;
		*node_SW = *qtpoint.node_SW;
		*node_SE = *qtpoint.node_SE;
	}
	if (!qtpoint.pointData.empty())
	{
		pointData = qtpoint.pointData;
	}
	return *this;
}

void quadtree_maker::split_tree()
{
	if (node_NW == nullptr)
	{
		//upper bounds
		float X0 = t_bounds.get_X1();
		float Y0 = t_bounds.get_Y1();
		//lower bounds
		float X1 = t_bounds.get_X2();
		float Y1 = t_bounds.get_Y2();
		//the midpoints
		float x_mid = (X0 + X1) / 2.0f;
		float y_mid = (Y0 + Y1) / 2.0f;

		node_NW = new quadtree_maker(QT_Data_Box(X0, Y0, x_mid, y_mid));
		node_NE = new quadtree_maker(QT_Data_Box(x_mid, Y0, X1, y_mid));
		node_SW = new quadtree_maker(QT_Data_Box(X0, y_mid, x_mid, Y1));
		node_SE = new quadtree_maker(QT_Data_Box(x_mid, y_mid, X1, Y1));
		while (!pointData.empty())
		{
			if (node_NW->place(pointData.back()));
			else if (node_NE->place(pointData.back()));
			else if (node_SW->place(pointData.back()));
			else if (node_SE->place(pointData.back()));

			pointData.pop_back();

		}

	}

}

void quadtree_maker::join()
{
	if (node_NW != nullptr)
	{
		if ((node_NW->pointData.empty()) &&
			(node_NE->pointData.empty()) &&
			(node_SW->pointData.empty()) &&
			(node_SE->pointData.empty()) &&

			(node_NW->node_NW == nullptr) &&
			(node_NE->node_NW == nullptr) &&
			(node_SW->node_NW == nullptr) &&
			(node_SE->node_NW == nullptr)) 
		{

			delete node_NW;
			delete node_NE;
			delete node_SW;
			delete node_SE;

			node_NW = nullptr;
			node_NE = nullptr;
			node_SW = nullptr;
			node_SE = nullptr;
		}
	}
}

bool quadtree_maker::place(QT_PointGetter* point)
{
	if (t_bounds.intersects_p(*point))
	{
		if (node_NW != nullptr)
		{

			if (node_NW->place(point)) {}
			else if (node_NE->place(point)) {}
			else if (node_SW->place(point)) {}
			else if (node_SE->place(point)) {}
		}
		else
		{
			for (auto& e : pointData)
			{
				if (e == point) return false;
			}
			pointData.push_back(point);

			if (*pointData.front() != *point)
			{
				split_tree();
			}


		}
		return true;
	}
	return false;
}

bool quadtree_maker::remove(QT_PointGetter* point)
{
	if (t_bounds.intersects_p(*point))
	{
		if (node_NW != nullptr)
		{
			if (node_NW->remove(point) ||
				node_NE->remove(point) ||
				node_SW->remove(point) ||
				node_SE->remove(point))
			{
				join();
				return true;
			}

		}

		for (auto i = pointData.begin(); i != pointData.end(); i++)
		{
			if (*i == point)
			{
				pointData.erase(i);
				return true;
			}

		}

	}
	return false;
}

void quadtree_maker::query_p(const QT_data_point& point, qtpoint_vec& points) const
{
	if (t_bounds.intersects_p(*(QT_PointGetter*)&point))
	{
		if (node_NW != nullptr)
		{
			qtpoint_vec NW; node_NW->query_p(point, NW);
			qtpoint_vec NE; node_NE->query_p(point, NE);
			qtpoint_vec SW; node_SW->query_p(point, SW);
			qtpoint_vec SE; node_NW->query_p(point, SE);

			points.reserve(NW.size() + NE.size() + SW.size() + SE.size());
			points.insert(points.end(), NW.begin(), NW.end());
			points.insert(points.end(), NE.begin(), NE.end());
			points.insert(points.end(), SW.begin(), SW.end());
			points.insert(points.end(), SE.begin(), SE.end());
		}
		if (!pointData.empty())
		{
			if (*pointData.back() == *(QT_PointGetter*)&point)
			{
				points = pointData;
			}

		}

	}
}

void quadtree_maker::query_b(const QT_Data_Box& bounds, qtpoint_vec& points) const
{
	if (t_bounds.intersects_b(bounds))
	{
		if (node_NW != nullptr)
		{
			qtpoint_vec NW; node_NW->query_b(bounds, NW);
			qtpoint_vec NE; node_NE->query_b(bounds, NE);
			qtpoint_vec SW; node_SW->query_b(bounds, SW);
			qtpoint_vec SE; node_NW->query_b(bounds, SE);

			points.reserve(NW.size() + NE.size() + SW.size() + SE.size());
			points.insert(points.end(), NW.begin(), NW.end());
			points.insert(points.end(), NE.begin(), NE.end());
			points.insert(points.end(), SW.begin(), SW.end());
			points.insert(points.end(), SE.begin(), SE.end());
		}
		if (!pointData.empty())
		{
			if (bounds.intersects_p(*pointData.back()))
			{
				points = pointData;
			}

		}

	}

}

void quadtree_maker::pull_p(const QT_data_point& point, qtpoint_vec& points)
{
	if (t_bounds.intersects_p(*(QT_PointGetter*)&point))
	{
		if (node_NW != nullptr)
		{
			qtpoint_vec NW; node_NW->pull_p(point, NW);
			qtpoint_vec NE; node_NE->pull_p(point, NE);
			qtpoint_vec SW; node_SW->pull_p(point, SW);
			qtpoint_vec SE; node_NW->pull_p(point, SE);

			points.reserve(NW.size() + NE.size() + SW.size() + SE.size());
			points.insert(points.end(), NW.begin(), NW.end());
			points.insert(points.end(), NE.begin(), NE.end());
			points.insert(points.end(), SW.begin(), SW.end());
			points.insert(points.end(), SE.begin(), SE.end());

			join();
		}
		if (!pointData.empty())
		{
			if (*pointData.back() == *(QT_PointGetter*)&point)
			{
				points = pointData;
				pointData.clear();
			}

		}

	}
}

void quadtree_maker::pull_b(const QT_Data_Box& bounds, qtpoint_vec& points)
{
	if (t_bounds.intersects_b(bounds))
	{
		if (node_NW != nullptr)
		{
			qtpoint_vec NW; node_NW->pull_b(bounds, NW);
			qtpoint_vec NE; node_NE->pull_b(bounds, NE);
			qtpoint_vec SW; node_SW->pull_b(bounds, SW);
			qtpoint_vec SE; node_NW->pull_b(bounds, SE);

			points.reserve(NW.size() + NE.size() + SW.size() + SE.size());
			points.insert(points.end(), NW.begin(), NW.end());
			points.insert(points.end(), NE.begin(), NE.end());
			points.insert(points.end(), SW.begin(), SW.end());
			points.insert(points.end(), SE.begin(), SE.end());

			join();
		}
		if (!pointData.empty())
		{
			if (bounds.intersects_p(*pointData.back()))
			{
				points = pointData;
				pointData.clear();
			}

		}

	}
}

void quadtree_maker::clear()
{
	if (node_NW != nullptr)
	{
		delete node_NW;
		delete node_NE;
		delete node_SW;
		delete node_SE;

		node_NW = nullptr;
		node_NE = nullptr;
		node_SW = nullptr;
		node_SE = nullptr;
	}


}






















#pragma endregion
//&& item->mX <= rNode->getX() + rNode->getW() && item->mY <= rNode->getX() + rNode->getH())

/*
bound_BOX::bound_BOX()
{
float width = 0.0f;
float height = 0.0f;

}

bound_BOX::bound_BOX(const float &a, const float &b)
{
float set_Width(a);
float set_Height(b);
}


QUAD_Node::QUAD_Node()
{
_Point_xy empty_Quad;
bound_BOX empty_QuadBox;
empty_Quad.x1 = 0.0f;
empty_Quad.y1 = 0.0f;
empty_QuadBox.get_height();
empty_QuadBox.get_width();


}

QUAD_Node::QUAD_Node(const _Point_xy &tp, const bound_BOX &t_box)
{
set_Point(tp);
set_Box(t_box);


}

Quadtree_::Quadtree_()
{
QUAD_Node Qnodetemp;
set_node(Qnodetemp);
}

Quadtree_::Quadtree_(const QUAD_Node &Qnode, const int t_level)
{
set_node(Qnode);



}
*/

/*
MinMax::MinMax(std::vector<float> d, std::vector<float> r)
{
d = domain;
r = range;
}

float MinMax::get_MinMax_X_value(std::vector<float> &a)
{

//auto large_element = max_element(a.begin(), a.end());
//max_valueX = a[large_element];
//return max_valueX;
min_max_Values tempx;

for (auto temp : a)
{
if (max_valueX < temp)
{
max_valueX = temp;
}
}


MinMax* tempvalues;
tempvalues = new MinMax(domain, range);

float t_x = tempvalues->get_MinMax_X_value(domain);


tempx.max = max_valueX;
//tempx.min = min_valueX;


}
class MinMax
{
public:

MinMax(std::vector<float> d, std::vector<float> r);

float get_MinMax_X_value(std::vector<float> &a);

float get_MinMax_Y_value(std::vector<float> a);



private:
float max_valueY;
float min_valueY;

float max_valueX;
float min_valueX;

};

// The objects that we want stored in the quadtree
struct Node_Quad
{
Point_XY pos_Quad;
float data_Quad;
Node_Quad(Point_XY _pos, float _data)
{
pos_Quad = _pos;
data_Quad = _data;
}
Node_Quad()
{
data_Quad = 0;
}
};



class Quadtree {
public:
	Quadtree(float x, float y, float width, float height, int level, int maxLevel);
	~Quadtree();

	void            AddObject(GameObject *object);
	std::vector<GameObject*> GetObjectsAt(float x, float y);
	void            Clear();

private:
	float           x;
	float           y;
	float           width;
	float           height;
	int             level;
	int             maxLevel;
	std::vector<GameObject*> objects;

	Quadtree *      parent;
	Quadtree *      NW;
	Quadtree *      NE;
	Quadtree *      SW;
	Quadtree *      SE;

	bool            contains(Quadtree *child, GameObject *object);
};

class Quadtree {

	enum Node {
		NW = 0,
		NE,
		SW,
		SE,
		NodeCount
	};

public:
	Quadtree();

	Quadtree(float left, float right, float top, float down, unsigned int numObjectsToGrow = 3);

	~Quadtree();

	void				AddObject(GameObject *object);

	void				Clear();

	std::vector<GameObject*>			GetObjectsAt(float x, float y);

private:
	double				left;

	double				right;

	double				top;

	double				down;

	unsigned int			numObjectsToGrow;

	std::vector<GameObject*>			objects;

	Quadtree * 			nodes;

	bool				isLeaf;

	bool				contains(GameObject *object);

	bool				contains(float x, float y);

	void				createLeaves();

	void				moveObjectsToLeaves();

};

*/




/*
struct kdNode *newKDnode(Point_xy Point_Coord[])
{
struct kdNode* temp = new kdNode;

for (int i = 0; i < 2; i++)
{
temp->point[i] = Point_Coord[i];
}

temp->left = temp->right = NULL;
return temp;
};

const int k = 2;


// A method to create a node of K D tree
struct Node* newNode(float arr[])
{
struct Node* temp = new Node;

for (int i = 0; i<k; i++)
temp->point[i] = arr[i];

temp->left = temp->right = NULL;
return temp;
}

// Inserts a new node and returns root of modified tree
// The parameter depth is used to decide axis of comparison
Node *insertRec(Node *root, float point[], unsigned depth)
{
// Tree is empty?
if (root == NULL)
return newNode(point);

// Calculate current dimension (cd) of comparison
unsigned cd = depth % k;

// Compare the new point with root on current dimension 'cd'
// and decide the left or right subtree
if (point[cd] < (root->point[cd]))
root->left = insertRec(root->left, point, depth + 1);
else
root->right = insertRec(root->right, point, depth + 1);

return root;
}

// A utility method to determine if two Points are same
// in K Dimensional space
bool arePointsSame(float point1[], float point2[])
{
// Compare individual pointinate values
for (int i = 0; i < k; ++i)
if (point1[i] != point2[i])
return false;

return true;
}

// Searches a Point represented by "point[]" in the K D tree.
// The parameter depth is used to determine current axis.
bool searchRec(Node* root, float point[], unsigned depth)
{
// Base cases
if (root == NULL)
return false;
if (arePointsSame(root->point, point))
return true;

// Current dimension is computed using current depth and total
// dimensions (k)
unsigned cd = depth % k;

// Compare point with root with respect to cd (Current dimension)
if (point[cd] < root->point[cd])
return searchRec(root->left, point, depth + 1);

return searchRec(root->right, point, depth + 1);
}

// Searches a Point in the K D tree. It mainly uses
// searchRec()
bool search(Node* root, float point[])
{
// Pass current depth as 0
return searchRec(root, point, 0);
}


for (int a = 0; a < number_of_objects; a++)
{
tree_points_vec.push_back(tree_Coord(point_tree_x[a], point_tree_y[a]));
//point_xy[a] = tree_points_vec[a];

}


/*
struct kdNode {
Point_xy* point;
kdNode *left;
kdNode *right;

};


/*
int Find_Median_index(const std::vector<Tree_Point>& P)
{
if (P.size() % 2 == 0)
return ((P.size() / 2 - 1) + (P.size() / 2)) / 2; // even
else
return P.size() / 2; // odd
}

Tree_Point Median_Point(const std::vector<Tree_Point> P_T)
{
float temppos = Find_Median_index(P_T);
return P_T[temppos];
}


float fristCooord()
{
for (int i = 0; i < domain.size(); i++)
{
return domain[i];
}
}

float secCooord()
{
for (int i = 0; i < domain.size(); i++)
{
return range[i];
}
}

Point_xy get_P()
{
float u = fristCooord();
float v = secCooord();
return Point_xy(u, v);
}


*/








#endif _DATA_TREE_NODE_H