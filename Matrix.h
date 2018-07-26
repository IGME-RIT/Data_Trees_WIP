#ifndef _MATRIX_H
#define _MATRIX_H
/*
class matrix {
	std::valarray<float> d;
	float rows, cols;
public:
	matrix(float r, float c) : d(0.0f, r*c), rows(r), cols(c) {}
	float& operator()(float r, float c) { return d[r*cols + c]; }
	float operator()(float r, float c) const { return d[r*cols + c]; }
	float size1() const { return rows; }
	float size2() const { return cols; }
};
matrix make_one(float r, float c) { return { r,c }; }
*/

/*
int main() {
auto m = make_one(5, 3);
m(2, 1) = 7; // demo storing a value
for (int r = 0; r != m.size1(); ++r) { // demo output
for (int c = 0; c != m.size2(); ++c) std::cout << m(r, c) << ' ';
std::cout << '\n';
}
}
*/

#endif _MATRIX_H
