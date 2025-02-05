#ifndef VEC_HPP
#define VEC_HPP

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <initializer_list>
#include <string>
#include <strstream>
#include <fstream>


using std::cout, std::endl, std::array, std::vector, std::initializer_list, std::string, std::strstream;
static long long int num_calls = 0;
const double PI = 3.141592653589793;

//enum class Axis { x, y, z };

template <typename T>
struct Vec { 
    T* data; // TODO don't use dynamic memory
    unsigned int length;
    Vec();
    Vec(int _length);
    Vec(std::initializer_list<T> args);
    Vec(const Vec&);
    Vec& operator=(const Vec& a);
    ~Vec();
    Vec slice(int start, int count) const;
    Vec resize(int new_length) const;
    Vec operator+(const Vec& a) const;
    Vec operator-(const Vec& a) const;
    Vec operator*(const T c) const;
    Vec operator/(const T c) const;
    bool operator==(const Vec& a) const;
    bool operator!=(const Vec& a) const;
    T mag() const;
    Vec unit() const;
    T dot(const Vec& a) const;
    Vec cross(const Vec& a) const;
    Vec lerp(const Vec& a, T amount);
};

struct Color {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t a;
    Color operator*(const Color c) const;
    Color operator*(const double k) const; 
};

template <typename T>
struct Tri {
    std::array<Vec<T>, 3> verts;
    Color color;
    bool draw_wire_frame;
    bool draw_face;
    bool draw_back_face;
    Tri();
    Tri(const Vec<T>& a, const Vec<T>& b, const Vec<T>& c, Color color,
            bool draw_wire_frame = true, bool draw_face = true, bool draw_back_face = false);
    Tri(const vector<Vec<T>>& verts_input, int start_index, Color color, 
            bool draw_wire_frame = true, bool draw_face = true, bool draw_back_face = false);
    Tri(std::initializer_list<Vec<T>> verts_input, Color color, 
            bool draw_wire_frame = true, bool draw_face = true, bool draw_back_face = false);
};

template <typename T>
struct Mat {
    T* data;
    unsigned int num_rows;
    unsigned int num_cols;

    Mat();
    Mat(int _num_rows, int _num_cols);
    Mat(std::initializer_list<Vec<T>> args);
    Mat(const Mat& _copy);
    Mat& operator=(const Mat& a);
    ~Mat();
    T get(int row, int col) const;
    void set(int row, int col, T value);
    void set_identity();
    Vec<T> get_row(int row) const;
    void set_row(int row, const Vec<T>& vec);
    void zero_row(int row);
    void add_to_row(int row, const Vec<T>& vec);
    //Mat slice(int 
    Mat operator*(const T a) const;
    Mat operator*(const Mat& a) const;
    Mat transpose() const;
    Mat sub_matrix(int _row, int _col) const;
    T minor(int row, int col) const;
    T cofactor(int row, int col) const;
    T det() const;  // determinant
    Mat adj() const;  // adjoint
    Mat inv() const;  // inverse
};

template <typename T>
struct Model {
    std::vector<Tri<T>> tris;
    Mat<T> trans;
    Mat<T> rot;
    Mat<T> scale;
    Color color;  // set all triangles to this color
    bool draw_wire_frame;
    bool draw_face;
    bool draw_back_face;
    Model();
    Model(std::initializer_list<Vec<T>> tris, Color color, bool draw_wire_frame = true, bool draw_back_face = false);
    bool load_from_file(string filename, Color color = {20, 20, 20, 20}, bool draw_wire_frame = true, 
            bool draw_back_face = false);
};

template <typename T>
struct PlaneLineSegmentIntersectResult {
    bool intersect;
    Vec<T> point;
};
//Vec<T> plane_line_segment_intersect(const Vec<T>& plane, const Vec<T>& p0, const Vec<T>& p1)

//// assume Vec<T> length is 3 or 4, but the number
//// of Vec<T>'s is always 3
//template <typename T>
//struct Triangle {
//    std::array<Vec<T>, 3> p;
//    Triangle();
//    ~Triangle();
//    Triangle(const Triangle&);
//    Triangle& operator=(const Triangle&);
//};
//
//template <typename T>
//struct Mesh {
//    std::vector<Triangle<T>> tris;
//    Mesh();
//    ~Mesh();
//    Mesh(const Mesh&);
//    Mesh& operator=(const Mesh&);
//};

template<typename T>
Vec<T> operator*(const T c, const Vec<T> a);

template<typename T>
Vec<T> operator/(const T c, const Vec<T> rhs);

// vector * matrix
template<typename T>
Vec<T> operator*(const Vec<T>& a, const Mat<T>& b);

// print vector
template<typename T>
std::ostream& operator<<(std::ostream& os, const Vec<T>& obj);

// print matrix
template<typename T>
std::ostream& operator<<(std::ostream& os, const Mat<T>& obj);

//struct Line {
//    double a, b, c;  // ax + by + c = 0
//    Line(double _a, double _b, double _c) : a(_a), b(_b), c(_c) {};
//};

// return line (in homogenous coords, with 3 dimension) that 
// passes through 2 points
template<typename T>
Vec<T> get_line_from_segment(const Vec<T>& pt0, const Vec<T>& pt1);

// line to line intersection. Lines expressed as homogenous coordinates
// with form ax + by + c = 0. element 0 = a, element 1 = b, etc
template<typename T>
Vec<T> line_line_intersection(const Vec<T>& a, const Vec<T>& b);



//
// definitions
//

template <typename T>
Vec<T>::Vec() 
{
    length = 0;
    data = nullptr;
}

// constructor
template <typename T>
Vec<T>::Vec(int _length) : length(_length) 
{ 
    assert(length > 0);
    data = new T[length];
    // probably unecessary?
    for (int i = 0; i < length; ++i) {
        data[i] = 0; 
    }
}

template <typename T>
Vec<T>::Vec(std::initializer_list<T> args) : length(args.size())
{
    int i = 0;
    data = new T[length];
    for(auto element : args) {
        data[i] = element; 
        i++;
    }
}

// copy constructor
template <typename T>
Vec<T>::Vec(const Vec& a) : length(a.length) 
{ 
    //assert(length > 0);
    data = new T[length];
    for (int i = 0; i < length; ++i) {
        data[i] = a.data[i]; 
    }
}

// copy assignment operator
template <typename T>
Vec<T>& Vec<T>::operator=(const Vec& a)
{
    if (this == &a) {
        return *this;
    }
    delete[] data;
    length = a.length;
    data = new T[length];
    for (int i = 0; i < length; ++i) {
        data[i] = a.data[i]; 
    }
    return *this;
}

// destructor
template <typename T>
Vec<T>::~Vec()
{
    delete[] data;
}

// return sub vector of length = count
template <typename T>
Vec<T> Vec<T>::slice(int start, int count) const
{
    assert(count > 0);
    assert(start + count <= length);
    Vec r(count);
    for (int i = start; i < start + count; ++i) { 
        r.data[i - start] = data[i]; 
    }
    return r;    
}

// returns resized vector. resizing begins at end of vector
template <typename T>
Vec<T> Vec<T>::resize(int new_length) const
{
    assert(new_length > 0);
    Vec r(new_length); 
    for (int i = 0; i < r.length; ++i) {
        r.data[i] = data[i];
    }
    return r;
}

template <typename T> 
Vec<T> Vec<T>::operator+(const Vec& a) const
{
    assert(length == a.length);
    Vec result(length);
    for (int i = 0; i < length; ++i) {
        result.data[i] = data[i] + a.data[i];
    }
    return result;
}

template <typename T> 
Vec<T> Vec<T>::operator-(const Vec& a) const 
{
    assert(length == a.length);
    Vec result(length);
    for (int i = 0; i < length; ++i) {
        result.data[i] = data[i] - a.data[i];
    }
    return result;
}

template <typename T> 
Vec<T> Vec<T>::operator*(const T c) const 
{
    Vec result(length);
    for (int i = 0; i < length; ++i) {
        result.data[i] = c * data[i];
    }
    return result;
}

template <typename T> 
Vec<T> Vec<T>::operator/(const T c) const 
{
    Vec result(length);
    for (int i = 0; i < length; ++i) {
        result.data[i] = data[i] / c;
    }
    return result;
}

template<typename T>
bool Vec<T>::operator==(const Vec& rhs) const
{
    assert(length == rhs.length);
    for(int i = 0; i < length; ++i) {
        if (data[i] != rhs.data[i]) {
            return false;
        }
    }
    return true;
}

template<typename T>
bool Vec<T>::operator!=(const Vec& rhs) const
{
    //return !(this* == rhs);
    assert(length == rhs.length);
    for(int i = 0; i < length; ++i) {
        if (data[i] != rhs.data[i]) {
            return true;
        }
    }
    return false;
}

template<typename T>
T Vec<T>::mag() const
{
    T temp = 0;
    for (int i = 0; i < length; ++i) {
        temp += data[i] * data[i];
    }
    return sqrt(temp);
}

template<typename T>
Vec<T> Vec<T>::unit() const
{
    Vec r(length);
    T m = this->mag();
    for (int i = 0; i < length; ++i) {
        r.data[i] = data[i] / m; 
    }
    return r;
}


template<typename T>
T Vec<T>::dot(const Vec& a) const
{
    assert(length == a.length);
    T r = 0;
    for (int i = 0; i < length; ++i) {
        r += data[i] * a.data[i];
    }
    return r;
}

class IntRotate {
    int num;
    int rotate_at;
public:
    IntRotate(int _num, int _rotate_at) : num(_num), rotate_at(_rotate_at) {
        assert(rotate_at > 0);
    }
    int inc() {
        num += 1;
        if (num >= rotate_at) {
            num = 0;
        }
        return num;
    }
    int dec() {
        num -= 1;
        if (num < 0) {
            num = rotate_at - 1;
        }
        return num;
    }
    int get() { return num; }
    void set(int _arg) { 
        assert( _arg >= 0 && _arg < rotate_at);
        num = _arg;
    }
};

// cross that accomodated/ignored 4th element
//template<typename T>
//Vec<T> Vec<T>::cross(const Vec& a) const
//{
//    //assert(length == 3 && length == a.length);
//    assert( (length == 3 || length == 4) && length == a.length);
//    Vec<T> i(length), j(length), k(length);
//    //i.x = 1.0f; j.y = 1.0f; k.z = 1.0f;
//    i.data[0] = 1.0f; j.data[1] = 1.0f; k.data[2] = 1.0f;
//    Vec<T> r1(length), r2(length);
//    //r1 = (y * a.z * i) + (z * a.x * j) + (x * a.y * k);
//    r1 = (data[1] * a.data[2] * i) 
//            + (data[2] * a.data[0] * j) 
//            + (data[0] * a.data[1] * k);
//    //cout << "r1 " << r1 << endl;
//    //r2 = -1.0 * (a.y * z * i) - (a.z * x * j) - (a.x * y * k);
//    r2 = -1.0 * (a.data[1] * data[2] * i) 
//            - (a.data[2] * data[0] * j) 
//            - (a.data[0] * data[1] * k);
//    //cout << "r2 " << r2 << endl;
//    //r = (y * a.z * i) + (z * a.x * j) + (x * a.y * k)
//    //        - (a.y * z * i) - (a.z * x * j) - (a.x * y * k);
//    auto r = r1 + r2;
//    //if (length == 4) {
//    //    r[3] = 1.0
//    return r1 + r2; 
//}

template<typename T>
Vec<T> Vec<T>::cross(const Vec& a) const
{
    //assert(length == 3 && length == a.length);
    assert( length == 3 && a.length == 3);
    Vec<T> i(3), j(3), k(3);
    //i.x = 1.0f; j.y = 1.0f; k.z = 1.0f;
    i.data[0] = 1.0f; j.data[1] = 1.0f; k.data[2] = 1.0f;
    Vec<T> r1(3), r2(3);
    //r1 = (y * a.z * i) + (z * a.x * j) + (x * a.y * k);
    r1 = (data[1] * a.data[2] * i) 
            + (data[2] * a.data[0] * j) 
            + (data[0] * a.data[1] * k);
    //cout << "r1 " << r1 << endl;
    //r2 = -1.0 * (a.y * z * i) - (a.z * x * j) - (a.x * y * k);
    r2 = -1.0 * (a.data[1] * data[2] * i) 
            - (a.data[2] * data[0] * j) 
            - (a.data[0] * data[1] * k);
    //cout << "r2 " << r2 << endl;
    //r = (y * a.z * i) + (z * a.x * j) + (x * a.y * k)
    //        - (a.y * z * i) - (a.z * x * j) - (a.x * y * k);
    return r1 + r2; 
}


template <typename T>    
Vec<T> Vec<T>::lerp(const Vec& a, T amount)
{
    assert(amount <= 1.0);
    assert(amount >= 0.0);
    return (1.0 - amount) * (*this) + amount * a; 
}

// TODO adjust alpha?
Color Color::operator*(const Color c) const {
    return {uint8_t(r * c.r), uint8_t(g * c.g), uint8_t(b * c.b), uint8_t(a) };
}

// TODO adjust alpha?
Color Color::operator*(const double k) const {
    return {uint8_t(r * k), uint8_t(g * k), uint8_t(b * k), uint8_t(a) };
}

template <typename T>
Tri<T>::Tri()
{
    verts = { Vec<T>(4), Vec<T>(4), Vec<T>(4)};
    //verts[0] = a;
    //verts[1] = b;
    //verts[2] = c;
    color = {100, 100, 100, 100};;
    draw_wire_frame = true;
    draw_face = true;
    draw_back_face = false;
}

template <typename T>
Tri<T>::Tri(
        const Vec<T>& a, 
        const Vec<T>& b, 
        const Vec<T>& c, 
        Color color_arg,
        bool draw_wire_frame_arg,
        bool draw_face_arg,
        bool draw_back_face_arg)
{
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
    color = color_arg;
    draw_wire_frame = draw_wire_frame_arg;
    draw_face = draw_face_arg;
    draw_back_face = draw_back_face_arg;
}

template<typename T>
Tri<T>::Tri(const vector<Vec<T>>& verts_input, 
        int start_index, 
        Color color_arg, 
        bool draw_wire_frame_arg,
        bool draw_face_arg,
        bool draw_back_face_arg)
{
    verts[0] = verts_input[start_index];
    verts[1] = verts_input[start_index + 1];
    verts[2] = verts_input[start_index + 2];
    color = color_arg;
    draw_wire_frame = draw_wire_frame_arg;
    draw_face = draw_face_arg;
    draw_back_face = draw_back_face_arg;
}

template <typename T>
Tri<T>::Tri(std::initializer_list<Vec<T>> verts_input, 
        Color color_arg, 
        bool draw_wire_frame_arg,
        bool draw_face_arg,
        bool draw_back_face_arg) 
{
    // expect 3 vertices 
    assert(verts_input.size() == 3);    
    int i = 0;
    for (const auto& vert : verts_input) {
        verts[i] = vert;
        i++;
    }
    color = color_arg;
    draw_wire_frame = draw_wire_frame_arg;
    draw_face = draw_face_arg;
    draw_back_face = draw_back_face_arg;
}

//template<typename T>
//Tri<T>::Tri(std::initializer_list<

// constructor

template<typename T>
Mat<T>::Mat() {
    // just make a 1x1 default matrix
    num_rows = 1;
    num_cols = 1;
    data = new T[1];
    data[0] = 1.0; // TODO or should be T{} ?
}

template<typename T>
Mat<T>::Mat(int _num_rows, int _num_cols) : num_rows(_num_rows), 
        num_cols(_num_cols)
{
    // each row vector is already zero, right?
    assert(num_rows > 0 && num_cols > 0);
    data = new T[num_rows * num_cols];
    for (int i = 0; i < num_rows * num_cols; ++i) {
        data[i] = 0;  //T(); 
    }
}

//Mat(std::initializer_list<std::initializer_list<Vec<T>>> args);
template<typename T>
Mat<T>::Mat(std::initializer_list<Vec<T>> args) {
    //int i = 0;
    //data = new T[length];
    //for(auto element : args) {
    //    data[i] = element; 
    //    i++;
    //}
    //cout << "Mat constructor with initializer_list" << endl;
    assert( args.size() != 0);
    for (auto& vec : args) {
        // get length of the vector
        num_cols = vec.length;     
        break;
    }
    num_rows = args.size();
    data = new T[num_rows * num_cols];
    int i = 0;
    for(auto& vec : args) {
        set_row(i, vec);
        i++;
    }
}

// copy constructor
template<typename T>
Mat<T>::Mat(const Mat& a)
{
    num_rows = a.num_rows; 
    num_cols = a.num_cols; 
    data = new T[num_rows * num_cols];
    for (int i = 0; i < num_rows * num_cols; ++i) {
        data[i] = a.data[i];
    }
}

// copy assignment
template<typename T>
Mat<T>& Mat<T>::operator=(const Mat& a)
{
    if (this == &a) { return *this; }
    delete[] data;
    num_rows = a.num_rows;
    num_cols = a.num_cols;
    data = new T[num_rows * num_cols];
    for (int i = 0; i < num_rows * num_cols; ++i) {
        data[i] = a.data[i];
    }
    return *this;
}

// destructor
template<typename T>
Mat<T>::~Mat()
{
    delete[] data;
}

template<typename T>
T Mat<T>::get(int row, int col) const
{
    assert( row >= 0 && row < num_rows 
        && col >= 0 && col < num_cols);
    //return rows[row].data[col];
    return data[row * num_cols + col]; 
}

template<typename T>
void Mat<T>::set(int row, int col, T value)
{
    assert( row >= 0 && row < num_rows 
        && col >= 0 && col < num_cols);
    //rows[row].data[col] = value;
    data[row * num_cols + col] = value;
}

template<typename T>
void Mat<T>::set_identity()
{
    assert(num_rows == num_cols);
    // zero everything first
    for (int i = 0; i < num_rows * num_cols; ++i) {
        data[i] = 0;
    }
    for (int i = 0; i < num_rows; ++i) {
        set(i, i, 1.0);
    }
}

// TODO test
template<typename T>
Vec<T> Mat<T>::get_row(int row) const {
    Vec<T> result(num_cols);
    for (int i = 0; i < num_cols; ++i) {
        result.data[i] = get(row, i);
    }
    return result;
}

// assign row vector to matrix
template<typename T>
void Mat<T>::set_row(int row, const Vec<T>& v)
{
    assert(v.length == num_cols); 
    for (int i = 0; i < v.length; ++i) {
        set(row, i, v.data[i]);
    }
}

template<typename T>
void Mat<T>::zero_row(int row) {
    for (int i = 0; i < num_cols; ++i) {
        set(row, i, T(0.0));
    }
}

// add vector to row in matrix
template<typename T>
void Mat<T>::add_to_row(int row, const Vec<T>& v)
{
    assert(v.length == num_cols); 
    for (int i = 0; i < v.length; ++i) {
        set(row, i, get(row, i) + v.data[i]);
    }
}

// scalar * matrix
template<typename T>
Mat<T> Mat<T>::operator*(const T b) const
{
    Mat<T> r(num_rows, num_cols);
    //for (int row = 0; row < num_cols; ++row) {
    //    for (int col = 0; col < num_rows; ++col) {
    //        r.set(row, col, sum);
    //    }
    //}
    for (int i = 0; i < num_rows * num_cols; ++i) {
        r.data[i] = b * data[i];
    }
    return r;
}

// matrix * matrix
template<typename T>
Mat<T> Mat<T>::operator*(const Mat<T>& b) const
{
    
    if(this->num_cols != b.num_rows) {
        cout << "CRASH! " << this->num_cols << ", " << b.num_rows << '\n';
    }

    assert(this->num_cols == b.num_rows);
    Mat<T> r(num_rows, b.num_cols);
    T sum = 0;
    for (int row = 0; row < this->num_cols; ++row) {
        for (int col = 0; col < b.num_rows; ++col) {
            sum = 0;
            for (int k = 0; k < this->num_cols; ++k) {
                sum += this->get(row, k) * b.get(k, col);
            }
            r.set(row, col, sum);
        }
    }
    return r;
}

template<typename T>
Mat<T> Mat<T>::transpose() const
{
    // num_rows, num_cols
    // m.rows , v.data
    Mat<T> result(num_cols, num_rows);
    // iterate through the elements of input matrix
    for (int row = 0; row < num_rows; ++row) {
        for (int col = 0; col < num_cols; ++col) {
            T input_element = get(row, col);
            result.set(col, row, input_element);
        }
    }
    return result;
}

// deletes a row and entry from input matrix and returns new matrix
template<typename T>
Mat<T> Mat<T>::sub_matrix(int row, int col) const 
{
    Mat<T> result(num_rows - 1, num_cols - 1);
   
    int result_index = 0;
    for (int i = 0; i < num_rows * num_cols; ++i) {
        // check if row should be skipped
        if (i >= row * num_cols && i < (row +1) * num_cols) {
            continue; // skip this row
        }
        // check if column should be skipped
        if (i % num_cols == col) {
            continue;
        }
        result.data[result_index] = data[i];
        result_index++;
    }
    return result; 
}

template<typename T>
T Mat<T>::minor(int row, int col) const 
{
    Mat<T> mat_new = sub_matrix(row, col);
    return mat_new.det();
    //return sub_matrix(row, col).det();
}

template<typename T>
T Mat<T>::cofactor(int row, int col) const
{
    // check if even or odd
    T c = ((row + col) % 2 == 0) ? T(1.0f) : T(-1.0f);
    return c * minor(row, col);
}

// determinant of matrix
template<typename T>
T Mat<T>::det() const
{
    // too much recursion? slow?
    assert(num_rows == num_cols);
    if (num_rows == 2) {
        return get(0,0) * get(1,1) - get(0, 1) * get(1, 0);
    }
    // choose top row when doing determinant calculation
    T sum = 0;
    for (int i = 0; i < num_cols; ++i) {
        sum += get(0, i) * cofactor(0, i);
    }
    return sum;
}

// adjoint matrix
template<typename T>
Mat<T> Mat<T>::adj() const  
{
    // classical adjoint method 
    Mat result(num_rows, num_cols);
    int cur_row = 0;
    int cur_col = 0;
    for (int i = 0; i < num_rows * num_cols; ++i) {
        // find indexes of matrix we are iterating thru (0 based indexing)
        cur_col = i % num_cols;
        cur_row = i / num_cols;
        result.set(cur_row, cur_col, cofactor(cur_row, cur_col));
    }
    return result.transpose();
}

// inverse matrix
template<typename T>
Mat<T> Mat<T>::inv() const 
{
    return (1.0 / det()) * adj();
}

//
// template functions
//

// scalar * matrix
template<typename T>
Mat<T> operator*(const T c, const Mat<T>& a)
{
    Mat<T> r(a.num_rows, a.num_cols);
    for (int row = 0; row < a.num_rows; ++row) {
        for (int col = 0; col < a.num_cols; ++col) {
            r.set(row, col, c * a.get(row, col) );
        }
    }
    return r;
}

//template <typename T>
//struct Model {
//    std::vector<Tri<T>> tris;
//    Mat<T> trans;
//    Mat<T> rot;
//    Mat<T> scale;
//    Model(std::initializer_list<Tri<T> tris);
//};

template <typename T>
Model<T>::Model()
{
    tris.clear();
    trans = Mat<T>(4, 4);
    rot = Mat<T>(4, 4);
    scale = Mat<T>(4, 4);
    trans.set_identity();
    rot.set_identity();
    scale.set_identity();
    color = {20, 20, 20, 20};
    draw_wire_frame = true;
    draw_face = true;
    draw_back_face = false;
}

//Model(std::initializer_list<Vec<T>> tris);
template <typename T>
Model<T>::Model(std::initializer_list<Vec<T>> verts, 
        Color color_arg, 
        bool draw_wire_frame_arg,
        bool draw_back_face_arg) 
{
    //template <typename T>
    //struct Model {
    //    std::vector<Tri<T>> tris;
    //    Mat<T> trans;
    //    Mat<T> rot;
    //    Mat<T> scale;
    //    Color color;
    //    Model(std::initializer_list<Vec<T>> tris);
    //    //Model(std::initializer_list<Tri<T>> tris);
    //};
    assert(verts.size() % 3 == 0);
    trans = Mat<T>(4, 4);
    rot = Mat<T>(4, 4);
    scale = Mat<T>(4, 4);
    trans.set_identity();
    rot.set_identity();
    scale.set_identity();
    color = color_arg;
    draw_wire_frame = draw_wire_frame_arg;
    draw_face = true;
    draw_back_face = draw_back_face_arg;
    auto it = verts.begin();
    int i = 0;
    while (it < verts.end()) {
        Vec<T> a; Vec<T> b; Vec<T> c;
        a = *it;
        b = *(it + 1);
        c = *(it + 2);
        //cout << "debugin' stuff: " << a << " " << b << " " << c << '\n';
        Tri<T> tri{ a, b, c, color, draw_wire_frame, true, draw_back_face}; 
        //Tri<T> tri{ *it, *(it + 1), *(it + 2), color, draw_wire_frame_arg}; 
        //tris[i] = tri;
        tris.push_back(tri);
        it += 3;
        i++;
    }
}

template <typename T>
bool Model<T>::load_from_file(string filename, 
        Color color_arg, 
        bool draw_wire_frame_arg,
        bool draw_back_face_arg)
{
    std::ifstream f(filename);

    if (!f.is_open()) {
        return false;
    }

    tris.clear();
    vector<Vec<T>> verts;
    color = color_arg;
    draw_wire_frame = draw_wire_frame_arg;
    draw_face = true;
    draw_back_face = draw_back_face_arg;
    
    while (!f.eof()) {
        char line[200];
        f.getline(line, 200);
        strstream s;
        s << line;

        char junk;
        if (line[0] == 'v') {
            Vec<T> v(4);
            s >> junk >> v.data[0] >> v.data[1] >> v.data[2];
            v.data[3] = 1.0;
            verts.push_back(v);
        } else if (line[0] == 'f') {
            int i[3];
            s >> junk >> i[0] >> i[1] >> i[2];
            i[0] -= 1; // file indexes from 1, instead of 0 (how this program indexes)
            i[1] -= 1; 
            i[2] -= 1;

            //  Tri(const Vec<T>& a, const Vec<T>& b, const Vec<T>& c, Color color,
            //          bool draw_wire_frame = true, bool draw_face = true);
            tris.push_back( Tri( verts[i[0]], verts[i[1]], verts[i[2]], color_arg, draw_wire_frame_arg, true, draw_back_face));
        } 

    }
}


//
// function templates
//

template<typename T>
Vec<T> operator*(const T c, const Vec<T> a)
{
    Vec<T> result(a.length);
    for (int i = 0; i < a.length; ++i) {
        result.data[i] = c * a.data[i];
    }
    return result;
}

template<typename T>
Vec<T> operator/(const T c, const Vec<T> rhs)
{
    Vec<T> result(rhs.length);
    for (int i = 0; i < rhs.length; ++i) {
        result.data[i] = c / rhs.data[i];
    }
    return result;

}

// row vector * matrix
template<typename T>
Vec<T> operator*(const Vec<T>& a, const Mat<T>& m)
{
    //assert(a.length == m.num_rows);
    //Vec<T> result(a.length);
    ////int result_index = 0;
    ////for (int i = 0; i < m.num_rows; ++i) {
    ////    result.data[result_index] = a 
    ////    result_index++;
    ////}
    ////T temp = 0.0;
    //for (int col = 0; col < m.num_cols; ++col) {
    //    //temp = 0.0;
    //    for (int i = col; i < m.num_cols * m.num_rows; i += m.num_cols) {
    //        // i is an index that can retrieve a matrix element
    //        // i will point to elements associated with the column vector
    //        result.data[col] += a.data[col] * m.data[i];
    //    }
    //}
    assert(a.length == m.num_rows);
    Vec<T> result(m.num_cols);
    for (int i = 0; i < result.length; ++i) {
        for (int j = 0; j < m.num_rows; ++j) {
            result.data[i] += a.data[j] * m.get(j, i); 
        }
    }
    return result;
    //assert(this->num_cols == b.num_rows);
    //Mat<T> r(num_rows, b.num_cols);
    //T sum = 0;
    //for (int row = 0; row < this->num_cols; ++row) {
    //    for (int col = 0; col < b.num_rows; ++col) {
    //        sum = 0;
    //        for (int k = 0; k < this->num_cols; ++k) {
    //            sum += this->get(row, k) * b.get(k, col);
    //        }
    //        r.set(row, col, sum);
    //    }
    //}
    //return r;
    //return result;
}

template <typename T> std::ostream& operator<<(std::ostream& os, const Vec<T>& obj) {
    os << "( ";
    for (int i = 0; i < obj.length; ++i) {
        os << obj.data[i] << " "; 
    }
    os << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Color& obj) {
    os << "( ";
    os << int(obj.r) << ", " << int(obj.g) << ", " << int(obj.b) << ", " << int(obj.a);
    os << ")";
    return os;
}

template <typename T> std::ostream& operator<<(std::ostream& os, const Mat<T>& m /*obj*/) {
    //os << "{ " << endl;
    //for (int i = 0; i < obj.length; ++i) {
    //    //os << "  ( " << obj.rows[i] << " )" << endl;
    //    os << "  " << obj.rows[i] << endl;
    //}
    //os << "}";
    //return os;
    
//    cout << "matrix is square of " 
//            << m.num_rows << " " 
//            << m.num_cols << endl;
    cout << "{ " << endl;
    for (int row = 0; row < m.num_rows; ++row) {
        for (int col = 0; col < m.num_cols; ++col) {
            cout << "  " <<  m.get(row, col) << " ";
        }
        cout << endl;
    }
    cout << "}";
    return os;
}

template <typename T>
void print_matrix(const Mat<T>& m) 
{
    
    cout << "matrix is square of " << m.num_rows << " " << m.num_cols << endl;
    cout << "{ " << endl;
    for (int row = 0; row < m.num_rows; ++row) {
        for (int col = 0; col < m.num_cols; ++col) {
            cout << m.get(row, col) << " ";
        }
        cout << endl;
    }
    cout << "}";
}

// return line that passes thru points pt0 and pt1.
// input: endpoints of segment. in (x, y) coordinates
// return: line expressed in homogenous coordinates ax + by + c = 0
template<typename T>
Vec<T> get_line_from_segment(const Vec<T>& pt0, const Vec<T>& pt1) {
    assert(pt0.length == pt1.length && pt0.length == 2);
    // make sure points are not identical
    assert(!(pt0.data[0] == pt1.data[0] && pt0.data[1] == pt1.data[1]));
    Vec<T> result(3);
    // find a. thus try to find slope, if slope infinity, a = 0
    T a = 0; T b = 0; T c = 0;
    if (pt1.data[0] - pt0.data[0] == 0.0) {
        // verticle line
        a = 1;
        b = 0;
        c = - pt0.data[0];
    } else {
        // not a verticle line
        b = 1;  // b * y
        T m = ( pt1.data[1] - pt0.data[1]) / (pt1.data[0] - pt0.data[0]);
        a = -m;
        // find c. thus find intercept. use pt0 and m
        c = m * pt0.data[0] - pt0.data[1];  // c = mx - y 
    }
    result.data[0] = a;
    result.data[1] = b;
    result.data[2] = c;
    return result;
}

// line to line intersection. Lines expressed as homogenous coordinates
// with form ax + by + c = 0. element 0 = a, element 1 = b, etc
// return: if return vector has 0 in 2nd element (c value), then no intersection
//         otherwise lines intersect at 0th and 1st element
template<typename T>
Vec<T> line_line_intersection(const Vec<T>& a, const Vec<T>& b) {
    assert(a.length == b.length && a.length == 3);
    Vec<T> result = a.cross(b);
    T w = result.data[2];
    if (w != 0.0) {
        return (1.0 / w) * result;    
    } else {
        return result; // no intersection
    }
}

template<typename T>
Vec<T> line_line_segment_intersection(const Vec<T>& line, const Vec<T>& pt0, const Vec<T>& pt1){
    assert(pt0.length == 2 && pt1.length == 2 && line.length == 3);
    // make sure points are not identical. so they can describe a line
    assert(!(pt0.data[0] == pt1.data[0] && pt0.data[1] == pt1.data[1]));
    Vec<T> segment_line(3);
    segment_line = get_line_from_segment(pt0, pt1);
    Vec<T> intersection_coords(3);  // use copy constructor instead?
    intersection_coords = line_line_intersection(line, segment_line); 
    if (intersection_coords.data[2] == 0.0) {
        // no intersection! Or lines could be the same?
        // put in large negative numbers as a warning
        // 0.0 at index 2 shows no intersection
        return Vec<T>{-1000000.0, -1000000.0, 0.0};
    }
    // intersection candidate found, now check if it exists on segment

    
    // check which parametric equation to use
    T t = 0.0; // parameter variable
    if (pt1.data[0] != pt0.data[0]) {  
        // can use parametric eq with x and t
        t = (intersection_coords.data[0] - pt0.data[0]) / (pt1.data[0] - pt0.data[0]);
    } else {
        // can use parametric eq with y and t
        t = (intersection_coords.data[1] - pt0.data[1]) / (pt1.data[1] - pt0.data[1]);
        //cout << "  parametric eq with y and t" << endl;
        //cout << "  intersection " << intersection_coords << endl;
        //cout << "  pt0 " << pt0 << endl;
        //cout << "  pt1 " << pt1 << endl;
        //cout << "  t = " << t << endl;
    }
    //if ( t >= 0.0 && t <= 1.0) {
    if ( t >= 0.0 && t <= 1.0) {
        // intersection_coords is a valid intersection between a line and a segment
    } else {
       // not a valid candidate
       intersection_coords.data[2] = 0.0;
    }
    // intersection if data[2] != 0.0 
    return intersection_coords;  
}

// Check if point is in a region defined by a line and another point
// input: line (ax + by + c = 0) defines 2 regions, element count 3
//        pt_region, selects a region, element count 2
//        pt_test, tested against pt_region, to see if both points 
//                 are in same region
template<typename T>
bool point_in_region(const Vec<T>& line, const Vec<T>& pt_region, const Vec<T>& pt_test) 
{
    // check dimensions
    assert(line.length == 3); //failed
    assert(pt_region.length == 2); 
    assert(pt_test.length == 2);
    // make sure pt_region isn't on line
    assert(line.data[0] * pt_region.data[0] + line.data[1] * pt_region.data[1] 
            + line.data[2] != 0.0); 
    // make sure pt_test isn't on line
    assert(line.data[0] * pt_test.data[0] + line.data[1] * pt_test.data[1] 
            + line.data[2] != 0.0); 
   
    T region_sign = line.data[0] * pt_region.data[0] + line.data[1] 
            * pt_region.data[1] + line.data[2]; 

    T test_sign = line.data[0] * pt_test.data[0] + line.data[1] 
            * pt_test.data[1] + line.data[2]; 

    // if both have the same sign, they're in the same region
    return (region_sign * test_sign > 0.0);
}

// clips polygon inside region
// 2 dimension clipping
template<typename T>
void clip_polygon(const std::vector<Vec<T>>& arg_input_vertices, 
        const std::vector<Vec<T>>& edges,
        const Vec<T>& pt_inside, 
        std::vector<Vec<T>>& output_vertices)
{
    using std::vector;
    // check input 
    assert(arg_input_vertices.size() > 0 && arg_input_vertices[0].length == 2); // input_vertices length = 2 
    assert(edges.size() > 0 && edges[0].length == 3); // lines expressed as homogenous coords so length = 3
    assert(pt_inside.length == 2);

    vector<Vec<T>> input_vertices(2 * arg_input_vertices.size() + 1); // minimize resizes
    output_vertices = arg_input_vertices;
    
    for(auto& edge : edges) {
        input_vertices = output_vertices;
        output_vertices.clear();
        for (int i = 0; i < input_vertices.size(); ++i) {
            Vec<T>& current = input_vertices[i];
            int prev_index = (i - 1 >= 0) 
                    ? i - 1 : input_vertices.size() - 1;
            
            Vec<T>& prev = input_vertices[prev_index];
            if (point_in_region(edge, pt_inside, current)) {
                if (!point_in_region(edge, pt_inside, prev)) {
                    // find intersection 
                    Vec intersection = line_line_segment_intersection(edge, prev, current);
                    assert(intersection.data[2] != 0.0); // there has to be a crossing
                    output_vertices.push_back(Vec{intersection.data[0], intersection.data[1]});
                }
                output_vertices.push_back(current);
            } else if (point_in_region(edge, pt_inside, prev)) {
                Vec intersection = line_line_segment_intersection(edge, prev, current);
                assert(intersection.data[2] != 0.0); // there has to be a crossing
                output_vertices.push_back(Vec{intersection.data[0], intersection.data[1]});
            }
        }
    }
}

// return line that passes thru points p0 and p1.
// input: endpoints of segment. in (x, y, z) coordinates
// return: line expressed as dir_vec * t + point, 6 numbers total  
template<typename T>
Vec<T> get_line_from_segment3d(const Vec<T>& p0, const Vec<T>& p1) {
    assert(p0.length == p1.length && p0.length == 3);
    // make sure points are not identical
    assert(p0 != p1);
    Vec<T> dir = p1 - p0;
    return Vec{dir.data[0], dir.data[1], dir.data[2], 
            p0.data[0], p0.data[1], p0.data[2]};
}

// gets distance point is from plane
// plane has length 4, point length 3
template <typename T>
T dist_from_plane(const Vec<T>& plane, const Vec<T>& point)
{
    assert(plane.length == 4);
    assert(point.length == 3);

    // if 0, point on plane. if > 0, point above plane in direction of normal
    return plane.data[0] * point.data[0] + plane.data[1] * point.data[1] 
            + plane.data[2] * point.data[2] + plane.data[3];
}

template <typename T>
bool point_in_region3d(const Vec<T>& plane, const Vec<T>& _p_region, const Vec<T>& _p_test)
{
    auto p_region = _p_region;
    auto p_test = _p_test;
    assert(plane.length == 4);
    if (p_region.length == 4) {
        p_region = p_region.slice(0, 3);
    }
    assert(p_region.length == 3);
    if (p_test.length == 4) {
        p_test = p_test.slice(0, 3);
    }
    assert(p_test.length == 3);
    //cout << "p_test.length == " << p_test.length << endl;
    //cout << "_p_test.length == " << _p_test.length << endl;
    
    T p_region_dist = dist_from_plane(plane, p_region);
    T p_test_dist = dist_from_plane(plane, p_test);

    // if same sign (+ + or - -) then two points are in the same region
    return p_region_dist * p_test_dist > 0.0;
    //// old below
    //assert(plane.length == 4);
    //assert(p_region.length == 3);
    //cout << "p_test.length == " << p_test.length << endl;
    //assert(p_test.length == 3);
    //
    //T p_region_dist = dist_from_plane(plane, p_region);
    //T p_test_dist = dist_from_plane(plane, p_test);

    //// if same sign (+ + or - -) then two points are in the same region
    //return p_region_dist * p_test_dist > 0.0;
}

// returns point where line and plane intersect,
// params: line dimension 6, (ax, by, cz) + (d, e, f)
//         plane dimension 4, ax + by + cz + w = 0
// return: returns point, dimension 4, if last element 0, no intersection
template <typename T>
Vec<T> line_plane_intersect(const Vec<T>& line, const Vec<T>& plane)
{
    assert(line.length == 6); assert(plane.length == 4); 
    // find parameter t 
    T denominator = plane.data[0] * line.data[0] + plane.data[1] * line.data[1] 
            + plane.data[2] * line.data[2];
    if (denominator == 0.0) {
        // no intersection, last number 0.0
        // large numbers give hint that not intersection found
        return Vec{-10000.0, -10000.0, -10000.0, 0.0}; //
    }
    T t = (-plane.data[0] *line.data[3] - plane.data[1] * line.data[4] 
            - plane.data[2] * line.data[5] - plane.data[3]) / denominator;
    //cout << "   >> line_plane_intersect " << endl;
    //cout << "      line: " << line << endl;
    //cout << "      plane: " << plane << endl;
    //cout << "      t: << " << t << endl;
    // last element set to 1.0 to show this is a valid coordinate
    Vec result = {t * line.data[0] + line.data[3], 
                  t * line.data[1] + line.data[4], 
                  t * line.data[2] + line.data[5], 
                  1.0 };
    auto temp = plane.data[0] * result.data[0] + plane.data[1] * result.data[1] + plane.data[2] * result.data[2] 
            + plane.data[3]; 
    if (temp <= 1.0e-9) { // TODO should negative temps be considered?
        temp = 0.0;
    }
    //cout << "      check answer (should be 0) " << temp  << endl;
    if (temp != 0.0) {
        cout << "temp " << temp << endl;
        assert (temp == 0.0);
    }
    return Vec{t * line.data[0] + line.data[3], t * line.data[1] + line.data[4], t * line.data[2] + line.data[5], 1.0};
}


template <typename T>
bool point_on_line3d(const Vec<T>& ptest, const Vec<T>& line)
{
    assert(ptest.length == 3);    
    assert(line.length == 6);    
    Vec<T> line_point = line.slice(3, 3);
    Vec<T> line_dir = line.slice(0, 3);
    if (ptest == line_point) { return true; }
    Vec<T> cross_prod = line_dir.cross(ptest - line_point);
    // if cross product magnitude == 0. then point is on the line
    return cross_prod.mag() == 0.0;
}

// p0 and p1 comprise the line segment.
template <typename T>
bool point_on_line_segment3d(const Vec<T>& ptest, const Vec<T>& p0, const Vec<T>& p1)
{
    assert(p0.length == 3 && p1.length == 3 & ptest.length == 3);    
    Vec line = get_line_from_segment3d(p0, p1);
    num_calls++;
    if (!point_on_line3d(ptest, line)) {
        //cout << "  point isn't on line (no need to check segment)\n";
        return false;
    }
    // still need to check if on segment, rather than full line
    Vec dir = line.slice(0, 3);
    Vec a_point = line.slice(3, 3);
    
    // solve for t, t should theoretically be the same for all 3 equations
    using std::array;
    array<T, 3> t;
    //cout << ">>point_on_line_segment3d"  << endl;
    //cout << "   ptest = " << ptest << endl;
    //cout << "   p1 = " << p1 << endl;
    //cout << "   p0 = " << p0 << endl;
    //cout << "   dir = " << dir << endl;
    //cout << "   a_point = " << a_point << endl;
    //cout << "   num_calls = " << num_calls << endl;
    //cout << "   printing out all t values for point_on_line_segment3d" << endl;
    for (int i = 0; i < 3; ++i) {
        if (dir.data[i] == 0.0) {
            t[i] = 0.0;//-0.0;
        } else {
            t[i] = (ptest.data[i] - a_point.data[i]) / dir.data[i];
        }
        //cout << "   " <<  t[i] << endl;
    }

    //Vec t = (1.0/dir) * (ptest - a_point);

    if ( t[0] >= 0.0 && t[0] <= 1.0 
            && t[1] >= 0.0 && t[1] <= 1.0 
            && t[2] >= 0.0 && t[2] <= 1.0 ) 
    {
        //cout << "   point IS on line segment\n";
        return true;
    } else {
        //cout << "   point IS NOT line segment\n";
        return false;
    }

}

// returns point where line segment and plane intersect,
// params: line dimension 6, (ax, by, cz) + (d, e, f)
//         plane dimension 4, ax + by + cz + w = 0

template <typename T>
Vec<T> plane_line_segment_intersect(const Vec<T>& plane, const Vec<T>& p0, const Vec<T>& p1)
{
    //cout << "plsi p0 and p1 length " << p0.length << ", " << p1.length << endl; 
    assert(p0.length == 3);
    assert(p1.length == 3);
    assert(plane.length == 4); 
    Vec line = get_line_from_segment3d(p0, p1);
    Vec intersection = line_plane_intersect(line, plane);

    if (intersection.data[3] == 0.0) {
        // line doesn't even intersect plane, let alone segment
        // 0.0 in last element signify no intersection
        return Vec{-10000.0, -10000.0, -10000.0, 0.0}; 
    }
    // does the intersection point lie on the line segment?
    // TODO problem with point_on_line_segment3d?
    if (point_on_line_segment3d(intersection.slice(0,3), p0, p1)) {
        return intersection;
    } else {
        // point not on segment, (although it is on line)
        intersection.data[3] = 0.0;
        return intersection;
    }
}

// clips polygon inside region
// 3 dimensional clipping
template<typename T>
void clip_polygon3d(
        const std::vector<Vec<T>>& arg_input_vertices, 
        const std::vector<Vec<T>>& planes,
        const Vec<T>& pt_inside, 
        std::vector<Vec<T>>& output_vertices)
{
    using std::vector;
    // check input 
    //3 and 4
    //cout << " hi there! " << arg_input_vertices.size() << " and " << arg_input_vertices[0].length << endl;
    assert(arg_input_vertices.size() > 0 && arg_input_vertices[0].length == 3); // input_vertices length = 3 
    assert(planes.size() > 0 && planes[0].length == 4); // plane edges have 4 dimensions
    assert(pt_inside.length == 3);
    vector<Vec<T>> input_vertices(2 * arg_input_vertices.size() + 1); // minimize resizes
    
    output_vertices = arg_input_vertices;
    
    for(auto& edge : planes) {
        input_vertices = output_vertices;
        output_vertices.clear();
        for (int i = 0; i < input_vertices.size(); ++i) {
            Vec<T>& current = input_vertices[i];
            int prev_index = (i - 1 >= 0) 
                    ? i - 1 : input_vertices.size() - 1;
            //cout << "current_index " << i << endl;
            //cout << "prev_index " << prev_index << endl;
            //cout << "current vertex " << input_vertices[i] << endl;
            //cout << "prev vertex " << input_vertices[prev_index] << endl;
            //cout << "edge " << edge << endl;
            //cout << endl;
            Vec<T>& prev = input_vertices[prev_index];
            if (point_in_region3d(edge, pt_inside, current)) {
                if (!point_in_region3d(edge, pt_inside, prev)) {
                    // find intersection 
                    Vec intersection = plane_line_segment_intersect(edge, prev, current);
                    //assert(intersection.data[3] != 0.0); // there has to be a crossing //TODO  failed?
                    if (intersection.data[3] == 0.0) {
                        // TODO problem in plane_line_segment_intersect?
                        //cout << ">>clip_polygon3d\n";
                        //cout << " assertion failed, current inside, prev outside? " << endl;
                        //cout << " i: " << i << endl;
                        //cout << " num_calls: " << num_calls << endl;
                        //cout << " edge: " << edge << endl;
                        //cout << " current: " << current << endl<< " prev: " << prev << endl;
                        //cout << " intersection: " << intersection << endl;
                        //cout << " set breakpoint here" << endl;
                        //TODO check //assert(intersection.data[3] != 0.0); // there has to be a crossing //TODO  failed?
                    }
                    output_vertices.push_back(intersection.slice(0, 3));
                    //output_vertices.push_back(Vec{intersection.data[0], intersection.data[1]});
                }
                output_vertices.push_back(current);
            } else if (point_in_region3d(edge, pt_inside, prev)) {
                Vec intersection = plane_line_segment_intersect(edge, prev, current);
                //assert(intersection.data[3] != 0.0); // there has to be a crossing
                if (intersection.data[3] == 0.0) {
                    //cout << ">>clip_polygon3d\n";
                    //cout << "assertion failed, current outside, prev inside? " << endl;
                    //cout << " num_calls: " << num_calls << endl;
                    //cout << " edge: " << edge << endl;
                    //cout << " current: " << current << endl << " prev: " << prev << endl;
                    //cout << " intersection: " << intersection << endl;
                    // TODO checkassert(intersection.data[3] != 0.0); // there has to be a crossing //TODO  failed?
                }
                output_vertices.push_back(intersection.slice(0, 3));
                //output_vertices.push_back(Vec{intersection.data[0], intersection.data[1]});
            }
        }
    }
}

// 3 dimensional clipping
template<typename T>
void clip_line3d(
        const vector<Vec<T>>& _arg_input_vertices, 
        const vector<Vec<T>>& planes,
        const Vec<T>& pt_inside, 
        vector<Vec<T>>& _output_vertices)
{
    //cout << "start clip_line3d\n";

    // kludge, convert homogeneous 4d coordinates to 3d
    vector<Vec<T>> arg_input_vertices{};
    vector<Vec<T>> output_vertices{};
    for (const auto& v : _arg_input_vertices) {
        arg_input_vertices.push_back(v.slice(0, 3));
    }
    // checking
    //for (const auto& v : arg_input_vertices) { 
    //    cout << "input vertices length " << v.length << endl;
    //}

    // check input 
    assert(arg_input_vertices.size() > 0 
            && arg_input_vertices[0].length == 3); // input_vertices length = 3 
    assert(planes.size() > 0 && planes[0].length == 4); // plane edges have 4 dimensions
    assert(pt_inside.length == 3);
    vector<Vec<T>> input_vertices(2 * arg_input_vertices.size() + 1); // minimize resizes
    // generally 2 Vec<T>'s per line, although a line could have more... always multiple
    // of 2 though
    
    output_vertices = arg_input_vertices;
    
    //// checking
    //for (const auto& v : output_vertices) { 
    //    cout << "output_vertices length " << v.length << endl;
    //}
    int count_loops = 0;
    for(auto& edge : planes) {
        //cout << "count_loops " << count_loops << endl; 
        count_loops++;
        input_vertices = output_vertices; 
        output_vertices.clear();
        for (int i = 0; i < input_vertices.size(); i = i + 2) {
            // TODO fix later for lines buffers with more than two points 
            Vec<T>& current = input_vertices[i + 1]; 
            Vec<T>& prev = input_vertices[i]; 
            //cout << "current length " << current.length << endl;
            //cout << "prev length " << prev.length << endl;

            if (point_in_region3d(edge, pt_inside, current) 
                    && point_in_region3d(edge, pt_inside, prev)) {
                output_vertices.push_back(prev);
                output_vertices.push_back(current);
            } else if (point_in_region3d(edge, pt_inside, current) 
                    && !point_in_region3d(edge, pt_inside, prev)) {
                Vec p = plane_line_segment_intersect(edge, prev, current);
                output_vertices.push_back(p.slice(0, 3));
                output_vertices.push_back(current);
            } else if (!point_in_region3d(edge, pt_inside, current) 
                    && point_in_region3d(edge, pt_inside, prev)) {
                output_vertices.push_back(prev);
                Vec p = plane_line_segment_intersect(edge, prev, current);
                output_vertices.push_back(p.slice(0, 3));
            }
    

        }
        
    }
    // more kludgy stuff. convert 3d coordinates to 4d homogeneous coordinates
    //auto temp_outputs = output_vertices;
    //output_vertices.clear();
    for (auto& temp : output_vertices) {
        _output_vertices.emplace_back( Vec<double>{
            temp.data[0], 
            temp.data[1], 
            temp.data[2], 
            1.0} ); 
    }
}

// p0 and p1 and 2d points. w is the 2nd element
// boundary is a line, ax + by + c = 0. elements are a, b, c
void clip_line_new(
        const Vec<double>& p0, 
        const Vec<double>& p1, 
        const Vec<double>& boundary, 
        vector<Vec<double>>& output_line_segment)
{
    assert(p0.length == 2 && p1.length == 2);
    assert(boundary.length == 3);
    
    // determine boundary

    // check if both points are outside boundary ( in a seperate function?)
    //if (p0.data[0] == 

}
    
//template <typename T>
//PlaneLineSegmentIntersectResult<T> plane_line_segment_intersect(const Vec<T>& plane, const Vec<T>& p0, const Vec<T>& p1)
//{
//    //cout << "plsi p0 and p1 length " << p0.length << ", " << p1.length << endl; 
//    assert(p0.length == 3);
//    assert(p1.length == 3);
//    assert(plane.length == 4); 
//    Vec line = get_line_from_segment3d(p0, p1);
//    Vec intersection = line_plane_intersect(line, plane);
//    PlaneLineSegmentIntersectResult<T> result;
//
//    if (intersection.data[3] == 0.0) {
//        // line doesn't even intersect plane, let alone segment
//        // 0.0 in last element signify no intersection
//        return Vec{-10000.0, -10000.0, -10000.0, 0.0}; 
//    }
//    // does the intersection point lie on the line segment?
//    // TODO problem with point_on_line_segment3d?
//    if (point_on_line_segment3d(intersection.slice(0,3), p0, p1)) {
//        return intersection;
//    } else {
//        // point not on segment, (although it is on line)
//        intersection.data[3] = 0.0;
//        return intersection;
//    }
//}

//// render triangle between 3 points
//template<typename T>
//void draw_triangle(const Vec<T>& a, const Vec<T>& b, const Vec<T>& c)
//{
//    // using SDL_RenderDrawLine calls to draw a triangle
//    // call line function from everypoint on the line b to c, back to a
//    // use integers to count pixels
//    
//    //find manhattan distance from b to c
//    Vec bc = b - c;
//    Vec bc_unit = bc.unit(); 
//    Vec cur = c;
//    for (size_t i = 0; i < bc.data[0] + bc.data[1] < ++i) {
//        SDL_RenderDrawLine(gsdl.renderer, a.data[0],  a.data[1],
//                cur.data[0], cur.data[1]);
//        cur += bc.unit();
//        //kludge
//        if (cur.magnitude() > bc.magnitude())
//        { break;}
//    }
//}

// rotate about y (rotate heading)
Mat<double> rotate_y(double rad) 
{
    // note! returns 4x4 matrix
    // doesn't rotate by 15 degrees, but by 7.5 degrees
    Mat<double> m(4, 4);
    //rad = (PI * 1.5)/360.0;
    m.set(0, 0, cos(rad));
    m.set(0, 2, -sin(rad));
    m.set(1, 1, 1.0);
    m.set(2, 0, sin(rad));
    m.set(2, 2, cos(rad));
    m.set(3, 3, 1.0);
    return m;
}

// rotate about x (rotate pitch)
Mat<double> rotate_x(double rad)
{
    // note! returns 4x4 matrix
    Mat<double> m(4, 4);
    //rad = (PI * 1.5)/360.0;
    m.set(0, 0, 1.0);
    m.set(1, 1, cos(rad));
    m.set(2, 1, -sin(rad));
    m.set(1, 2, sin(rad));
    m.set(2, 2, cos(rad));
    m.set(3, 3, 1.0);
    return m;
}

#endif
