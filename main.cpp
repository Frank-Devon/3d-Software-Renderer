#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <exception>
#include "vec.hpp"
#include "test.hpp"
#include "sdl_init.hpp"

//ctags -R --languages=C++ --extra=+fq --exclude=.git .
using std::cout, std::endl, std::vector, std::array;

// vertexes of object
#define NUM_VERTICES 9
#define NUM_MESHES 3
#define ROTATE_FIRST 0
#define RENDER_LINES 1
#define NUM_LINES 22

array<vector<Vec<double>>, NUM_MESHES> models  =
{
    // object 0
    vector{
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0}
    }
    ,
        // object 1
        vector{
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0,  1.0, -1.0, 1.0},
        Vec{ 1.0,  1.0,  1.0, 1.0},
        Vec{-1.0,  1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{-1.0, -1.0, -1.0, 1.0},
        Vec{ 1.0, -1.0,  1.0, 1.0},
        Vec{-1.0, -1.0,  1.0, 1.0}

        }
        ,
        // coordinate marker
        vector{
            // x axis
            Vec{0.0, 0.0, 0.0, 1.0 },
            Vec{1.0, 0.0, 0.0, 1.0 },
            Vec{1.0, 0.0, 0.0, 1.0 },
            
            // y axis
            Vec{0.0, 0.0, 0.0, 1.0 },
            Vec{0.0, 1.0, 0.0, 1.0 },
            Vec{0.0, 1.0, 0.0, 1.0 },
            
            // z axis
            Vec{0.0, 0.0, 0.0, 1.0 },
            Vec{0.0, 0.0, 1.0, 1.0 },
            Vec{0.0, 0.0, 1.0, 1.0 }
        }
        //vector{
        //    Vec{0.0, 0.0, 0.0, 1.0 },
        //    Vec{0.0, 1.0, 0.0, 1.0 },
        //    Vec{0.0, 0.0, 1.0, 1.0 },
        //    
        //    Vec{0.0, 0.0, 0.0, 1.0 },
        //    Vec{0.0, 1.0, 0.0, 1.0 },
        //    Vec{1.0, 0.0, 0.0, 1.0 }
        //}
};
// vertices drawn to screen
array<vector<Vec<double>>, NUM_MESHES> model_outputs = 
{
    vector<Vec<double>>{},
    vector<Vec<double>>{},
    vector<Vec<double>>{}
};

array<vector<Vec<double>>, NUM_MESHES> model_outputs_clipped;

array<vector<Vec<double>>, NUM_LINES> lines_world;
array<vector<Vec<double>>, NUM_LINES> lines_outputs; 
array<vector<Vec<double>>, NUM_LINES> lines_outputs_clipped;
int lines_outputs_clipped_count = 0;

int win_res_x = 1000;
int win_res_y = 1000;
// meshes
//Control_Select control_select = Control_Select::model_0;
enum Control_Select : unsigned int { m_0 = 0, m_1 };
Control_Select control_select = m_0;
bool print_clip_vertices = false; // TODO delete?
bool print_info = false;
std::array<std::array<Vec<double>, NUM_VERTICES>, NUM_MESHES> meshes; 

//// model to world matrices
std::array<Mat<double>, NUM_MESHES> model_transforms =  { 
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 4.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 7.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
};

// initial positions of objects here
std::array<Mat<double>, NUM_MESHES> model_translate =  { 
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{-2.0, 0.0, 6.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{2.0, 0.0, 9.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
};
// initial rotation = none
std::array<Mat<double>, NUM_MESHES> model_rotate =  { 
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
};

// initial scale = 1.0 aka none
std::array<Mat<double>, NUM_MESHES> model_scale =  { 
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
    ,
    Mat{  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    }
    
};

Mat<double> world_to_camera = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};

Mat<double> world_to_camera_translate = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, -3.0, 11.5, 1.0}
};

Mat<double> world_to_camera_rotate = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, -1.0, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};

// translation matrices
Mat<double> translate_x_plus = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.1, 0.0, 0.0, 1.0}
};
Mat<double> translate_x_minus = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{-0.1, 0.0, 0.0, 1.0}
};
Mat<double> translate_y_plus = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.1, 0.0, 1.0}
};
Mat<double> translate_y_minus = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, -0.1, 0.0, 1.0}
};
Mat<double> translate_z_plus = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.1, 1.0}
};
Mat<double> translate_z_minus = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, -0.1, 1.0}
};
// end translation matrices
Mat<double> scale_enlarge = {  
    Vec{1.05, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.05, 0.0, 0.0},
    Vec{0.0, 0.0, 1.05, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};
Mat<double> scale_shrink = {  
    Vec{0.95, 0.0, 0.0, 0.0}, 
    Vec{0.0, 0.95, 0.0, 0.0},
    Vec{0.0, 0.0, 0.95, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};


Mat<double> ortho_matrix(4, 4);
Mat<double> clip_matrix(4, 4); // not used yet, to be used for cliping vertices outside view frustrum.
Mat<double> turn_left(4, 4);          
Mat<double> turn_right(4, 4);          

// clipped vertices go here
vector<Vec<double>> output_vertices;
vector<int> output_vertices_per_polygon;
vector<Vec<double>> temp_vertices;
vector<Vec<double>> clip_planes = {
    Vec{1.0, 0.0, 0.0, -90.0},
    Vec{1.0, 0.0, 0.0, 90.0},
    Vec{0.0, 1.0, 0.0, -90.0},
    Vec{0.0, 1.0, 0.0, 90.0},
    Vec{0.0, 0.0, 1.0, -90.0},
    Vec{0.0, 0.0, 1.0, 90.0}
};
//vector<Vec<double>> clip_planes = {
//    Vec{1.0, 0.0, 0.0, -1.0},
//    Vec{1.0, 0.0, 0.0, 1.0},
//    Vec{0.0, 1.0, 0.0, -1.0},
//    Vec{0.0, 1.0, 0.0, 1.0},
//    Vec{0.0, 0.0, 1.0, -1.0},
//    Vec{0.0, 0.0, 1.0, 1.0}
//};

// render triangle between 3 points
template<typename T>
void draw_triangle(const Vec<T>& a, const Vec<T>& b, const Vec<T>& c)
{
    // using SDL_RenderDrawLine calls to draw a triangle
    // call line function from everypoint on the line b to c, back to a
    // use integers to count pixels
    
    //find manhattan distance from b to c
    Vec bc = b - c;
    Vec bc_unit = bc.unit(); 
    Vec cur = c;
    for (size_t i = 0; i < abs(bc.data[0]) + abs(bc.data[1]);  ++i) {
        SDL_RenderDrawLine(gsdl.renderer, a.data[0],  a.data[1],
                cur.data[0], cur.data[1]);
        cur = cur +  bc.unit();
        //kludge
        if ((cur - c).mag() > bc.mag())
        { break;}
    }
}

void test_matrix()
{
    cout << "test_matrix start\n" ;
    double turn1 = 0.4637; // 26.57 degrees CCW
    double turn2 = 1.5708; // 90 degrees CCW
    Mat<double> T1 = {  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{sqrt(5), 0.0, 0.0, 1.0}
    };
    Mat<double> T2 = {  
        Vec{1.0, 0.0, 0.0, 0.0}, 
        Vec{0.0, 1.0, 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{2.0 * sqrt(5), 0.0, 0.0, 1.0}
    };
    Mat<double> R1 = {  // rotate 26.57deg CCW
        Vec{cos(turn1), sin(turn1), 0.0, 0.0}, 
        Vec{-sin(turn1), cos(turn1), 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    };
    Mat<double> R2 = { // rotate 90.00deg CCW 
        Vec{cos(turn2), sin(turn2), 0.0, 0.0}, 
        Vec{-sin(turn2), cos(turn2), 0.0, 0.0},
        Vec{0.0, 0.0, 1.0, 0.0},
        Vec{0.0, 0.0, 0.0, 1.0}
    };
    Vec<double> v{1.0, 1.0, 0.0, 1.0}; 
    auto vT1 = v * T1;
    auto vR1 = v * R1;
    auto T1R1 = T1 * R1; 
    auto R1T1 = R1 * T1;
    auto vT1R1 = v * T1 * R1;
    auto vR1T1 = v * R1 * T1;
    auto vT2R2 = v * T2 * R2;
    auto vT1R1T2R2 = v * T1 * R1 * T2 * R2;

    cout  << "vT1  \n" << vT1   << endl;
    cout  << "vR1  \n" << vR1   << endl;
    cout  << "T1R1 \n" << T1R1  << endl;
    cout  << "R1T1 \n" << R1T1  << endl;
    cout  << "vT1R1\n" << vT1R1 << endl;
    cout  << "vR1T1\n" << vR1T1 << endl;
    cout  << "vT2R2\n" << vT2R2 << endl;
    cout << "vT1R1T2R2\n" << vT1R1T2R2 << endl;
    cout << "test_matrix end\n" ;
}

//Mat<double> ConstructCamera(const Vec<double>& pos, const Vec<double>& look_at,
//        const Vec<double>& up)
//{
//    
//}

void construct_grid() 
{
    // lines_world
    assert(lines_world.size() == 22);
    //for(vector<Vec<double>>& points : lines_world) {
    // construct lines with constant z values (appears horizontal)
    int i = 0; // index of lines_world
    for (int j = -5; j < 6; ++j) {
        double x_start = -5.0;
        double x_end   = 5.0;
        lines_world[i].push_back(Vec<double>{ x_start, 0.0, (double)j, 1.0});
        lines_world[i].push_back(Vec<double>{ x_end, 0.0, (double)j, 1.0});
        i++;
    }
    // construct lines with constant x values
    for (int j = -5; j < 6; ++j) {
        double z_start = -5.0;
        double z_end   = 5.0;
        lines_world[i].push_back(Vec<double>{ (double)j, 0.0, z_start, 1.0});
        lines_world[i].push_back(Vec<double>{ (double)j, 0.0, z_end, 1.0});
        i++;
    }
}

int main() 
{
    test_matrix();
    int index = 0;
    // initialize model_outputs?
    for(auto& model : models) {
        model_outputs[index] = model; 
        index++;
    }
    construct_grid();
    //cout << "testing grid\n";
    //for(auto x : lines_world) {
    //    // display lines coords
    //    cout << x[0] << ", " << x[1] << endl;
    //}
    // initialize lines_outputs
    for(auto& line : lines_outputs){
        line.push_back(Vec<double>{});
        line.push_back(Vec<double>{});
    }
    // initialize lines_outputs_clipped
    for(auto& line : lines_outputs_clipped){
        line.push_back(Vec<double>{});
        line.push_back(Vec<double>{});
    }
    //index = 0;
    //for(auto& line : lines_world) {
    //    lines_outputs_clipped[index] = line; 
    //    index++;
    //}
    //getchar();


    // doesn't rotate by 15 degrees, but by 7.5 degrees
    double rad = (PI * 1.5)/360.0;
    turn_left.set(0, 0, cos(rad));
    turn_left.set(0, 2, -sin(rad));
    turn_left.set(1, 1, 1.0);
    turn_left.set(2, 0, sin(rad));
    turn_left.set(2, 2, cos(rad));
    turn_left.set(3, 3, 1.0);
    rad = - rad;
    turn_right.set(0, 0, cos(rad));
    turn_right.set(0, 2, -sin(rad));
    turn_right.set(1, 1, 1.0);
    turn_right.set(2, 0, sin(rad));
    turn_right.set(2, 2, cos(rad));
    turn_right.set(3, 3, 1.0);

    auto test0 = turn_left * turn_left * turn_left;
    cout << "test shiet" << endl << test0 << endl;
    // camera space to clip space
    double near = 1.00;  // 1.0
    double far = 20000000.0;//200000.0;  // 20.0
    double projection_width = 0.6;
    double projection_height = 0.6;
    clip_matrix.set(0, 0, near );
    clip_matrix.set(1, 1, near );
    clip_matrix.set(2, 2, (near + far));
    clip_matrix.set(2, 3, 1.0 );
    clip_matrix.set(3, 2, -1.0*( near * far) ); 
    //clip_matrix.set(0, 0, near / (projection_width * 2));
    //clip_matrix.set(1, 1, near / (projection_height * 2));
    //clip_matrix.set(2, 2, (near + far)/(far - near));
    //clip_matrix.set(2, 3, 1.0 );
    //clip_matrix.set(3, 2, (- 2 * near * far) / (far - near));
    cout << "here's the clip matrix" << endl;
    cout << clip_matrix << endl;
    // ortho matrix set up
    double r = 1.0, l = -1.0, t = 1.0, b = -1.0;
    ortho_matrix.set(0, 0, 2.0/(r-l));
    ortho_matrix.set(3, 0, - (r+l)/(r-l));
    ortho_matrix.set(1, 1, 2.0/(t-b));
    ortho_matrix.set(3, 1, - (t+b)/(t-b));
    ortho_matrix.set(2, 2, 2.0/(near-far));
    ortho_matrix.set(3, 2, - (near+far)/(near-far));
    ortho_matrix.set(3, 3, 1.0);



    using std::cout, std::endl;
    //test1();


    // TODO
    // Important assumption: input vertices count == output vertices count.
    // this only true when there's no clipping. fix later 
    //cam.reserve

    //getchar();
    bool exit = false;    
    cout << "gsdl status: " << gsdl.created_ok << endl;
    SDL_Event e;
    Uint64 time_frame_start;
    Uint64 time_frame_duration = 0;

    bool posted_info = true;
    while (!exit) {
        time_frame_start = SDL_GetTicks64();
        while (SDL_PollEvent(&e) != 0) {
            switch (e.type) {
                case SDL_QUIT:
                    exit = true;
                    break;
                case SDL_KEYDOWN:
                    // rotate selected object
                    if (e.key.keysym.sym == SDLK_q) {
                        model_rotate[control_select] = model_rotate[control_select] * turn_left;
                    }
                    if (e.key.keysym.sym == SDLK_e) {
                        model_rotate[control_select] = model_rotate[control_select] * turn_right;
                    }

                    if (e.key.keysym.sym == SDLK_1) {
                        cout << "1 pressed" << endl;
                        control_select = m_0;
                    }
                    if (e.key.keysym.sym == SDLK_2) {
                        cout << "2 pressed" << endl;
                        control_select = m_1;
                    }
                    // translate selected model
                    if (e.key.keysym.sym == SDLK_a) { // left
                        model_translate[control_select].add_to_row(3, Vec{-0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_d) { // right
                        model_translate[control_select].add_to_row(3, Vec{0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_w) { // forward
                        model_translate[control_select].add_to_row(3, Vec{0.00,0.0,0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_s) { // backward
                        model_translate[control_select].add_to_row(3, Vec{0.00,0.0,-0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_z) {  // below
                        model_translate[control_select].add_to_row(3, Vec{0.00,-0.08,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_x) { // above
                        model_translate[control_select].add_to_row(3, Vec{0.00,0.08,0.0,0.0});
                    }
                    // translate camera
                    if (e.key.keysym.sym == SDLK_j) { // left
                        auto m = translate_x_plus * world_to_camera_rotate.transpose();
                        //cout << "m\n" << m << endl;
                        world_to_camera_translate.add_to_row(3, m.get_row(3));
                        world_to_camera_translate.set(3, 3, 1.0);
                        //cout << "world_to_camera_translate\n" << world_to_camera_translate << endl;
                    }
                    if (e.key.keysym.sym == SDLK_l) { // right
                        //world_to_camera_translate.add_to_row( 3, Vec{-0.08,0.0,0.0,0.0});
                        auto m = translate_x_minus * world_to_camera_rotate.transpose();

                        //auto m = translate_x_plus * world_to_camera_rotate;
                        world_to_camera_translate.add_to_row(3, m.get_row(3));
                        world_to_camera_translate.set(3, 3, 1.0);
                        //cout << "world_to_camera_translate: \n" 
                        //        << world_to_camera_translate << endl;
                    }
                    if (e.key.keysym.sym == SDLK_i) { // forward
                        //world_to_camera_translate.add_to_row( 3, Vec{0.0,0.0,-0.16,0.0});
                        auto m = translate_z_plus * world_to_camera_rotate.transpose();
                        //cout << "m\n" << m << endl;
                        world_to_camera_translate.add_to_row(3, m.get_row(3));
                        world_to_camera_translate.set(3, 3, 1.0);
                        //cout << "world_to_camera_translate\n" << world_to_camera_translate << endl;
                    }
                    if (e.key.keysym.sym == SDLK_k) { // back
                        //world_to_camera_translate.add_to_row( 3, Vec{0.0,0.0,0.16,0.0});
                        auto m = translate_z_minus * world_to_camera_rotate.transpose();
                        //cout << "m\n" << m << endl;
                        world_to_camera_translate.add_to_row(3, m.get_row(3));
                        world_to_camera_translate.set(3, 3, 1.0);
                        //cout << "world_to_camera_translate\n" << world_to_camera_translate << endl;
                    }
                    if (e.key.keysym.sym == SDLK_n) {
                        world_to_camera_translate.add_to_row( 
                                3, Vec{0.0,0.16,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_m) {
                        world_to_camera_translate.add_to_row( 
                                3, Vec{0.0,-0.16,0.0,0.0});
                    }
                    // rotate camera
                    if (e.key.keysym.sym == SDLK_u) {
                        world_to_camera_rotate = world_to_camera_rotate * turn_left;
                        //Vec<double> translation = world_to_camera.get_row(3);
                        //world_to_camera.zero_row(3);
                        //world_to_camera = turn_left * world_to_camera;
                        //world_to_camera.add_to_row(3, translation);
                    }
                    if (e.key.keysym.sym == SDLK_o) {
                        world_to_camera_rotate = world_to_camera_rotate * turn_right;
                        //Vec<double> translation = world_to_camera.get_row(3);
                        //world_to_camera.zero_row(3);
                        //world_to_camera = turn_right * world_to_camera;
                        //world_to_camera.add_to_row(3, translation);
                    }
                    // transform (scale) selected model 
                    if (e.key.keysym.sym == SDLK_c) {
                        model_scale[control_select] = model_scale[control_select] * scale_shrink;
                    }
                    if (e.key.keysym.sym == SDLK_v) {
                        model_scale[control_select] = model_scale[control_select] * scale_enlarge;
                    }
                    // print debug info
                    if (e.key.keysym.sym == SDLK_p) {
                        //// print object in world coordinates
                        int model_index = 0;

                        cout << "start debug" << endl;
                        print_info = true;
                        //cout << "control_select " << control_select << endl;
                        //model_transforms[control_select] = model_scale[control_select] 
                        //    * model_rotate[control_select] * model_translate[control_select];
                        //cout << "model_transform" << endl;
                        //cout << model_transforms[control_select] << endl;
                        //cout << "model_scale" << endl;
                        //cout << model_scale[control_select] << endl;
                        //cout << "model_rotate" << endl;
                        //cout << model_rotate[control_select] << endl;
                        //cout << "model_translate" << endl;
                        //cout << model_translate[control_select] << endl;
                        //cout << "scale_shrink" << endl;
                        //cout << scale_shrink<< endl;
                        cout << "lines_outputs_clipped_count = " << lines_outputs_clipped_count<< endl;
                        //for (int i = 0; i < lines_outputs_clipped_count; i++) {
                        //    cout << "clip_lines " << lines_outputs_clipped[i][0] << ", " 
                        //            << lines_outputs_clipped[i][1] << endl;
                        //}
                        print_clip_vertices = true;

                        //cout << "world_to_camera" << endl;
                        //cout << world_to_camera << endl;
                        //cout << "world_to_camera_translate" << endl;
                        //cout << world_to_camera_translate << endl;
                        //cout << "world_to_camera_rotate" << endl;
                        //cout << world_to_camera_rotate << endl;
                        //cout << "turn_left" << endl;
                        //cout << turn_left << endl;
                        cout << "clip_matrix" << endl;
                        cout << clip_matrix << endl;

                    }
                    break;
                default:
                    break;
            }
        }
        SDL_SetRenderDrawColor(gsdl.renderer, 20, 20, 20, 20);
        SDL_RenderClear(gsdl.renderer);
        SDL_SetRenderDrawColor(gsdl.renderer, 120, 120, 120, 120);
        //SDL_RenderDrawLine(gsdl.renderer, 400, 400, 800, 800);

        //
        // rendering pipeline
        //

        //
        //model to world transform
        //
        if (!posted_info) { cout << "world space vertices" << endl; }
        int model_index = 0;
        for (auto& output : model_outputs) {
            // compute model_transform for given translation, rotation and scale
            model_transforms[model_index] = model_scale[model_index] 
                * model_rotate[model_index] * model_translate[model_index];
            for (int i = 0; i < output.size(); ++i) {
                output[i] = models[model_index][i] * model_transforms[model_index];
            }
            model_index++;
        }
        //
        //world to camera (triangles)
        //
        if (!posted_info) { cout << "camera space vertices" << endl; }
        model_index = 0;
        // translate first!
        auto v = -1.0 * world_to_camera_translate.get_row(3);
        v.data[3] = 1.0;
        // old way below
        world_to_camera = world_to_camera_translate * world_to_camera_rotate;
        // new way below
        //world_to_camera = world_to_camera_translate.transpose() 
        //        * world_to_camera_rotate.transpose();
        for (auto& model : model_outputs) {
            for (int i = 0; i < model.size(); ++i) {
                model[i] = model[i] * world_to_camera;
            }
            model_index++;
        }
        //
        //world to camera (lines)
        //
        model_index = 0;
        world_to_camera = world_to_camera_translate * world_to_camera_rotate;
        for (auto& output : lines_outputs) {
            for (int i = 0; i < output.size(); ++i) {
                output[i] = lines_world[model_index][i] * world_to_camera; 
            }
            model_index++;
        }
        if (print_info) {
            cout << "printing world coordinates of lines\n";
            for (auto& output : lines_world) {
                cout << output[0] << ", " << output[1] << endl;
            }
            cout << "printing view/camera coordinates of lines\n";
            for (auto& output : lines_outputs) {
                cout << output[0] << ", " << output[1] << endl;
            }
        }
       
        //
        //camera to clip (triangles)
        //
        model_index = 0;
        for (auto& model : model_outputs) {
            for (int i = 0; i < model.size(); ++i) {
                model[i] = model[i]  * clip_matrix * ortho_matrix;
            }
            model_index++;
        }
        //
        //camera to clip (lines)
        //
        model_index = 0;
        for (auto& output : lines_outputs) {
            for (int i = 0; i < output.size(); ++i) {
                output[i] = output[i] * clip_matrix * ortho_matrix;
            }
            model_index++;
        }
        if (print_info) {
            cout << "printing clip coordinates of lines\n";
            for (auto& output : lines_outputs) {
                cout << output[0] << ", " << output[1] << endl;
            }
        }
        //if (print_clip_vertices) {
        //    print_clip_vertices = false;
        //    for(const auto& output : lines_outputs) {
        //        cout << output[0] << ", " << output[1] << endl;
        //    }
        //}
        //// homogenous divide (triangles)
        //if (!posted_info) { cout << "homogenous divide vertices" << endl; }
        //model_index = 0;
        //for (auto& model : model_outputs) {
        //    for (int i = 0; i < model.size(); ++i) {
        //        double w = model[i].data[3];
        //        if (w != 0.0) {
        //            model[i] = (1.0/w) * model[i]; 
        //        } else {
        //            model[i].data[0] = 0.0; 
        //            model[i].data[1] = 0.0; 
        //            model[i].data[2] = 0.0; 
        //        }
        //    }
        //    //model_index++;
        //}
        //// homogenous divide (lines)
        //if (!posted_info) { cout << "homogenous divide vertices" << endl; }
        //model_index = 0;
        //for (auto& output : lines_outputs) {
        //    for (int i = 0; i < output.size(); ++i) {
        //        double w = output[i].data[3];
        //        if (w != 0.0) {
        //            output[i] = (1.0/w) * output[i]; 
        //        } else {
        //            output[i].data[0] = 0.0; 
        //            output[i].data[1] = 0.0; 
        //            output[i].data[2] = 0.0; 
        //        }
        //    }
        //    //model_index++;
        //}

        // 
        // clip out of bounds vertices (lines)
        //
        //lines_outputs_clipped = lines_outputs;
        model_index = 0;
        // clear previously clipped vertices
        for (auto& line_clipped : lines_outputs_clipped)
        {
            line_clipped.clear();
        }
        //array<vector<Vec<double>>, NUM_LINES> lines_outputs; 
        for (const vector<Vec<double>>& output : lines_outputs) {
        //for (const auto& output : lines_outputs) {
          
            //clip_plane = Vec<double>{
            //if (output[0].data[3] != output[1].data[3]) {
            //    cout << "asdf, " << output[0].data[3] << ", " << output[1].data[3] << endl;
            //    assert(output[0].data[3] == output[1].data[3]);
            //}

            clip_line3d(output, clip_planes, Vec{0.0, 0.0, 0.0}, 
                    lines_outputs_clipped[model_index]);
            int vertices_returned = lines_outputs_clipped[model_index].size(); 
            if (vertices_returned > 2) {
                assert(" line clip return's > 2 vertices\n");
            } else if (vertices_returned == 2) {
                model_index++;
            } else if (vertices_returned == 1) {
                assert(" line clip return's 1 vertices\n");
            } else if (vertices_returned == 0) {
                // fine, but don't increment model_index;
            } else {
                assert(" line clip returns wrong num of vertices\n");
            }
        }
        lines_outputs_clipped_count = model_index;

        //
        // homogenous divide (triangles)
        //
        if (!posted_info) { cout << "homogenous divide vertices" << endl; }
        model_index = 0;
        for (auto& model : model_outputs) {
            for (int i = 0; i < model.size(); ++i) {
                double w = model[i].data[3];
                if (w != 0.0) {
                    model[i] = (1.0/w) * model[i]; 
                } else {
                    model[i].data[0] = 0.0; 
                    model[i].data[1] = 0.0; 
                    model[i].data[2] = 0.0; 
                }
            }
            //model_index++;
        }
        //
        // homogenous divide (lines)
        //
        if (!posted_info) { cout << "homogenous divide vertices" << endl; }
        model_index = 0;
        for (auto& output : lines_outputs_clipped) {
            for (int i = 0; i < output.size(); ++i) {
                double w = output[i].data[3];
                if (w != 0.0) {
                    output[i] = (1.0/w) * output[i]; 
                } else {
                    output[i].data[0] = 0.0; 
                    output[i].data[1] = 0.0; 
                    output[i].data[2] = 0.0; 
                }
            }
            //model_index++;
        }
        //
        // project to screen (triangles)
        //
        if (!posted_info) { cout << "screen vertices" << endl; }
        for (auto& model : model_outputs) {
            for (int i = 0; i < model.size(); ++i) {
                model[i].data[0] = model[i].data[0] * win_res_x
                    / (model[i].data[2] * 2.0) + 1000.0/2.0;
                model[i].data[1] = -1.0 * model[i].data[1] * win_res_y
                    / (model[i].data[2] * 2.0) + 1000.0/2.0;
            }
        }
        //
        // project to screen (lines)
        //
        if (!posted_info) { cout << "screen vertices" << endl; }
        for (auto& output : lines_outputs_clipped) {
            for (int i = 0; i < output.size(); ++i) {
                output[i].data[0] = output[i].data[0] * win_res_x
                    / (output[i].data[2] * 2.0) + 1000.0/2.0;
                output[i].data[1] = -1.0 * output[i].data[1] * win_res_y
                    / (output[i].data[2] * 2.0) + 1000.0/2.0;
            }
        }

        //
        // render lines
        //
        // const int origin_index = 2;
        // bool dont_render_origin = false;
        if (RENDER_LINES) {
        int i = 0;
        for (auto& output : lines_outputs_clipped) {
            if(i < lines_outputs_clipped_count) {
                i++;
            } else {
                break;
            }
            //cout << " allo2" << endl;
            //cout << "output.size() == " << output.size() << endl;
            //cout << "lines_outputs.size() == " << lines_outputs.size() << endl;
            for (int i_line = 0; i_line < output.size(); i_line += 2) {
                SDL_SetRenderDrawColor(gsdl.renderer, 50, 100, 50, 20);
                // kludge, if any vertex out of bounds, break out of this drawing iteration
                //for (int i = i_tri; i < i_tri + 3; ++i) {
                //    if (model[i].data[0] <= 0 || model[i].data[0] >= win_res_x) {
                //        goto dont_draw;   
                //    }
                //    if (model[i].data[1] <= 0 || model[i].data[1] >= win_res_y) {
                //        goto dont_draw;
                //    }
                //}
                //if (model_index != origin_index) {
                //    //draw_triangle(model[i_tri], model[i_tri + 1], model[i_tri + 2]);
                //}
                // i_tri is index to begining of triangle
                //for (int i = i_tri; i < i_tri + 3; ++i) {
                //    // draws wireframe
                //    // i at vertex of triangle
                //    int j = (i + 1 < i_tri + 3) ? i + 1 : i_tri;
                //    //SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
                //    SDL_RenderDrawLine(gsdl.renderer, 
                //            output[i].data[0], output[i].data[1],
                //            output[j].data[0], output[j].data[1]);
                //}
                SDL_RenderDrawLine(gsdl.renderer, 
                        output[i_line].data[0], output[i_line].data[1],
                        output[i_line+1].data[0], output[i_line+1].data[1]);
                //cout << "allo" << endl;
            }
            dont_draw_line:;
            model_index++;
        } 
        }

        //
        // render triangles ( as wireframe or filled in)
        //
        model_index = 0;
        const int origin_index = 2;
        bool dont_render_origin = false;
        for (auto& model : model_outputs) {
            for (int i_tri = 0; i_tri < model.size(); i_tri += 3) {
                //Vec<double> triangle[3];
                SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
                if (model_index == origin_index && dont_render_origin) {
                    goto dont_draw;
                }
                // kludge, if any vertex out of bounds, break out of this drawing iteration
                for (int i = i_tri; i < i_tri + 3; ++i) {
                    if (model[i].data[0] <= 0 || model[i].data[0] >= win_res_x) {
                        goto dont_draw;   
                    }
                    if (model[i].data[1] <= 0 || model[i].data[1] >= win_res_y) {
                        goto dont_draw;
                    }
                }
                if (model_index != origin_index) {
                    //draw_triangle(model[i_tri], model[i_tri + 1], model[i_tri + 2]);
                }
                // i_tri is index to begining of triangle
                for (int i = i_tri; i < i_tri + 3; ++i) {
                    // draws wireframe
                    // i at vertex of triangle
                    int j = (i + 1 < i_tri + 3) ? i + 1 : i_tri;
                    //SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
                    if (model_index == origin_index) {
                        // give origin marker a special color
                        if (i_tri == 0) {
                            // x axis red
                            SDL_SetRenderDrawColor(gsdl.renderer, 250, 0, 0, 20);
                        } else if (i_tri == 3) {
                            // y axis green
                            SDL_SetRenderDrawColor(gsdl.renderer, 0, 250, 0, 20);
                        } else {
                            // z axis blue
                            SDL_SetRenderDrawColor(gsdl.renderer, 0, 0, 250, 20); 
                        }
                    } else {
                        SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
                    }
                    SDL_RenderDrawLine(gsdl.renderer, 
                            model[i].data[0], model[i].data[1],
                            model[j].data[0], model[j].data[1]);
                }
            }
            dont_draw:;
            model_index++;
        }


        if (posted_info == false) {
            posted_info = true;
        }
        print_info = false;


        SDL_RenderPresent(gsdl.renderer);
        time_frame_duration = SDL_GetTicks64() - time_frame_start;
        if (time_frame_duration < 1000/30) {
            SDL_Delay(1000/30 - time_frame_duration);
        }
    }
    SDL_DestroyWindow(gsdl.window);
    SDL_Quit();
}

