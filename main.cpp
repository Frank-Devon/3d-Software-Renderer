#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <exception>
#include <string>
#include <sstream>
#include <ios>
#include <istream>
#include <algorithm>

#include "vec.hpp"
#include "test.hpp"
#include "sdl_init.hpp"

//ctags -R --languages=C++ --extra=+fq --exclude=.git .
using std::cout, std::endl, std::vector, std::array, std::string;

// vertexes of object
#define NUM_VERTICES 9
#define NUM_MESHES 3
#define ROTATE_FIRST 0
#define RENDER_LINES 0  // renders grid lines
#define NUM_LINES 22
bool BACK_FACE_CULLING = true;

Vec<double> model_0_pos;
Vec<double> model_1_pos;

//array<Model<double>, 3> models = { 
vector<Model<double>> models = { 
    Model {
        {
            // Front face (2 triangles, 6 vertices)
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            
            // Back face (2 triangles, 6 vertices)
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            Vec{-1.0,  1.0,  1.0, 1.0},  // Back top left
            
            // Left face (2 triangles, 6 vertices)
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            Vec{-1.0,  1.0,  1.0, 1.0},  // Back top left
            
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{-1.0,  1.0,  1.0, 1.0},  // Back top left
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            
            // Right face (2 triangles, 6 vertices)
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            
            // Bottom face (2 triangles, 6 vertices)
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            
            // Top face (2 triangles, 6 vertices)
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            Vec{-1.0,  1.0,  1.0, 1.0}   // Back top left
            
        },  
        Color{ 20, 0, 0, 200},
        true
    },

    Model {
        {
            // Front face (2 triangles, 6 vertices)
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            
            // Back face (2 triangles, 6 vertices)
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            Vec{-1.0,  1.0,  1.0, 1.0},  // Back top left
            
            // Left face (2 triangles, 6 vertices)
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            Vec{-1.0,  1.0,  1.0, 1.0},  // Back top left
            
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{-1.0,  1.0,  1.0, 1.0},  // Back top left
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            
            // Right face (2 triangles, 6 vertices)
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            
            // Bottom face (2 triangles, 6 vertices)
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{ 1.0, -1.0, -1.0, 1.0},  // Front bottom right
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            
            Vec{-1.0, -1.0, -1.0, 1.0},  // Front bottom left
            Vec{ 1.0, -1.0,  1.0, 1.0},  // Back bottom right
            Vec{-1.0, -1.0,  1.0, 1.0},  // Back bottom left
            
            // Top face (2 triangles, 6 vertices)
            Vec{ 1.0,  1.0, -1.0, 1.0},  // Front top right
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            
            Vec{ 1.0,  1.0,  1.0, 1.0},  // Back top right
            Vec{-1.0,  1.0, -1.0, 1.0},  // Front top left
            Vec{-1.0,  1.0,  1.0, 1.0}   // Back top left
            
        },  
        Color{ 0, 20, 0, 200},
        true
    },

    Model {
        {
            // x axis
            Vec{0.0, 0.0, -0.1, 1.0 },
            Vec{2.0, 0.0, 0.0, 1.0 },
            Vec{0.0, 0.0, 0.1, 1.0 },
            
            // y axis
            Vec{0.0, 0.0, -0.1, 1.0 },
            Vec{0.0, 2.0, 0.0, 1.0 },
            Vec{0.0, 0.0, 0.1, 1.0 },
            
            // z axis
            Vec{-0.1, 0.0, 0.0, 1.0 },
            Vec{0.0, 0.0, 2.0, 1.0 },
            Vec{0.1, 0.0, 0.0, 1.0 }
            
        },  
        Color{ 200, 200, 200, 200},
        true,
        true
    },
};

// vertices drawn to screen
array<vector<Vec<double>>, NUM_MESHES> model_outputs = // TODO DELETE
{
    vector<Vec<double>>{},
    vector<Vec<double>>{},
    vector<Vec<double>>{}
};

auto models_world = models; // world space of models 
auto models_culled = models; // world space of models with back face culling applied
auto models_view = models; // view space of models
auto models_clipped = models; // clip space of models
auto models_ndc = models;  // normalized device coordinates of models
auto models_screen = models;  // screen space
auto models_out = models;  // perhaps could use a single variable instead of all the above

// basic lighting
Vec<double> light_dir = {-1.0, -1.0, 1.0};

array<vector<Vec<double>>, NUM_MESHES> model_outputs_clipped;

array<vector<Vec<double>>, NUM_LINES> lines_world;
array<vector<Vec<double>>, NUM_LINES> lines_outputs; 
array<vector<Vec<double>>, NUM_LINES> lines_outputs_clipped;
int lines_outputs_clipped_count = 0;
Model<double>* model_selected = nullptr; // model selected for transformations

int win_res_x = 1000;
int win_res_y = 1000;
// meshes
//Control_Select control_select = Control_Select::model_0;
enum Control_Select : unsigned int { m_0 = 0, m_1 };
Control_Select control_select = m_0; // TODO delete
bool print_clip_vertices = false; // TODO delete?
bool print_info = false;
std::array<std::array<Vec<double>, NUM_VERTICES>, NUM_MESHES> meshes; 


double cam_yaw = 0.0;
double cam_pitch = 0.0;

// important! determines yaw
Mat<double> cam_yaw_matrix = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};

// important! determines pitch
Mat<double> cam_pitch_matrix = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};

// should really be called point_at matrix
Mat<double> world_to_camera = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};

// still used to hold position
Mat<double> world_to_camera_translate = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 3.0, -10.0, 1.0}
};

// depreciated
Mat<double> world_to_camera_rotate = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.00, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};

// the ACTUAL world_to_camera matrix
Mat<double> look_at = {
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};

//
// probably overkill making these matrices
//
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


//
// user input variables
//
bool middle_mouse_down = false;
Vec<double> middle_mouse_down_pos_last_frame = {0, 0};
Vec<double> middle_mouse_down_delta = {0, 0};

//
// user input variables END
//

Mat<double> ortho_matrix(4, 4);
Mat<double> clip_matrix(4, 4); // not used yet, to be used for cliping vertices outside view frustrum.
Mat<double> turn_left(4, 4);          
Mat<double> turn_right(4, 4);          
Mat<double> turn_down(4, 4);          
Mat<double> turn_up(4, 4);          


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

// TODO Raster the triangle from by calling draw pixel function
// render triangle between 3 points
template<typename T>
//void draw_triangle(const Vec<T>& a, const Vec<T>& b, const Vec<T>& c)
void draw_triangle(const Tri<T>& tri) 
{
    //    typedef struct SDL_Color
    //    {
    //        Uint8 r;
    //        Uint8 g;
    //        Uint8 b;
    //        Uint8 a;
    //    } SDL_Color;
    //SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
    SDL_Vertex verts[3];
    SDL_Color color = { tri.color.r, tri.color.g, tri.color.b, tri.color.a };
    verts[0].position.x = tri.verts[0].data[0];
    verts[0].position.y = tri.verts[0].data[1];
    verts[0].color = color;
    verts[1].position.x = tri.verts[1].data[0];
    verts[1].position.y = tri.verts[1].data[1];
    verts[1].color = color;
    verts[2].position.x = tri.verts[2].data[0];
    verts[2].position.y = tri.verts[2].data[1];
    verts[2].color = color;
    SDL_RenderGeometry( gsdl.renderer, nullptr, &verts[0], 3, nullptr, 0 );
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

void draw_text(string str, int x, int y) {
    SDL_Color White = {255, 255, 255};
    SDL_Surface* surfaceMessage =
        TTF_RenderText_Solid(gsdl.font, str.c_str(), White);
    // now you can convert it into a texture
    SDL_Texture* Message = SDL_CreateTextureFromSurface(gsdl.renderer, surfaceMessage);
    SDL_Rect Message_rect; //create a rect
    Message_rect.x = x;  //controls the rect's x coordinate
    Message_rect.y = y; // controls the rect's y coordinte
    Message_rect.w = str.length() * 20; // controls the width of the rect
    Message_rect.h = 20; // controls the height of the rect
    SDL_RenderCopy(gsdl.renderer, Message, NULL, &Message_rect);
    
    SDL_FreeSurface(surfaceMessage);
    SDL_DestroyTexture(Message);
}


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

template<typename T>
Mat<T> look_at_matrix(Vec<T>& _pos, Vec<T>& _target, Vec<T>& _up) 
{
    assert(_pos.length == 3 && _target.length == 3 && _up.length == 3);
    // new forward
    Vec<T> forward = _target - _pos;
    forward = forward.unit();
    // new up
    Vec<T> a = forward.dot(_up) * forward; 
    Vec<T> up = _up - a;
    up = up.unit();
    // new right
    Vec<T> right = up.cross(forward);
    right = right.unit();
    // construct 3x3 R matrix
    Mat<T> R(3, 3);
    R.set_row(0, right);
    R.set_row(1, up);
    R.set_row(2, forward);
    //R.set(3, 3, 1.0);
    R = R.transpose(); // R is now inverted
    
    Vec<T> t = (- 1.0 * _pos) * R;  // translation portion

    // construct matrix
    Mat<T> m(4, 4);
    m.set(0, 0, R.get(0, 0));
    m.set(0, 1, R.get(0, 1));
    m.set(0, 2, R.get(0, 2));
    m.set(1, 0, R.get(1, 0));
    m.set(1, 1, R.get(1, 1));
    m.set(1, 2, R.get(1, 2));
    m.set(2, 0, R.get(2, 0));
    m.set(2, 1, R.get(2, 1));
    m.set(2, 2, R.get(2, 2));

    //set translation
    m.set(3, 0, t.data[0]);
    m.set(3, 1, t.data[1]);
    m.set(3, 2, t.data[2]);
    m.set(3, 3, 1.0);
    return m;
}


void rotate( double yaw, double pitch, Mat<double>& m)
{
    float cosYaw = cos(yaw);
    float sinYaw = sin(yaw);
    float cosPitch = cos(pitch);
    float sinPitch = sin(pitch);
   
    // first number is row, second is column
    m.set(0, 0, cosYaw); // m[0][0] = cosYaw;
    m.set(0, 1, -sinYaw * cosPitch); //   m[0][1] = -sinYaw * cosPitch;
    m.set(0, 2, sinYaw * sinPitch); // m[0][2] = sinYaw * sinPitch;
    
    m.set(1, 0, sinYaw); //m[1][0] = sinYaw;
    m.set(1, 1, cosYaw * cosPitch);  //m[1][1] = cosYaw * cosPitch;
    m.set(1, 2, -cosYaw * sinPitch); //m[1][2] = -cosYaw * sinPitch;
    
    m.set(2, 1, sinPitch); //m[2][1] = sinPitch;
    m.set(2, 1, cosPitch); //m[2][2] = cosPitch;
    m.set(3, 3, 1.0); //m[2][2] = cosPitch;
}

void test_lerp( ) 
{
    // test lerp
    Vec<double> a = { 1.0, 1.0 };
    Vec<double> b = { 5.0, 3.0 };
    cout << "testing lerp\n";
    cout << a.lerp(b, 0.0) << endl << a.lerp(b, 0.5) << endl << a.lerp(b, 1.0) << endl;
}

int main() 
{
    //test_matrix();
    light_dir = light_dir.unit(); // normalize light direction
    models[0].trans.set_row(3, Vec<double>{-9.0, 7.0, 17.0, 1.0});
    models[1].trans.set_row(3, Vec<double>{-2.0, 11.5, 19.0, 1.0});
    world_to_camera_translate.set_row(3, Vec<double>{0.0, 17.0, 2.2, 1.0});
    models.reserve(20);
    bool load_success = models[0].load_from_file("VideoShip.obj", {0, 0, 1, 1}, false);
    load_success = models[1].load_from_file("teapot.obj", {0, 0, 1, 1}, false);
    load_success = models[2].load_from_file("mountains.obj", {0, 0, 1, 1}, false);
    int index = 0;
    construct_grid();
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
    turn_right.set(0, 0, cos(rad));
    turn_right.set(0, 2, -sin(rad));
    turn_right.set(1, 1, 1.0);
    turn_right.set(2, 0, sin(rad));
    turn_right.set(2, 2, cos(rad));
    turn_right.set(3, 3, 1.0);
    rad = - rad;
    turn_left.set(0, 0, cos(rad));
    turn_left.set(0, 2, -sin(rad));
    turn_left.set(1, 1, 1.0);
    turn_left.set(2, 0, sin(rad));
    turn_left.set(2, 2, cos(rad));
    turn_left.set(3, 3, 1.0);

    rad = - rad;
    turn_up.set(0, 0, 1.0);
    turn_up.set(1, 1, cos(rad));
    turn_up.set(2, 1, -sin(rad));
    turn_up.set(1, 2, sin(rad));
    turn_up.set(2, 2, cos(rad));
    turn_up.set(3, 3, 1.0);
    rad = - rad;
    turn_down.set(0, 0, 1.0);
    turn_down.set(1, 1, cos(rad));
    turn_down.set(2, 1, -sin(rad));
    turn_down.set(1, 2, sin(rad));
    turn_down.set(2, 2, cos(rad));
    turn_down.set(3, 3, 1.0);

    auto test0 = turn_left * turn_left * turn_left;
    cout << "test shiet" << endl << test0 << endl;
    // camera space to clip space
    //double near = 1.00;  // 1.0


    //
    //  calculating orthographic and prespective matrices
    //
    double near = 1.00;  // 1.0
    double far = 20000.0;//200000.0;  // 20.0
    double projection_width = 0.6;
    double projection_height = 0.6;
    clip_matrix.set(0, 0,  near );
    clip_matrix.set(1, 1,  near );
    //clip_matrix.set(0, 0, near );
    //clip_matrix.set(1, 1, near );
    clip_matrix.set(2, 2, (near + far));
    clip_matrix.set(2, 3, 1.0 );
    clip_matrix.set(3, 2, -2.0*( near * far) ); 
    //clip_matrix.set(0, 0, near / (projection_width * 2));
    //clip_matrix.set(1, 1, near / (projection_height * 2));
    //clip_matrix.set(2, 2, (near + far)/(far - near));
    //clip_matrix.set(2, 3, 1.0 );
    //clip_matrix.set(3, 2, (- 2 * near * far) / (far - near));
    cout << "here's the clip matrix" << endl;
    cout << clip_matrix << endl;
    // ortho matrix set up
    //double r = 1.0, l = -1.0, t = 1.0, b = -1.0;
    double r = 0.3, l = -0.3, t = 0.3, b = -0.3;
    ortho_matrix.set(0, 0, 2.0/(r-l));
    ortho_matrix.set(3, 0, - (r+l)/(r-l));
    ortho_matrix.set(1, 1, 2.0/(t-b));
    ortho_matrix.set(3, 1, - (t+b)/(t-b));
    //ortho_matrix.set(2, 2, 2.0/(near-far));
    //ortho_matrix.set(3, 2, - (near+far)/(near-far));
    ortho_matrix.set(2, 2, - 2.0/(near-far));
    ortho_matrix.set(3, 2, - (near+far)/(near-far));
    ortho_matrix.set(3, 3, 1.0);
    //ortho_matrix.set_identity();



    using std::cout, std::endl;
    //getchar();
    bool exit = false;    
    cout << "gsdl status: " << gsdl.created_ok << endl;
    SDL_Event e;
    Uint64 time_frame_start;
    Uint64 time_frame_duration = 0;
    model_selected = &models[0]; // model selected for transformations

    bool posted_info = true;
    while (!exit) {
        time_frame_start = SDL_GetTicks64();
        while (SDL_PollEvent(&e) != 0) {
            switch (e.type) {
                case SDL_QUIT:
                    exit = true;
                    break;
                case SDL_MOUSEBUTTONDOWN:
                    if (e.button.button == SDL_BUTTON_MIDDLE) {
                        //std::cout << "Middle mouse button pressed down." << std::endl;
                        middle_mouse_down = true;
                        int x, y;
                        SDL_GetMouseState(&x, &y);
                        middle_mouse_down_pos_last_frame.data[0] = x;
                        middle_mouse_down_pos_last_frame.data[1] = y;
                    }
                    break;
                case SDL_MOUSEBUTTONUP:
                    if (e.button.button == SDL_BUTTON_MIDDLE) {
                        std::cout << "Middle mouse button released." << std::endl;
                        middle_mouse_down = false;

                    }
                    break;
                case SDL_KEYDOWN:
                    // rotate selected object
                    if (e.key.keysym.sym == SDLK_q) {
                        model_selected->rot = model_selected->rot * turn_left;
                    }
                    if (e.key.keysym.sym == SDLK_e) {
                        model_selected->rot = model_selected->rot * turn_right;
                    }
                    if (e.key.keysym.sym == SDLK_r) {
                        model_selected->rot = model_selected->rot * turn_up;
                    }
                    if (e.key.keysym.sym == SDLK_f) {
                        model_selected->rot = model_selected->rot * turn_down;
                    }

                    if (e.key.keysym.sym == SDLK_1) {
                        cout << "1 pressed" << endl;
                        model_selected = &models[0];
                    }
                    if (e.key.keysym.sym == SDLK_2) {
                        cout << "2 pressed" << endl;
                        model_selected = &models[1];
                    }
                    if (e.key.keysym.sym == SDLK_3) {
                        cout << "3 pressed" << endl;
                        if (models.size() > 3) { // exclude index 3, the axis!
                            model_selected = &models[3];
                        }
                    }
                    if (e.key.keysym.sym == SDLK_4) {
                        cout << "4 pressed" << endl;
                        if (models.size() > 4) { // exclude index 3, the axis!
                            model_selected = &models[4];
                        }
                    }
                    // transform selected model
                    if (e.key.keysym.sym == SDLK_a) { // left
                        model_selected->trans.add_to_row(3, Vec{-0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_d) { // right
                        model_selected->trans.add_to_row(3, Vec{0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_w) { // forward
                        model_selected->trans.add_to_row(3, Vec{0.00,0.0,0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_s) { // backward
                        model_selected->trans.add_to_row(3, Vec{0.00,0.0,-0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_z) {  // below
                        model_selected->trans.add_to_row(3, Vec{0.0,-0.08,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_x) { // above
                        model_selected->trans.add_to_row(3, Vec{0.0,0.08,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_c) { // shrink
                        model_selected->scale = model_selected->scale * scale_shrink;
                    }
                    if (e.key.keysym.sym == SDLK_v) { // enlarge
                        model_selected->scale = model_selected->scale * scale_enlarge;
                    }
                    // transform camera
                    if (e.key.keysym.sym == SDLK_j) { // left
                        //auto m = translate_x_minus * world_to_camera_rotate;
                        //world_to_camera_translate.add_to_row(3, m.get_row(3));
                        //world_to_camera_translate.set(3, 3, 1.0);
                        ////world_to_camera_translate.add_to_row( 
                        ////        3, Vec{-0.08,0.0,0.0,0.0});
                        auto m = translate_x_minus * rotate_y(cam_yaw);
                        world_to_camera_translate.add_to_row( 3, m.get_row(3));
                        world_to_camera_translate.set(3, 3, 1.0);
                    }
                    if (e.key.keysym.sym == SDLK_l) { // right
                        auto m = translate_x_plus * rotate_y(cam_yaw);
                        world_to_camera_translate.add_to_row( 3, m.get_row(3));
                        world_to_camera_translate.set(3, 3, 1.0);
                    }
                    if (e.key.keysym.sym == SDLK_i) { // forward
                        Vec<double> look_dir = {0.0, 0.0, 1.0, 1.0}; 
                        look_dir = look_dir * rotate_x(cam_pitch) * rotate_y(cam_yaw);
                        auto t = 0.08 * look_dir;
                        world_to_camera_translate.add_to_row(3, t);
                        world_to_camera_translate.set(3, 3, 1.0);
                    }
                    if (e.key.keysym.sym == SDLK_k) { // back
                        Vec<double> look_dir = {0.0, 0.0, 1.0, 1.0}; 
                        look_dir = look_dir * rotate_x(cam_pitch) * rotate_y(cam_yaw);
                        auto t = -0.08 * look_dir;
                        world_to_camera_translate.add_to_row(3, t);
                        world_to_camera_translate.set(3, 3, 1.0);
                    }
                    if (e.key.keysym.sym == SDLK_m) { // up
                        world_to_camera_translate.add_to_row( 
                                3, Vec{0.0,0.16,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_n) { // down
                        world_to_camera_translate.add_to_row( 
                                3, Vec{0.0,-0.16,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_u) { // rotate left
                        //world_to_camera_rotate = world_to_camera_rotate * turn_left;
                        //look_dir = look_dir * turn_left;
                        cam_yaw -= 0.035;
                        //cam_yaw_matrix = cam_yaw_matrix * turn_left;
                    }
                    if (e.key.keysym.sym == SDLK_o) { // rotate right
                        cam_yaw += 0.035;
                    }
                    if (e.key.keysym.sym == SDLK_LEFTBRACKET) { // rotate camera up 
                        //world_to_camera_rotate = world_to_camera_rotate * turn_up;
                        //look_dir = look_dir * turn_up;
                        cam_pitch -= 0.035;
                        //cam_pitch_matrix = cam_pitch_matrix * turn_up;
                    }
                    if (e.key.keysym.sym == SDLK_RIGHTBRACKET) { // rotate camera up 
                        cam_pitch += 0.035;
                    }
                    
                    // toggle back face culling
                    if (e.key.keysym.sym == SDLK_0) {
                        BACK_FACE_CULLING = !BACK_FACE_CULLING;
                    }
                    // add new ship
                    if (e.key.keysym.sym == SDLK_7) {
                        cout << "spawn new ship disabled\n";
                        //Model<double> m;
                        //static double z_pos = 4.0;
                        //m.trans.set(3, 2, z_pos);
                        //z_pos += 4.5;
                        //m.load_from_file("VideoShip.obj");
                        //models.push_back(m);
                    }
                    

                    // toggle drawing wire frames
                    if (e.key.keysym.sym == SDLK_8) {
                        for (auto& model : models) {
                            model.draw_wire_frame = !model.draw_wire_frame;
                        }
                    }
                    // toggle drawing faces
                    if (e.key.keysym.sym == SDLK_9) {
                        for (auto& model : models) {
                            model.draw_face = !model.draw_face;
                        }
                    }
                    // print debug info
                    if (e.key.keysym.sym == SDLK_p) {
                        //cout << "look_at matrix\n" << look_at << '\n';
                        cout << "models.size() " << models.size() << '\n';
                        cout << "printing look_at matrix\n" << look_at << '\n';
                        cout << "printing world_to_camera_rotate\n" << world_to_camera_rotate << '\n';
                        cout << "printing world_to_camera_translate\n" << world_to_camera_translate << '\n';
                    }
                    break;
                default:
                    break;
            }
        }
        
        if (middle_mouse_down) {
            // get current mouse position
            int x, y;
            SDL_GetMouseState(&x, &y);
            // calculate mouse delta
            middle_mouse_down_delta.data[0] = x - middle_mouse_down_pos_last_frame.data[0];
            middle_mouse_down_delta.data[1] = y - middle_mouse_down_pos_last_frame.data[1];
            // calculate new yaw
            cam_yaw = cam_yaw + middle_mouse_down_delta.data[0] * PI / 1000.0;
            // calculate new pitch
            cam_pitch = cam_pitch + middle_mouse_down_delta.data[1] * PI / 1000.0;
            // store last frame mouse pos
            middle_mouse_down_pos_last_frame.data[0] = x;
            middle_mouse_down_pos_last_frame.data[1] = y;
        } else {
            middle_mouse_down_pos_last_frame.data[0] = 0;
            middle_mouse_down_pos_last_frame.data[1] = 0;
        }


        SDL_SetRenderDrawColor(gsdl.renderer, 20, 20, 20, 20);
        SDL_RenderClear(gsdl.renderer);
        SDL_SetRenderDrawColor(gsdl.renderer, 120, 120, 120, 120);
        //SDL_RenderDrawLine(gsdl.renderer, 400, 400, 800, 800);

        //
        // rendering pipeline
        //
        int model_index = 0;

        //initialize output variable models_world
        models_world = models;
        models_out = models;

        //
        // model to world transform
        //
        for (int j = 0; j < models.size(); ++j) { // model index
            for (int i = 0; i < models[j].tris.size(); ++i) { // tri index
                models_world[j].tris[i].verts[0] = models[j].tris[i].verts[0] 
                        * models[j].scale * models[j].rot * models[j].trans;
                models_world[j].tris[i].verts[1] = models[j].tris[i].verts[1] 
                        * models[j].scale * models[j].rot * models[j].trans;
                models_world[j].tris[i].verts[2] = models[j].tris[i].verts[2] 
                        * models[j].scale * models[j].rot * models[j].trans;
                
                //compute lighting
                models_world[j].tris[i].color = models[j].color; // set initial color
                Vec a = (models_world[j].tris[i].verts[1] - models_world[j].tris[i].verts[0]).slice(0, 3);
                Vec b = (models_world[j].tris[i].verts[2] - models_world[j].tris[i].verts[0]).slice(0, 3);
                Vec c = -1.0 * a.cross(b);
                c = c.unit();
                double dp = c.dot( light_dir); 
                dp = dp * 255;
                int dp_int = int(dp); 
                dp_int = std::max( 0, dp_int);
                models_world[j].tris[i].color = models[j].color * dp_int ; // set initial color
                
                //misc data
                models_world[j].tris[i].draw_wire_frame = models[j].draw_wire_frame;
                models_world[j].tris[i].draw_face = models[j].draw_face;
            }
        }

        // get average position of model 0 for debug purposes
        model_0_pos = Vec<double>(4);
        for (int i = 0; i < models_world[0].tris.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                model_0_pos = model_0_pos + models_world[0].tris[i].verts[j];
            }
        }
        model_0_pos = model_0_pos / (models_world[0].tris.size() * 3); 
        // get average position of model 1 for debug purposes
        model_1_pos = Vec<double>(4);
        for (int i = 0; i < models_world[1].tris.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                model_1_pos = model_1_pos + models_world[1].tris[i].verts[j];
            }
        }
        model_1_pos = model_1_pos / (models_world[1].tris.size() * 3); 

        //
        // triangles: backface culling, in world space
        //
        
        models_culled = models_world;
        for (auto& m : models_culled ) {
            m.tris.clear();
        }
        for (int j = 0; j < models_world.size(); ++j) {
            for (int i = 0; i < models_world[j].tris.size(); ++i) {
                if (models_world[j].tris[i].draw_back_face) {
                    models_culled[j].tris.push_back( models_world[j].tris[i]);
                    continue;
                }
                // get normal of triangle
                Vec<double> line0 = (models_world[j].tris[i].verts[1] 
                        - models_world[j].tris[i].verts[0]).slice(0, 3);
                Vec<double> line1 = (models_world[j].tris[i].verts[2] 
                        - models_world[j].tris[i].verts[0]).slice(0, 3);
                Vec<double> norm = (line0.cross(line1)).unit();
                Vec<double> point = models_world[j].tris[i].verts[0].slice(0, 3); 
                Vec<double> point_cam = (world_to_camera_translate.get_row(3)).slice(0, 3);
                point = point - point_cam;
                double dp = norm.dot(point); 
                if (!BACK_FACE_CULLING || dp < 0.0) {
                    models_culled[j].tris.push_back( models_world[j].tris[i]);
                }
            }
        }

        //
        //world to camera (triangles)
        //

        // create Point_At matrix
        // TODO ideally, these should all be 3 dimensional vectors (use 3x3 rot matrix)
        models_out = models_culled;
        Vec<double> up = {0.0, 1.0, 0.0};
        Vec<double> cam_pos = (world_to_camera_translate.get_row(3)).slice(0, 3);
        Vec<double> look_dir = {0.0, 0.0, 1.0, 1.0}; 
        look_dir = look_dir * rotate_x(cam_pitch) * rotate_y(cam_yaw);
        Vec<double> target = cam_pos + look_dir.slice(0, 3);
        look_at = look_at_matrix(cam_pos, target, up);
       
        // old view code
        //models_out = models_culled;
        //auto A = world_to_camera_rotate.get_row(0);
        //A = A.slice(0, 3);
        //auto B = world_to_camera_rotate.get_row(1);
        //B = B.slice(0, 3);
        //auto C = world_to_camera_rotate.get_row(2);
        //C = C.slice(0, 3);
        //auto T = world_to_camera_translate.get_row(3);
        //T = T.slice(0, 3);
        //look_at = world_to_camera_rotate.transpose();
        ////look_at = world_to_camera_rotate;
        //look_at.set(3, 0, - A.dot(T));
        //look_at.set(3, 1, - B.dot(T));
        //look_at.set(3, 2, - C.dot(T));
        //look_at.set(3, 3, 1.0);
         
        static bool asdf1 = false;
        if (asdf1 == false) {
            asdf1 = true;
            cout << "printing look_at matrix\n" << look_at << '\n';
            cout << "printing world_to_camera_rotate\n" << world_to_camera_rotate << '\n';
        }
        
        for (auto& model : models_out) {
            for (int i = 0; i < model.tris.size(); ++i) {
                model.tris[i].verts[0] = model.tris[i].verts[0] * look_at;
                model.tris[i].verts[1] = model.tris[i].verts[1] * look_at;
                model.tris[i].verts[2] = model.tris[i].verts[2] * look_at;
            }
        }
        models_view = models_out;
        //
        //world to camera (lines)
        //
        model_index = 0;
        for (auto& output : lines_outputs) {
            for (int i = 0; i < output.size(); ++i) {
                output[i] = lines_world[model_index][i] * look_at;
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
        models_out = models_view;
        model_index = 0;
        for (auto& model : models_out) {
            for (auto& tri : model.tris) {
                for (auto& vert : tri.verts) {
                    vert = vert * clip_matrix * ortho_matrix;
                }
            }
        }
        models_clipped = models_out;
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

        // 
        // clip out of bounds vertices (triangles)
        //
        //model_index = 0;
        ////model_outputs_clipped = model_outputs;
        //// input is: models_clipped, output is:  
        auto models_clipped_out = models_clipped;
        int red_before = models_clipped[0].color.r;
        int red_after = models_clipped_out[0].color.r;

        static bool asdf1234 = false;
        if (!asdf1234) {
            asdf1234 = true;
            cout << "red before = " << red_before << '\n';
            cout << "red after = " << red_after << '\n';
        }


        for (auto& m : models_clipped_out) {
            m.tris.clear();
        }
        //clip near z plane
        int i_model = 0;
        vector<Vec<double>> clipped_points;
        clipped_points.reserve(6);
        for (const auto& model : models_clipped) {
            //models_clipped_out[i_model].tris.reserve(model.tris.size() * 3 * 2);
            models_clipped_out[i_model].tris.reserve(model.tris.size() * 3 * 2);
            //assert(models_clipped_out[i_model].tris.size())
            assert(model.tris[0].verts[0].length > 0 && model.tris[0].verts[0].length <= 4); 
            //if ( (model.tris.size() > 0 && model.tris.size() < 1000000) == false) {
            //    cout << "crashing! " << i_model << ", " << model.tris.size() << '\n';
            //}
            //assert(model.tris.size() > 0 && model.tris.size() < 1000000); //TODO crash here occasionally
            for (const auto& tri : model.tris) {
                

                //// check for simplest case where all vec's are inside or all outside clipping plane
                // TODO perspective and ortho mactrices have inaccuries.
                // left, right, bottom and top clipping happens inside viewing
                // near clipping works fine
                //// bottom clip course
                //if ( tri.verts[0].data[1] < -1.0 * tri.verts[0].data[3]
                //        && tri.verts[1].data[1] < -1.0 * tri.verts[1].data[3]
                //        && tri.verts[2].data[1] < -1.0 * tri.verts[2].data[3]) {
                //    
                //    continue;
                //}
                //// left clip course 
                //if ( tri.verts[0].data[0] < -1.0 * tri.verts[0].data[3]
                //        && tri.verts[1].data[0] < -1.0 * tri.verts[1].data[3]
                //        && tri.verts[2].data[0] < -1.0 * tri.verts[2].data[3]) {
                //    
                //    continue;
                //}
                //// right clip course 
                //if ( tri.verts[0].data[0] > 1.0 * tri.verts[0].data[3]
                //        && tri.verts[1].data[0] > 1.0 * tri.verts[1].data[3]
                //        && tri.verts[2].data[0] > 1.0 * tri.verts[2].data[3]) {
                //    
                //    continue;
                //}





                // check for simplest case where all vec's are inside or all outside clipping plane
                // check if all outside
                if ( tri.verts[0].data[2] > tri.verts[0].data[3]
                        && tri.verts[1].data[2] > tri.verts[1].data[3]
                        && tri.verts[2].data[2] > tri.verts[2].data[3]) {
                    models_clipped_out[i_model].tris.push_back( tri );
                    continue;  

                } else if ( tri.verts[0].data[2] < tri.verts[0].data[3]
                        && tri.verts[1].data[2] < tri.verts[1].data[3]
                        && tri.verts[2].data[2] < tri.verts[2].data[3]) {
                    
                    continue;
                }
                Vec<double> vert_prev = tri.verts[2]; // last vertex is [3], current is [0]
                bool prev_inside = vert_prev.data[2] >  1.0 * vert_prev.data[3]; // check if inside near plane 
                //bool prev_inside = vert_prev.data[2] <  vert_prev.data[3]; // check if inside near plane 
                // [2] is Z value!!!!  [3] is W value!!!!
                int verts_added = 0;
                clipped_points.clear(); 
                for (auto& vert_cur : tri.verts) {
                    bool cur_inside = vert_cur.data[2] >  1.0 * vert_cur.data[3];
                    //bool cur_inside = vert_cur.data[2] < vert_cur.data[3];

                    if ( (!prev_inside) != (!cur_inside) ) {
                        double lerp_factor = (vert_prev.data[3] - vert_prev.data[2]) 
                                / ( ( vert_prev.data[3] - vert_prev.data[2]) 
                                - ( vert_cur.data[3] - vert_cur.data[2]) );
                        clipped_points.push_back( vert_prev.lerp(vert_cur, lerp_factor));
                        verts_added++; // should just use size of vector?
                    }

                    if (cur_inside) {
                        clipped_points.push_back(vert_cur);
                        verts_added++;
                    }

                    prev_inside = cur_inside;
                    vert_prev = vert_cur;
                }
                // 0, 3, 4 vertices should be generated, so, 0, 1, 2 triangles should be generated
                //cout << "verts added = " << verts_added << '\n';
                assert(verts_added == 0 || verts_added == 3 || verts_added == 4);
                if (verts_added == 0) {
                    //do nothing
                } else if ( verts_added == 3) {
                    models_clipped_out[i_model].tris.push_back( 
                            Tri<double>(clipped_points[0], 
                            clipped_points[1], 
                            clipped_points[2], 
                            tri.color, 
                            tri.draw_wire_frame,
                            tri.draw_face));
                            //models_clipped_out[i_model].color));
                } else if ( verts_added == 4) {
                    //cout << "clipping happened!\n";
                    models_clipped_out[i_model].tris.push_back( 
                            Tri<double>(clipped_points[0], 
                            clipped_points[1], 
                            clipped_points[2], 
                            tri.color,
                            tri.draw_wire_frame,
                            tri.draw_face));
                    models_clipped_out[i_model].tris.push_back( 
                            Tri<double>(clipped_points[0], 
                            clipped_points[2], 
                            clipped_points[3], 
                            tri.color,
                            tri.draw_wire_frame,
                            tri.draw_face));
                }
            }
            models_clipped_out[i_model].color = model.color;

            i_model++;
        }


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
        models_out = models_clipped_out;
        if (!posted_info) { cout << "homogenous divide vertices" << endl; }
        model_index = 0;
        for (auto& model : models_out) {
            for (auto& tri : model.tris) {
                for (auto& vert : tri.verts) {
                    double w = vert.data[3];
                    vert = (1.0/w) * vert; 
                }
            }
        }
        models_ndc = models_out;
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
        }
        //
        // project to screen (triangles)
        //
        if (!posted_info) { cout << "screen vertices" << endl; }
        models_out = models_ndc;
        for (auto& model : models_out) {
            for (auto& tri : model.tris) {
                for (auto& vert : tri.verts) {
                    vert.data[0] = 1.0 * vert.data[0] * win_res_x 
                            /(vert.data[2] * 2.0 ) + 1000.0/2.0;
                    vert.data[1] = -1.0 * vert.data[1] * win_res_y 
                            /(vert.data[2] * 2.0 ) + 1000.0/2.0;
                }
            }
        }
        models_screen = models_out;
        //
        // project to screen (lines)
        //
        if (!posted_info) { cout << "screen vertices" << endl; }
        for (auto& output : lines_outputs_clipped) {
            for (int i = 0; i < output.size(); ++i) {
                output[i].data[0] = 1.0 * output[i].data[0] * win_res_x
                    / (output[i].data[2] * 2.0) + 1000.0/2.0;
                output[i].data[1] = -1.0 * output[i].data[1] * win_res_y
                    / (output[i].data[2] * 2.0) + 1000.0/2.0;
            }
        }

        // draw debug stuff before rastering

        //
        // raster lines
        //
        // const int origin_index = 2;
        // bool dont_render_origin = false;
        if (RENDER_LINES) {
        int i = 0;
        SDL_SetRenderDrawColor(gsdl.renderer, 0, 50, 0, 20);
        for (auto& output : lines_outputs_clipped) {
            if(i < lines_outputs_clipped_count) {
                i++;
            } else {
                break;
            }
            for (int i_line = 0; i_line < output.size(); i_line += 2) {
                SDL_RenderDrawLine(gsdl.renderer, 
                        output[i_line].data[0], output[i_line].data[1],
                        output[i_line+1].data[0], output[i_line+1].data[1]);
            }
            //dont_draw_line:;
            model_index++;
        } 
        }

        //
        // raster triangles ( as wireframe or filled in)
        //
        // but first sort them (fill up Tri's buffer)
        //vector<Tri<double>> tris;
        //for (int j = 0; j < model_outputs_clipped.size(); ++j) {
        //    Color c;
        //    bool draw_wire_frame = true;
        //    
        //    for (int i = 0; i < model_outputs_clipped[j].size(); i += 3) {
        //        if (j != 2) {
        //            c = {50, 50, 50, 100};
        //        } else {
        //            draw_wire_frame = false;
        //            if (i == 0) {
        //                c = {100, 0, 0, 100};
        //            } else if (i == 1) {
        //                c = {0, 100, 0, 100};
        //            } else {
        //                c = {0, 0, 100, 100};
        //            }
        //        }
        //        Tri<double> tri{model_outputs_clipped[j], i, c, draw_wire_frame};
        //        tris.push_back(tri);
        //    }
        //}
        //vector1.insert( vector1.end(), vector2.begin(), vector2.end() );
        //vector<Tri<double>> tris = models_screen[1].tris;
        //tris.insert(tris.end(), models_screen[1].tris.begin(), models_screen[1].tris.end());
        //tris.insert(tris.end(), models_screen[2].tris.begin(), models_screen[2].tris.end());
       


        vector<Tri<double>> tris = models_screen[0].tris;
        //tris.insert(tris.end(), models_screen[1].tris.begin(), models_screen[1].tris.end());
        //tris.insert(tris.end(), models_screen[2].tris.begin(), models_screen[2].tris.end());


        for (int i = 1; i < models_screen.size(); ++i) {
            tris.insert( tris.end(), models_screen[i].tris.begin(), models_screen[i].tris.end());
        }

        sort(tris.begin(), tris.end(), [](Tri<double>& tri0, Tri<double>& tri1) {
            double z0 = (tri0.verts[0].data[2] + tri0.verts[1].data[2] 
                + tri0.verts[2].data[2]) / 3.0;
            double z1 = (tri1.verts[0].data[2] + tri1.verts[1].data[2] 
                + tri1.verts[2].data[2]) / 3.0;
            return z0 > z1;
        });
        for (auto& tri : tris) {
            if (tri.draw_face) {
                draw_triangle(tri); 
            }
            // ignore for now
            //if (!tri.draw_wire_frame) {
            //    continue;
            //}
            // now draw wireframe
            // i_tri is index to begining of triangle
            if (tri.draw_wire_frame) { 
                for (int i = 0; i < 3; ++i) {
                    // draws wireframe
                    // i at vertex of triangle
                    int j = (i + 1 < 3) ? i + 1 : 0;
                    SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
                     
                    SDL_RenderDrawLine(gsdl.renderer, 
                            tri.verts[i].data[0], tri.verts[i].data[1],
                            tri.verts[j].data[0], tri.verts[j].data[1]);
                }
            }
        }

        
        // as TTF_RenderText_Solid could only be used on
        // SDL_Surface then you have to create the surface first
        std::ostringstream sstr{};
        sstr << "ship.pos : " << model_0_pos;
        draw_text(sstr.str(), 0, 0);
        sstr.str("");
        sstr.clear();
        sstr << "teapot.pos : " << model_1_pos;
        draw_text(sstr.str(), 0, 20);
        sstr.str("");
        sstr.clear();
        sstr << "cam : " << look_at.get_row(0) ;
        draw_text(sstr.str(), 0, 40);
        sstr.str("");
        sstr.clear();
        sstr << "cam : " << look_at.get_row(1) ;
        draw_text(sstr.str(), 0, 60);
        sstr.str("");
        sstr.clear();
        sstr << "cam : " << look_at.get_row(2) ;
        draw_text(sstr.str(), 0, 80);
        sstr.str("");
        sstr.clear();
        sstr << "cam : " << look_at.get_row(3) ;
        draw_text(sstr.str(), 0, 100);
        sstr.str("");
        sstr.clear();
        sstr << "cam world: " << world_to_camera_translate.get_row(3) ;
        draw_text(sstr.str(), 0, 120);
        //sstr.str("");
        //sstr.clear();
        //sstr << "BF cull: " << BACK_FACE_CULLING;
        //draw_text(sstr.str(), 0, 140);
        sstr.str("");
        sstr.clear();
        sstr << "look_dir: " << look_dir;
        draw_text(sstr.str(), 0, 140);
        sstr.str("");
        sstr.clear();
        sstr << "box0 tri[0].color" << models_world[0].tris[0].color;
        draw_text(sstr.str(), 0, 160);

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

