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
#include <thread>
#include <utility>
#include <stdexcept>

#include "vec.hpp"
//#include "test.hpp"
#include "sdl_init.hpp"

//ctags -R --languages=C++ --extra=+fq --exclude=.git .
using std::cout, std::endl, std::vector, std::array, std::string, std::swap,
      std::thread, std::pair;

// vertexes of object
#define NUM_VERTICES 9
#define NUM_MESHES 3
#define ROTATE_FIRST 0
#define RENDER_LINES 0  // renders grid lines
#define NUM_LINES 22
#define NUM_THREADS 8
#define NUM_ATTRB_PER_VERT 12
bool BACK_FACE_CULLING = true;
bool Custom_Draw_Line = true;
bool Custom_Draw_Triangle = true;

Vec<double> model_0_pos;
Vec<double> model_1_pos;
vector<Model<double>> models;

// vertices drawn to screen
array<vector<Vec<double>>, NUM_MESHES> model_outputs = // TODO DELETE
{
    vector<Vec<double>>{},
    vector<Vec<double>>{},
    vector<Vec<double>>{}
};


// TODO delete these?
auto models_world = models; // world space of models 
auto models_culled = models; // world space of models with back face culling applied
auto models_view = models; // view space of models
auto models_clipped = models; // clip space of models
auto models_ndc = models;  // normalized device coordinates of models
auto models_screen = models;  // screen space
auto models_out = models;  // perhaps could use a single variable instead of all the above

//
// new verts variables
//

// inputs 
vector<Vec4> pos_in;  // concat'd list of all input positions
vector<Vec4> norms_in;
vector<Vec4> norm_faces_in;
vector<Vec4> colors_in; // per vertex
vector<Vec3> tex_in;  // 3rd element is 'w'
vector<unsigned int> indices;
vector<Material*> materials_in;
vector<Mat4x4*> model_world_transforms;
vector<Mat4x4*> model_world_rotations; // per vertex rotations
vector<Mat4x4*> model_world_rotations_faces; // per face/tri rotation
Vec3 cam_pos;
double cam_yaw = 0.0;
double cam_pitch = 0.0;

vector<Vec4> pos_world;   // per vert
vector<Vec4> norm_world;  // per vert
vector<Vec4> norm_faces_world;  // per tri
vector<Vec4> color_faces_world; // per tri
vector<bool> tris_culled;

vector<Vec4> pos_pre_clip;
vector<Vec4> pos_clip; 
vector<unsigned int> indices_clip;
vector<bool> tris_clipped;

vector<Vec4> pos_ndc;
vector<Vec4> pos_screen;
//
// end new verts variables
//

// basic lighting
Vec<double, 3> light_dir = {-1.0, -1.0, 1.0};
LightDirectional<double> light_directional { 
        Vec3{-1.0, 0.0, 1.0}, Vec4{1.0, 1.0, 1.0, 1.0}, 1.0};
//Vec4 light_ambient = {0.1, 0.1, 0.1, 1.0};
LightAmbient<double> light_ambient { Vec4{1.0, 1.0, 1.0, 1.0}, 0.1};

int lines_outputs_clipped_count = 0;
Mesh<double>* mesh_selected = nullptr; // model selected for transformations

const int win_res_x = 1000;
const int win_res_y = 1000;
float depth_buffer[win_res_y][win_res_x];
// meshes
//Control_Select control_select = Control_Select::model_0;
enum Control_Select : unsigned int { m_0 = 0, m_1 };
Control_Select control_select = m_0; // TODO delete
bool print_clip_vertices = false; // TODO delete?
bool print_info = false;
std::array<std::array<Vec<double>, NUM_VERTICES>, NUM_MESHES> meshes; 



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
//Mat<double> world_to_camera = {  
//    Vec{1.0, 0.0, 0.0, 0.0}, 
//    Vec{0.0, 1.0, 0.0, 0.0},
//    Vec{0.0, 0.0, 1.0, 0.0},
//    Vec{0.0, 0.0, 0.0, 1.0}
//};

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

// all 4x4 matrices
Mat<double> ortho_matrix;
Mat<double> clip_matrix; // not used yet, to be used for cliping vertices outside view frustrum.
Mat<double> world_to_camera_matrix;
Mat<double> world_to_clip_matrix;
Mat<double> turn_left;
Mat<double> turn_right;
Mat<double> turn_down;
Mat<double> turn_up;


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

void draw_line(int x0, int y0, int x1, int y1)
{
    int dx = abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int error = dx + dy;
    int e2;

    while (true) {
        //plot(x0, y0);
        SDL_RenderDrawPoint(gsdl.renderer, x0, y0);

        e2 = 2 * error;
        if (e2 >= dy) {
            if (x0 == x1) { break; }
            error = error + dy;
            x0 = x0 + sx;
        }
        if (e2 <= dx) {
            if (y0 == y1) { break; }
            error = error + dx;
            y0 = y0 + sy;
        }
    }
}

void draw_triangle(	    int x1, int y1, float u1, float v1, float w1, Vec4 n1,
						int x2, int y2, float u2, float v2, float w2, Vec4 n2,
						int x3, int y3, float u3, float v3, float w3, Vec4 n3,
                        //float light_intensity,
                        LightDirectional<double>& _light_directional,
                        LightAmbient<double>& _light_ambient,
	                    SDL_Surface* _surface = nullptr)
{
    // sort y's in ascending order
	if (y2 < y1) {
		swap(y1, y2);
		swap(x1, x2);
		swap(u1, u2);
		swap(v1, v2);
		swap(w1, w2);
        swap(n1, n2);
	}
	if (y3 < y1) {
		swap(y1, y3);
		swap(x1, x3);
		swap(u1, u3);
		swap(v1, v3);
		swap(w1, w3);
        swap(n1, n3);
	}
	if (y3 < y2) {
		swap(y2, y3);
		swap(x2, x3);
		swap(u2, u3);
		swap(v2, v3);
		swap(w2, w3);
        swap(n2, n3);
	}

	int dy1 = y2 - y1;
	int dx1 = x2 - x1;
	float dv1 = v2 - v1;
	float du1 = u2 - u1;
	float dw1 = w2 - w1;
    Vec4 dn1 = n2 - n1; // vec
     
	int dy2 = y3 - y1;
	int dx2 = x3 - x1;
	float dv2 = v3 - v1;
	float du2 = u3 - u1;
	float dw2 = w3 - w1;
    Vec4 dn2 = n3 - n1;

	float tex_u, tex_v, tex_w;
    Vec4 norm;

	float dax_step = 0, dbx_step = 0,
		dw1_step = 0, dw2_step = 0,
		du1_step = 0, dv1_step = 0,
		du2_step = 0, dv2_step = 0;

    Vec4 dn1_step, dn2_step;

	if (dy1) dax_step = dx1 / (float)abs(dy1);
	if (dy2) dbx_step = dx2 / (float)abs(dy2);

	if (dy1) du1_step = du1 / (float)abs(dy1);
	if (dy1) dv1_step = dv1 / (float)abs(dy1);
	if (dy1) dw1_step = dw1 / (float)abs(dy1);

	if (dy2) du2_step = du2 / (float)abs(dy2);
	if (dy2) dv2_step = dv2 / (float)abs(dy2);
	if (dy2) dw2_step = dw2 / (float)abs(dy2);

    if (dy1) dn1_step = dn1 / (float)abs(dy1); 
    if (dy2) dn2_step = dn2 / (float)abs(dy2); 

	if (dy1)
	{
		for (int i = y1; i <= y2; i++)
		{
			int ax = x1 + (float)(i - y1) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;

			float tex_su = u1 + (float)(i - y1) * du1_step;
			float tex_sv = v1 + (float)(i - y1) * dv1_step;
			float tex_sw = w1 + (float)(i - y1) * dw1_step;

			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = w1 + (float)(i - y1) * dw2_step;

            Vec4 norm_s = n1 + (double)(i - y1) * dn1_step;
            Vec4 norm_e = n1 + (double)(i - y1) * dn2_step;

			if (ax > bx) // ax is to be smaller than bx
			{
				swap(ax, bx);
				swap(tex_sw, tex_ew);
				swap(tex_su, tex_eu);
				swap(tex_sv, tex_ev);
                swap(norm_s, norm_e);
			}

			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;
            norm = norm_s;

			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;

			for (int j = ax; j < bx; j++)
			{
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
                norm  = (1.0 - t) * norm_s + (double)t * norm_e;
				//if (tex_w > depth_buffer[i*ScreenWidth() + j])
                // proper clipping would eliminate this check
                if (i >= 0 && i < win_res_y && j >= 0 && j < win_res_x) {
				    if (tex_w < depth_buffer[i][j]) // originally <
				    {
                        Uint32 pixel = *((Uint32*)_surface->pixels 
                            + (unsigned int)((tex_v * _surface->h)) * _surface->w 
                            + (unsigned int)((tex_u * _surface->w)) );
                        //Uint8 r, g, b, a;
                        Vec<Uint8, 4> c;//Color<Uint8> c;  // initially stores diffuse, then stores output
                        SDL_GetRGBA(pixel, _surface->format, &c.r, &c.g, &c.b, &c.a);
                        
                        Vec3 norm3d = slice<double, 3, 4>(norm);
                        norm3d.normalize();
                        double dp = std::max(0.0, norm3d.dot( -1.0 * _light_directional.dir));

                        c.r = (Uint8)(dp * (double)c.r * _light_directional.color.r);
                        c.g = (Uint8)(dp * (double)c.g * _light_directional.color.g);
                        c.b = (Uint8)(dp * (double)c.b * _light_directional.color.b);

                        SDL_SetRenderDrawColor(gsdl.renderer, c.r, c.g, c.b, 255);  //20);

                        SDL_RenderDrawPoint(gsdl.renderer, j, i);
                        depth_buffer[i][j] = tex_w; 
                        
				    }
                }
				t += tstep;
			}
		}
	}

	dy1 = y3 - y2;
	dx1 = x3 - x2;
	dv1 = v3 - v2;
	du1 = u3 - u2;
	dw1 = w3 - w2;
    dn1 = n3 - n2; 

	if (dy1) dax_step = dx1 / (float)abs(dy1);
	if (dy2) dbx_step = dx2 / (float)abs(dy2); // redundant?

	du1_step = 0, dv1_step = 0;
	if (dy1) du1_step = du1 / (float)abs(dy1);
	if (dy1) dv1_step = dv1 / (float)abs(dy1);
	if (dy1) dw1_step = dw1 / (float)abs(dy1);

    if (dy1) dn1_step = dn1 / (float)abs(dy1);

	if (dy1)
	{
		for (int i = y2; i <= y3; i++)
		{
			int ax = x2 + (float)(i - y2) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;

			float tex_su = u2 + (float)(i - y2) * du1_step;
			float tex_sv = v2 + (float)(i - y2) * dv1_step;
			float tex_sw = w2 + (float)(i - y2) * dw1_step;

			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = w1 + (float)(i - y1) * dw2_step;

            Vec4 norm_s = n2 + (double)(i - y2) * dn1_step; // (n3 - n2) / ( y3 - y2)
            Vec4 norm_e = n1 + (double)(i - y1) * dn2_step; 

			if (ax > bx)
			{
				swap(ax, bx);
				swap(tex_sw, tex_ew);
				swap(tex_su, tex_eu);
				swap(tex_sv, tex_ev);
                swap(norm_s, norm_e);
			}

			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;
            norm = norm_s;

			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;

			for (int j = ax; j < bx; j++)
			{
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
                norm = (1.0 - t) * norm_s + (double)t * norm_e;

                if (i >=0 && i < win_res_y && j >= 0 && j < win_res_x) {
				    if (tex_w < depth_buffer[i][j]) // originally <
				    {
                        Uint32 pixel = *((Uint32*)_surface->pixels 
                            + (unsigned int)((tex_v * _surface->h)) * _surface->w 
                            + (unsigned int)(((tex_u * _surface->w))) );

                        //Uint8 r, g, b, a;
                        //SDL_GetRGBA(pixel, _surface->format, &r, &g, &b, &a);
                        //Uint8 r, g, b, a;
                        Vec<Uint8, 4> c;
                        SDL_GetRGBA(pixel, _surface->format, &c.r, &c.g, &c.b, &c.a);

                        Vec3 norm3d = slice<double, 3, 4>(norm);
                        norm3d.normalize();
                        double dp = norm3d.dot( -1.0 * _light_directional.dir);
                        dp = std::max(0.0, dp);
                        dp = 0.0 + dp * 1.0; // ambient term + diffuse term?
                        //dp = 0.0;
                        c.r = (Uint8)( dp * (double)c.r * _light_directional.color.r);
                        c.g = (Uint8)(dp * (double)c.g * _light_directional.color.g);
                        c.b = (Uint8)(dp * (double)c.b * _light_directional.color.b);

                        SDL_SetRenderDrawColor(gsdl.renderer, c.r, c.g, c.b, 255);

                        SDL_RenderDrawPoint(gsdl.renderer, j, i);
                        depth_buffer[i][j] = tex_w; // pDepthBuffer[i*ScreenWidth() + j] = tex_w;
				    }
                }
				t += tstep;
			}
		}
	}
}

void model_world_transform(
        const vector<Vec<double, 4>>& _verts, 
        vector<Vec<double, 4>>& _verts_out,
        //Mat<double, 4, 4>& transform, 
        const int _start, const int _end) 
{
    for (int i = _start; i < _start + _end; ++i) {
        ////verts[i] = verts[i] * transform; 
        //_verts_out[i] = _verts[i] * models[verts_i[i]].srt;
    }
}

template<typename T>
Mat<T, 4, 4> look_at_matrix(Vec<T, 3>& _pos, Vec<T, 3>& _target, Vec<T, 3>& _up) 
{
    //assert(_pos.length == 3 && _target.length == 3 && _up.length == 3);
    // new forward
    Vec<T, 3> forward = _target - _pos;
    forward.normalize(); 
    // new up
    Vec<T, 3> a = forward.dot(_up) * forward; 
    Vec<T, 3> up = _up - a;
    up.normalize(); 
    // new right
    Vec<T, 3> right = up.cross(forward);
    right.normalize(); 
    // construct 3x3 R matrix
    Mat<T, 3, 3> R;
    R.set_row(0, right);
    R.set_row(1, up);
    R.set_row(2, 1.0 * forward);
    R = R.transpose(); // R is now inverted
    
    Vec<T, 3> t = ( -1.0 * _pos) * R;  // translation portion

    // construct matrix
    Mat<T> m; // 4x4
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
    Vec<double, 2> a = { 1.0, 1.0 };
    Vec<double, 2> b = { 5.0, 3.0 };
    //( 1 1 0 0 )
    //( 3 2 0 0 )
    //( 5 3 0 0 )


    cout << "testing lerp\n";
    cout << a.lerp(b, 0.0) << endl << a.lerp(b, 0.5) << endl << a.lerp(b, 1.0) << endl;
}

void split_range( const int _range, 
                  const int _splits, 
                  vector<pair<int, int>> _pairs) 
{
    _pairs.clear();
    //int models_size = 34;
    //vector<pair<int, int>> i_ranges;
    int v_per_t = _range / _splits;
    int v_per_t_remain = _range % _splits;
    for (int i = 0; i < _splits; ++i) {
        _pairs.emplace_back(i * v_per_t, (i + 1) * v_per_t);
    }
    _pairs[_splits - 1].second += v_per_t_remain;
   
    for (int i = 0; i < _splits; ++i) {
        cout << _pairs[i].first << ", " << _pairs[i].second << '\n';
    }

}

void clear_vertices( )
{
    // input buffers
    pos_in.clear();
    norms_in.clear();
    norm_faces_in.clear();
    colors_in.clear();
    tex_in.clear();
    materials_in.clear();
    model_world_transforms.clear();
    model_world_rotations.clear();
    model_world_rotations_faces.clear();
    indices.clear();
    // world 
    pos_world.clear(); 
    norm_world.clear(); 
    norm_faces_world.clear();
    color_faces_world.clear();

    pos_pre_clip.clear();
    pos_clip.clear();
    indices_clip.clear();
    tris_clipped.clear();

    pos_ndc.clear();
    pos_screen.clear();
}

// done usually once per frame (if scene is very dynamic)
void init_vertices( )
{
    clear_vertices(); 
    //int max_index = 0;
    for (Model<double>& model : models) {
        assert(model.meshes.size() == 1);
        for (Mesh<double>& mesh : model.meshes) {
            //int v_count = mesh.pos.size();
            int max_index = pos_in.size();
            for (Vec4 p : mesh.pos) { 
                pos_in.emplace_back(p); 
                model_world_transforms.push_back(&mesh.srt);
            }
            for (Vec4 n : mesh.norms) { 
                norms_in.emplace_back(n); 
                model_world_rotations.push_back(&mesh.rot);
            }
            for (Vec4 nf : mesh.norm_faces) { 
                norm_faces_in.emplace_back(nf); 
                model_world_rotations_faces.push_back(&mesh.rot);
                materials_in.push_back(&mesh.mat);
            }
            for (Vec4 c : mesh.colors) { colors_in.emplace_back(c); }
            for (Vec3 t : mesh.texs) { tex_in.emplace_back(t); }
            //int n_verts = indices.size() / NUM_ATTRB_PER_VERT; // 12
            assert( indices.size() % NUM_ATTRB_PER_VERT == 0 );
            // pos, norm, tex, color
            mesh.srt = mesh.scale * mesh.rot * mesh.trans;
            
            // find max index
            for ( 
                unsigned int i = 0; 
                i < mesh.indices.size();
                i += NUM_ATTRB_PER_VERT )
            {
                // pos vertex indices
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 0]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 1]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 2]);
                // norm vertex indices 
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 3]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 4]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 5]);
                // tex vertex indices 
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 6]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 7]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 8]);
                // color vertex indices
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 9]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 10]);
                indices.push_back( max_index /** 3*/ + mesh.indices[i + 11]);
            }
        }
    }
    // initialize output buffers
    pos_world = pos_in; // inefficient!
    //pos_world.clear();
    norm_world = norms_in; //inefficient
    //norm_world.clear();
    //norm_faces_world[i_tri] = norm_faces_in[i_tri] * 
    //        *model_world_rotations[i_tri];
    int num_faces = indices.size() / NUM_ATTRB_PER_VERT;
    assert( indices.size() % NUM_ATTRB_PER_VERT == 0 );
    norm_faces_world.resize( num_faces );
    tris_culled.resize( num_faces );
    color_faces_world.resize( num_faces);
    pos_pre_clip.resize( pos_world.size() );
    pos_clip.resize( pos_world.size() ); //newer
    indices_clip.reserve( indices.size() * 2 ); // RESERVE!?
    tris_clipped.resize( num_faces );
    pos_ndc.resize( pos_world.size() );
    pos_screen.resize( pos_world.size() );
    
    
}

void model_to_world(int _start_v, int _end_v) 
{
    // uniform inputs: pos_in, norms_in, model_world_transforms, color_faces
    // maybe don't need _b_index
    for (int i = _start_v; i < _end_v; ++i) {
        pos_world[i] = pos_in[i] * (*(model_world_transforms[i]));
        // sizes: 1, 1, 2
        //norm_world[i] = norms_in[i] * (*(model_world_rotations[i]));
    }
    

//    for (int i = 0; i < norms_in.size(); i++) {
//        norm_world[i] = norms_in[i] * (*(model_world_rotations[i]));
//    }

}

void model_to_world_normals(int _start_n, int _end_n)
{
    for (int i = _start_n; i < _end_n; i++) {
        norm_world[i] = norms_in[i] * (*(model_world_rotations[i]));
    }
}

void model_to_world_part2(int _start_tri, int _end_tri) 
{
    //vector<Vec4> norm_faces_in;
    //vector<Vec4> norm_faces_buffers;
    for (int i_tri = _start_tri; i_tri < _end_tri; ++i_tri) {
        // compute world normal
        norm_faces_world[i_tri] = norm_faces_in[i_tri] * 
                *  model_world_rotations_faces[i_tri];
        norm_faces_world[i_tri].normalize();
        ////// back face culling, set tris_culled[i_tri]
        ////// cam_pos
        ////Vec4 p_tri = pos_world[indices[i_tri * 12] + 0]; // point on triangle
        ////Vec4 p_cam = world_to_camera_translate.get_row(3);
        ////Vec4 v = p_tri - p_cam;
        ////v.w = 0;
        ////v.normalize();
        ////double dp0 = v.dot(norm_faces_world[i_tri]);
        ////tris_culled[i_tri] = dp0 < 0.0;
        // back face culling, set tris_culled[i_tri]
        // cam_pos
        Vec3 p_tri = slice<double, 3, 4>(pos_world[indices[i_tri * 12] + 0]); // point on triangle
        Vec3 v = p_tri - cam_pos;
        v.normalize();
        double dp0 = v.dot(slice<double, 3, 4>(norm_faces_world[i_tri]));
        tris_culled[i_tri] = dp0 >= 0.3; // TODO something wrong  //0.100;
        
        // compute lighting
        //Vec<double, 3> c = -1.0 * a0.cross4d(b0);
        //c.normalize(); //c = c.unit();
        Vec3 c = slice<double, 3, 4>(norm_faces_world[i_tri]);
        c = -1.0 * c;
        //double dp = c.dot(slice<double, 3, 4>(light_directional.dir));
        double dp = c.dot(light_directional.dir);
        dp = std::max(0.0, dp);
        (color_faces_world[i_tri]).r = dp * light_directional.color.r;
        (color_faces_world[i_tri]).g = dp * light_directional.color.g;
        (color_faces_world[i_tri]).b = dp * light_directional.color.b;
        //(color_faces_world[i_tri]).a = light_directional.b;

        //dp = dp * 255;
        //int dp_int = int(dp); 
        //dp_int = std::max( 0, dp_int);
        //models_world[j].tris[i].color = models[j].color * dp_int ; // set initial color

         
    }
}

void world_to_clip(int _start_v, int _end_v) {
    // _start and _end represent vertices, not faces
    for (int i = _start_v; i < _end_v; ++i) {
        pos_clip[i] = pos_world[i] * world_to_clip_matrix;
        //pos_pre_clip[i] = pos_world[i] * world_to_clip_matrix;
        ////pos_pre_clip[i] = pos_world[i] * look_at_matrix * clip_matrix;
    }
}

//vector<Vec4> pos_pre_clip;
//vector<Vec4> pos_clip; 
//vector<unsigned int> indices_clip;


void clip(int _start_tri, int _end_tri) 
{
    indices_clip = indices; 
    for (int i = _start_tri * NUM_ATTRB_PER_VERT;
            i < _end_tri * NUM_ATTRB_PER_VERT;
            i += NUM_ATTRB_PER_VERT) 
    {

        //if (tris_culled[i/NUM_ATTRB_PER_VERT]) {
        //    continue;
        //}

        //// check for simplest case where all vec's are inside or all outside clipping plane
        // bottom clip course
        if (    pos_clip[indices[i + 0]].y < - pos_clip[indices[i + 0]].w
             && pos_clip[indices[i + 1]].y < - pos_clip[indices[i + 1]].w
             && pos_clip[indices[i + 2]].y < - pos_clip[indices[i + 2]].w )
        {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // CLIPPED!
            continue;
        }
        // top clip course 
        if (    pos_clip[indices[i + 0]].y > pos_clip[indices[i + 0]].w
             && pos_clip[indices[i + 1]].y > pos_clip[indices[i + 1]].w
             && pos_clip[indices[i + 2]].y > pos_clip[indices[i + 2]].w )
        {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // CLIPPED!
            continue;
        }
        // left clip course 
        if (    pos_clip[indices[i + 0]].x < - pos_clip[indices[i + 0]].w
             && pos_clip[indices[i + 1]].x < - pos_clip[indices[i + 1]].w
             && pos_clip[indices[i + 2]].x < - pos_clip[indices[i + 2]].w )
        {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // CLIPPED!
            continue;
        }
        // right clip course 
        if (    pos_clip[indices[i + 0]].x > pos_clip[indices[i + 0]].w
             && pos_clip[indices[i + 1]].x > pos_clip[indices[i + 1]].w
             && pos_clip[indices[i + 2]].x > pos_clip[indices[i + 2]].w )
        {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // CLIPPED!
            continue;
        }
        // near clip course 
        if ( pos_clip[indices[i + 0]].z >= - pos_clip[indices[i + 0]].w &&
             pos_clip[indices[i + 1]].z >= - pos_clip[indices[i + 1]].w &&
             pos_clip[indices[i + 2]].z >= - pos_clip[indices[i + 2]].w ) 
        {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 0; // NOT CLIPPED!
            continue;
        } 
        else if ( pos_clip[indices[i + 0]].z < - pos_clip[indices[i + 0]].w &&
                  pos_clip[indices[i + 1]].z < - pos_clip[indices[i + 1]].w &&
                  pos_clip[indices[i + 2]].z < - pos_clip[indices[i + 2]].w ) 
        {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // CLIPPED!
            continue;
        } 

        Vec4 pos_prev = pos_clip[indices[i + 2]];
        bool prev_inside = pos_prev.z >= -pos_prev.w;
        Vec3 tex_prev = tex_in[indices[i + 8]]; // tex coordinate
        
        //int verts_added = 0;
        vector<Vec4> pos_new{};
        vector<Vec3> tex_new{}; // new texture coordinates
        pos_new.reserve(4);
        tex_new.reserve(4);
        //pos_new.clear();
        //texture lerp? and vertex normal lerp needed?
        for (int j = 0; j < 3; ++j) {
            Vec4 pos_cur = pos_clip[indices[i + j]];
            bool cur_inside = pos_cur.z >= -pos_cur.w;
            Vec3 tex_cur = tex_in[indices[i + j + 6]];
            if ( (!prev_inside) != (!cur_inside) ) {

                double lerp_factor = (-pos_prev.w - pos_prev.z) 
                        / ( ( -pos_prev.w - pos_prev.z) 
                        - ( -pos_cur.w - pos_cur.z) );
                pos_new.push_back( pos_prev.lerp(pos_cur, lerp_factor));
                tex_new.push_back( tex_prev.lerp(tex_cur, lerp_factor));
            }
            if ( cur_inside ) {
                pos_new.push_back(pos_cur);
                tex_new.push_back(tex_cur);
            }
            prev_inside = cur_inside;
            pos_prev = pos_cur;
            tex_prev = tex_cur;
        }

        //cout << "DEBUG pos_new.size() = " << pos_new.size() << endl;
        assert(pos_new.size() == 0 || pos_new.size() == 3 || pos_new.size() == 4);
        if (pos_new.size() == 0) {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // CLIPPED!
        } else if ( pos_new.size() == 3) {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // old tri clipped!
            pos_clip.push_back(pos_new[0]);
            pos_clip.push_back(pos_new[1]);
            pos_clip.push_back(pos_new[2]);
            tex_in.push_back(tex_new[0]);
            tex_in.push_back(tex_new[1]);
            tex_in.push_back(tex_new[2]);
            // *END NEW*
            color_faces_world.push_back( color_faces_world[i/NUM_ATTRB_PER_VERT] );
            tris_clipped.push_back(0); // Nnew tri is not clipped!
            // set new indices
            //TODO CHECK WASTED MEM
            // i is index into indices array. 12 elements for each triangle
            indices_clip.resize( indices_clip.size() + NUM_ATTRB_PER_VERT);  
            indices_clip[indices_clip.size() - 12 + 0] = pos_clip.size() - 3; //pos
            indices_clip[indices_clip.size() - 12 + 1] = pos_clip.size() - 2; //pos
            indices_clip[indices_clip.size() - 12 + 2] = pos_clip.size() - 1; //pos
            indices_clip[indices_clip.size() - 12 + 6] = tex_in.size() - 3; //tex
            indices_clip[indices_clip.size() - 12 + 7] = tex_in.size() - 2; //tex
            indices_clip[indices_clip.size() - 12 + 8] = tex_in.size() - 1; //tex
            materials_in.push_back(materials_in[i/NUM_ATTRB_PER_VERT]);
            //indices_clip[indices_clip

        } else if ( pos_new.size() == 4) {
            tris_clipped[i/NUM_ATTRB_PER_VERT] = 1;  // old tri clipped!
            pos_clip.push_back(pos_new[0]);
            pos_clip.push_back(pos_new[1]);
            pos_clip.push_back(pos_new[2]);
            pos_clip.push_back(pos_new[3]);
            tex_in.push_back(tex_new[0]);
            tex_in.push_back(tex_new[1]);
            tex_in.push_back(tex_new[2]);
            tex_in.push_back(tex_new[3]);
            color_faces_world.push_back( color_faces_world[i/NUM_ATTRB_PER_VERT] );
            color_faces_world.push_back( color_faces_world[i/NUM_ATTRB_PER_VERT] );
            tris_clipped.push_back(0); // New tri #1 is not clipped!
            tris_clipped.push_back(0); // New tri #2 is not clipped!
            //TODO CHECK WASTED MEM
            indices_clip.resize( indices_clip.size() + 2 * NUM_ATTRB_PER_VERT);  
            indices_clip[indices_clip.size() - 24 + 0] = pos_clip.size() - 4; //pos
            indices_clip[indices_clip.size() - 24 + 1] = pos_clip.size() - 3; //pos
            indices_clip[indices_clip.size() - 24 + 2] = pos_clip.size() - 2; //pos
            indices_clip[indices_clip.size() - 24 + 6] = tex_in.size() - 4; //tex
            indices_clip[indices_clip.size() - 24 + 7] = tex_in.size() - 3; //tex
            indices_clip[indices_clip.size() - 24 + 8] = tex_in.size() - 2; //tex
            
            // 2nd triangle
            indices_clip[indices_clip.size() - 12 + 0] = pos_clip.size() - 4; //pos
            indices_clip[indices_clip.size() - 12 + 1] = pos_clip.size() - 2; //pos
            indices_clip[indices_clip.size() - 12 + 2] = pos_clip.size() - 1; //pos
            indices_clip[indices_clip.size() - 12 + 6] = tex_in.size() - 4; //tex
            indices_clip[indices_clip.size() - 12 + 7] = tex_in.size() - 2; //tex
            indices_clip[indices_clip.size() - 12 + 8] = tex_in.size() - 1; //tex
            materials_in.push_back(materials_in[i/NUM_ATTRB_PER_VERT]);
            materials_in.push_back(materials_in[i/NUM_ATTRB_PER_VERT]);
        } 
    }
}

void homogenous_divide(int _start_v, int _end_v) {
    // _start and _end represent vertices
    double w = 0.0;
    for (int i = _start_v; i < _end_v; ++i) {
        w = pos_clip[i].w;
        if ( w == 0.0) {
            continue;
        }
        pos_ndc[i] = (1.0 / w) * pos_clip[i] ;
        pos_ndc[i].w = w;
        // keep "old" w in original place, for texture coordinate perspective correction
    }
}

void map_to_window(int _start_v, int _end_v)
{
    //if (!posted_info) { cout << "screen vertices" << endl; }
    //models_out = models_ndc;
    for (int i = _start_v; i < _end_v; ++i) {
        pos_screen[i].x =  1.0 * pos_ndc[i].x * win_res_x 
                / 2.0  + ((double)win_res_x)/2.0;
        pos_screen[i].y =  -1.0 * pos_ndc[i].y * win_res_y 
                / 2.0  + ((double)win_res_y)/2.0;
        pos_screen[i].z = pos_ndc[i].z;
        pos_screen[i].w = pos_ndc[i].w; //1.0; // for texturing?
    }
}

void raster_triangles(int _start_tri, int _end_tri)
{
//    for (int i = 0; i < win_res_y; ++i) {
//        for (int j = 0; j < win_res_x; ++j) {
//            depth_buffer[i][j] = 2.0f;
//        }
//    }
    for (int i = _start_tri * NUM_ATTRB_PER_VERT; 
             i < _end_tri * NUM_ATTRB_PER_VERT; 
             i += NUM_ATTRB_PER_VERT) 
    {
        if (tris_clipped[i/NUM_ATTRB_PER_VERT]) {
            continue; // don't draw clipped triangles
        }
        //if (tris_culled[i/NUM_ATTRB_PER_VERT]) {
        //    continue;
        //}

        Vec4 c = color_faces_world[i / 12];
        assert (c.r <= 1.0);
        assert (c.g <= 1.0);
        assert (c.b <= 1.0);
        c.r = std::max(0.0, c.r);
        c.g = std::max(0.0, c.g);
        c.b = std::max(0.0, c.b);
        //assert (c.r >= 0.0);
        //assert (c.g >= 0.0);
        //assert (c.b >= 0.0);
        SDL_SetRenderDrawColor(gsdl.renderer, 
                int (c.r * 255), //tri.color.r, 
                int (c.g * 255), //tri.color.g,
                int (c.b * 255), //tri.color.b, 
                255);  //20);




// func signature
//void draw_triangle(	int x1, int y1, float w1, float u1, float v1, 
//						int x2, int y2, float w2, float u2, float v2,
//						int x3, int y3, float w3, float u3, float v3,
//	                    SDL_Surface* tex = nullptr)
//

        double w0 = pos_screen[indices_clip[i + 0]].w;
        double w1 = pos_screen[indices_clip[i + 1]].w;
        double w2 = pos_screen[indices_clip[i + 2]].w;

        double z0 = pos_screen[indices_clip[i + 0]].z;
        double z1 = pos_screen[indices_clip[i + 1]].z;
        double z2 = pos_screen[indices_clip[i + 2]].z;
        //draw_triangle(
        auto ab0   =     pos_screen[indices_clip[i + 0]].x; 
        auto ab1   =     pos_screen[indices_clip[i + 0]].y; 
        auto ab2   =     tex_in[indices_clip[i + 6]].x ; // TODO indices_clip[i+6] needs init?
        auto ab3   =     tex_in[indices_clip[i + 6]].y ;
        auto ab4   =     z0; //w0 

        auto ab5   =     pos_screen[indices_clip[i + 1]].x; 
        auto ab6   =     pos_screen[indices_clip[i + 1]].y; 
        auto ab7   =     tex_in[indices_clip[i + 7]].x ; // TODO indices_clip[i+7] needs init?
        auto ab8   =     tex_in[indices_clip[i + 7]].y ;
        auto ab9   =     z1;

        auto ab10  =      pos_screen[indices_clip[i + 2]].x; 
        auto ab11  =      pos_screen[indices_clip[i + 2]].y; 
        auto ab12  =      tex_in[indices_clip[i + 8]].x ; // TODO indices_clip[i+8] needs init?
        auto ab13  =      tex_in[indices_clip[i + 8]].y ;
        auto ab14  =      z2;
        auto ab15  =      color_faces_world[i/NUM_ATTRB_PER_VERT].r; // ugly
        auto ab16  =      materials_in[i/12]->texture_diffuse;

        Vec4 norm_mean = (
                norm_world[indices_clip[i + 3]] +
                norm_world[indices_clip[i + 4]] +
                norm_world[indices_clip[i + 5]] ) / 3.0;


        draw_triangle(
                pos_screen[indices_clip[i + 0]].x, 
                pos_screen[indices_clip[i + 0]].y, 
                tex_in[indices_clip[i + 6]].x , // TODO indices_clip[i+6] needs init?
                tex_in[indices_clip[i + 6]].y ,
                z0, //w0 
                norm_world[indices_clip[i + 3]], // Vec4

                      
                pos_screen[indices_clip[i + 1]].x, 
                pos_screen[indices_clip[i + 1]].y, 
                tex_in[indices_clip[i + 7]].x , // TODO indices_clip[i+7] needs init?
                tex_in[indices_clip[i + 7]].y ,
                z1,
                norm_world[indices_clip[i + 4]], // Vec4
               
                pos_screen[indices_clip[i + 2]].x, 
                pos_screen[indices_clip[i + 2]].y, 
                tex_in[indices_clip[i + 8]].x , // TODO indices_clip[i+8] needs init?
                tex_in[indices_clip[i + 8]].y ,
                z2,
                norm_world[indices_clip[i + 5]], // Vec4
                //color_faces_world[i/NUM_ATTRB_PER_VERT].r, // TODO: ugly light intensity
                light_directional,
                light_ambient,
                materials_in[i/12]->texture_diffuse);

        //draw_triangle(tri.verts[0].x, tri.verts[0].y, tri.verts[0].z,
        //              tri.verts[1].x, tri.verts[1].y, tri.verts[1].z,
        //              tri.verts[2].x, tri.verts[2].y, tri.verts[2].z);

    }
//        for (auto& tri : tris) {
//            if (tri.draw_face && Custom_Draw_Triangle) {
//                // uses custom rasterizer
//
//                //SDL_Color color = { tri.color.r, tri.color.g, tri.color.b, tri.color.a };
//                SDL_SetRenderDrawColor(gsdl.renderer, 
//                    tri.color.r, tri.color.g, tri.color.b ,20);
//                draw_triangle(tri.verts[0].x, tri.verts[0].y, tri.verts[0].z,
//                              tri.verts[1].x, tri.verts[1].y, tri.verts[1].z,
//                              tri.verts[2].x, tri.verts[2].y, tri.verts[2].z);
//
//
//                 
//            } else if (tri.draw_face && !Custom_Draw_Triangle) {
//                draw_triangle(tri); // uses SDL rasterizer
//            }
//            // ignore for now
//            //if (!tri.draw_wire_frame) {
//            //    continue;
//            //}
//            // now draw wireframe
//            // i_tri is index to begining of triangle
//            if (tri.draw_wire_frame) { 
//                for (int i = 0; i < 3; ++i) {
//                    // draws wireframe
//                    // i at vertex of triangle
//                    int j = (i + 1 < 3) ? i + 1 : 0;
//                    SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
//                     
//                    //SDL_RenderDrawLine(gsdl.renderer, 
//                    //        tri.verts[i].data[0], tri.verts[i].data[1],
//                    //        tri.verts[j].data[0], tri.verts[j].data[1]);
//                    if (Custom_Draw_Line) {
//                        draw_line( tri.verts[i].data[0], tri.verts[i].data[1],
//                                   tri.verts[j].data[0], tri.verts[j].data[1]);
//                    } else {
//                        SDL_RenderDrawLine(gsdl.renderer, 
//                                tri.verts[i].data[0], tri.verts[i].data[1],
//                                tri.verts[j].data[0], tri.verts[j].data[1]);
//                    }
//                }
//            }
//        }
//
}

int main() 
{
    // quick test
    cout << "hello world!!!!!!" << endl;
    test_lerp();
    models.reserve(20);
    light_dir.normalize();  // normalize light direction
    //Model<double> model0{"VideoShip.obj", {0, 0, 1, 1}, false}; 
    //Model<double> model1{"teapot.obj", {0, 0, 1, 1}, false}; 
    // flat_export.obj
    //Model<double> model0{"monkey.obj", {0, 0, 1, 1}, false};  
    
    //
    //
    //Model<double> model1{"texture_test.obj", {0, 0, 1, 1}, false};  
    Model<double> model1{"torus_texture.obj", {0, 0, 1, 1}, false};  
    //Model<double> model1{"uvsphere0.obj", {0, 0, 1, 1}, false};  
    //Model<double> model1{"icosphere0.obj", {0, 0, 1, 1}, false};  

    //
    //


    //models.emplace_back(model0);
    models.emplace_back(model1);
    //models[1].meshes[0].trans.set_row(3, Vec<double>{-0.0, 17.0, 3.5, 1.0});
    //models[0].meshes[0].trans.set_row(3, Vec<double>{-0.0, 17.0, 6.5, 1.0});
    models[0].meshes[0].trans.set_row(3, Vec<double>{0.0, 0.0, 0.0, 1.0});

    cam_pos = Vec<double, 3>{0.0, 0.0, -5.0};



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
    cout << "test turn" << endl << test0 << endl;
    // camera space to clip space
    //double near = 1.00;  // 1.0


    //
    //  calculating orthographic and prespective matrices
    //
    double near = 1.00000;  // 1.0
    double far = 20000.0;//20000.0;//200000.0;  // 20.0
    //double projection_width = 0.6;
    //double projection_height = 0.6;
    clip_matrix.set(0, 0,  near );
    clip_matrix.set(1, 1,  near );
    clip_matrix.set(2, 2, (near + far)/(far - near));
    clip_matrix.set(2, 3, 1.0 );
    clip_matrix.set(3, 2, -2.0 * ( near * far)/(far - near) );  // 2.0!!!!!!! TODO
    ////clip_matrix.set(0, 0,  1.0 );
    ////clip_matrix.set(1, 1,  1.0 );
    ////clip_matrix.set(2, 2, (near + far)/(far - near));
    ////clip_matrix.set(2, 3, 1.0 );
    ////clip_matrix.set(3, 2, -2.0*( near * far)/(far - near) ); 
    //clip_matrix.set(0, 0, near / (projection_width * 2));
    //clip_matrix.set(1, 1, near / (projection_height * 2));
    //clip_matrix.set(2, 2, (near + far)/(far - near));
    //clip_matrix.set(2, 3, 1.0 );
    //clip_matrix.set(3, 2, (- 2 * near * far) / (far - near));
    cout << "here's the clip matrix" << endl;
    cout << clip_matrix << endl;
    // ortho matrix set up
    //double r = 3.0, l = -3.0, t = 3.0, b = -3.0;
    double r = 1.0, l = -1.0, t = 1.0, b = -1.0;
    ortho_matrix.set(0, 0, 2.0/(r-l));
    ortho_matrix.set(3, 0, - (r+l)/(r-l));
    ortho_matrix.set(1, 1, 2.0/(t-b));
    ortho_matrix.set(3, 1, - (t+b)/(t-b));
    ortho_matrix.set(2, 2, -2.0/(far-near));
    ortho_matrix.set(3, 2, - (near+far)/(near-far));
    ortho_matrix.set(3, 3, 1.0);

    //
    //  CLIP MATRIX TEST
    //
    cout << "clip matrix test\n";
    cout << "clip matrix test\n";
    cout << "clip matrix test\n";
    Vec<double, 4> v0 = { 0.0, 0.0, 1.0, 1.0}; // ( 0 0 -1 1 )
    cout << v0 * clip_matrix << endl; 
    cout << Vec<double>{0.0, 0.0, far, 1.0} *clip_matrix <<endl;//( 0 0 2000 2000 )
    cout << Vec<double>{0.0, 0.0,(far - near)/2, 1.0} *clip_matrix <<endl;//( 0 0 2000 2000 )
    cout << "------\n";
    cout << Vec<double>{0.0, 0.0,0.0, 1.0} *clip_matrix <<endl;//( 0 0 2000 2000 )
    cout << Vec<double>{0.0, 0.0,near * 1.1, 1.0} *clip_matrix <<endl;  //( 0 0 2000 2000 )
    cout << Vec<double>{0.0, 0.0,near * 2.0, 1.0} *clip_matrix <<endl;  //( 0 0 2000 2000 )
    cout << Vec<double>{0.0, 0.0,near * 3.0, 1.0} *clip_matrix <<endl;  //( 0 0 2000 2000 )

    Vec<double, 4> v1 {0.0, 0.0, -1.0, 1.0}; 
    Vec<double, 4> v2 {0.0, 0.0, 7.0, 1.0}; 
    Vec<double, 4> v3 {0.0, 0.0, 3.0, 1.0}; 
    //Vec<double, 4> v2 {0., 0., 3.0, 1.0}; 
    cout << "v1: " << v1 << endl << "v2 : " << v2 << endl;
    v1 = v1 * clip_matrix;
    v2 = v2 * clip_matrix;
    cout << "v1*clip: " << v1 << endl << "v2*clip : " << v2 << endl;
    cout << "v1 inside(0): " << (v1.z >= -v1.w) << endl;
    cout << "v2 inside(1): " << (v2.z >= -v2.w) << endl;
    //double lerp_factor = (v1.z + v1.w) 
    //       / ( ( v1.z + v1.w - v2.w + v2.z) );
    double lerp_factor = (-v1.w - v1.z) / ( ( -v1.w - v1.z) - ( -v2.w - v2.z) );
    //?works? double lerp_factor = (-v1.w - v1.z) / ( ( -v1.w - v1.z) - ( -v2.w - v2.z) );
    cout << "lerp_factor = " << lerp_factor << endl;
    cout << "lerp v1 to v2 = " << v1.lerp(v2, lerp_factor) << endl;
                    

    cout << " end clip matrix test\n";
    cout << " end clip matrix test\n";
    cout << " end clip matrix test\n";


    using std::cout, std::endl;
    //getchar();
    bool exit = false;    
    cout << "gsdl status: " << gsdl.created_ok << endl;
    SDL_Event e;
    Uint64 time_frame_start;
    Uint64 time_frame_duration = 0;
    mesh_selected = &(models[0].meshes[0]); // model selected for transformations

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
                        mesh_selected->rot = mesh_selected->rot * turn_left;
                    }
                    if (e.key.keysym.sym == SDLK_e) {
                        mesh_selected->rot = mesh_selected->rot * turn_right;
                    }
                    if (e.key.keysym.sym == SDLK_r) {
                        mesh_selected->rot = mesh_selected->rot * turn_up;
                    }
                    if (e.key.keysym.sym == SDLK_f) {
                        mesh_selected->rot = mesh_selected->rot * turn_down;
                    }

                    if (e.key.keysym.sym == SDLK_1) {
                        cout << "1 pressed" << endl;
                        mesh_selected = &(models[0].meshes[0]);
                    }
                    if (e.key.keysym.sym == SDLK_2) {
                        cout << "2 pressed" << endl;
                        mesh_selected = &(models[1].meshes[0]);
                    }
                    if (e.key.keysym.sym == SDLK_3) {
                        cout << "3 pressed" << endl;
                        if (models.size() > 3) { // exclude index 3, the axis!
                            mesh_selected = &(models[3].meshes[0]);
                        }
                    }
                    if (e.key.keysym.sym == SDLK_4) {
                        cout << "4 pressed" << endl;
                        if (models.size() > 4) { // exclude index 3, the axis!
                            mesh_selected = &(models[4].meshes[0]);
                        }
                    }
                    // transform selected model
                    if (e.key.keysym.sym == SDLK_a) { // left
                        mesh_selected->trans.add_to_row(3, Vec{-0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_d) { // right
                        mesh_selected->trans.add_to_row(3, Vec{0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_w) { // forward
                        mesh_selected->trans.add_to_row(3, Vec{0.00,0.0,0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_s) { // backward
                        mesh_selected->trans.add_to_row(3, Vec{0.00,0.0,-0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_z) {  // below
                        mesh_selected->trans.add_to_row(3, Vec{0.0,-0.08,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_x) { // above
                        mesh_selected->trans.add_to_row(3, Vec{0.0,0.08,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_c) { // shrink
                        mesh_selected->scale = mesh_selected->scale * scale_shrink;
                    }
                    if (e.key.keysym.sym == SDLK_v) { // enlarge
                        mesh_selected->scale = mesh_selected->scale * scale_enlarge;
                    }
                    // transform camera
                    if (e.key.keysym.sym == SDLK_j) { // left
                        //auto m = translate_x_minus * rotate_y(cam_yaw);
                        //world_to_camera_translate.add_to_row( 3, m.get_row(3));
                        //world_to_camera_translate.set(3, 3, 1.0);
                        Vec3 x = { -0.1, 0.0, 0.0}; 
                        x = x * rotate_y<double, 3, 3>(cam_yaw);
                        cam_pos += x;

                    }
                    if (e.key.keysym.sym == SDLK_l) { // right
                        //auto m = translate_x_plus * rotate_y(cam_yaw);
                        //world_to_camera_translate.add_to_row( 3, m.get_row(3));
                        //world_to_camera_translate.set(3, 3, 1.0);
                        Vec3 x = { 0.1, 0.0, 0.0}; 
                        x = x * rotate_y<double, 3, 3>(cam_yaw);
                        cam_pos += x;
                    }
                    if (e.key.keysym.sym == SDLK_i) { // forward
                        //Vec4 look_dir = {0.0, 0.0, 1.0, 1.0}; 
                        //look_dir = look_dir * rotate_x(cam_pitch) * rotate_y(cam_yaw);
                        //Vec4 t = 0.08 * look_dir;
                        //world_to_camera_translate.add_to_row(3, t);
                        //world_to_camera_translate.set(3, 3, 1.0);
                        
                        Vec3 look_dir = {0.0, 0.0, 1.0}; 
                        look_dir = look_dir * rotate_x<double, 3, 3>(cam_pitch) * rotate_y<double, 3, 3>(cam_yaw);
                        look_dir = 0.08 * look_dir;
                        cam_pos += look_dir;
                        
                    }
                    if (e.key.keysym.sym == SDLK_k) { // back
                        //Vec<double> look_dir = {0.0, 0.0, 1.0, 1.0}; 
                        //look_dir = look_dir * rotate_x(cam_pitch) * rotate_y(cam_yaw);
                        //auto t = -0.08 * look_dir;
                        //world_to_camera_translate.add_to_row(3, t);
                        //world_to_camera_translate.set(3, 3, 1.0);
                        
                        Vec3 look_dir = {0.0, 0.0, 1.0}; 
                        look_dir = look_dir * rotate_x<double, 3, 3>(cam_pitch) * rotate_y<double, 3, 3>(cam_yaw);
                        look_dir = -0.08 * look_dir;
                        cam_pos += look_dir;

                    }
                    if (e.key.keysym.sym == SDLK_m) { // up
                        //world_to_camera_translate.add_to_row( 
                        //        3, Vec{0.0,0.16,0.0,0.0});
                        cam_pos += Vec3{0.0, 0.16, 0.0};
                    }
                    if (e.key.keysym.sym == SDLK_n) { // down
                        //world_to_camera_translate.add_to_row( 
                        //        3, Vec{0.0, -0.16, 0.0, 0.0});
                        cam_pos += Vec3{0.0, -0.16, 0.0};
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
                    if (e.key.keysym.sym == SDLK_6) {
                        cout << "spawn new ship disabled\n";
                        //Model<double> m;
                        //static double z_pos = 4.0;
                        //m.trans.set(3, 2, z_pos);
                        //z_pos += 4.5;
                        //m.load_from_file("VideoShip.obj");
                        //models.push_back(m);
                    }
                    

                    // toggle draw line code: SDL vs custom 
                    if (e.key.keysym.sym == SDLK_6) {
                        //Custom_Draw_Line = !Custom_Draw_Line;
                    }
                    // toggle draw triangle code: SDL vs custom
                    if (e.key.keysym.sym == SDLK_7) {
                        //Custom_Draw_Triangle = !Custom_Draw_Triangle;
                    }
                    // toggle drawing wire frames
                    if (e.key.keysym.sym == SDLK_8) {
                        //for (auto& model : models) {
                        //    //model.draw_wire_frame = !model.draw_wire_frame;
                        //}
                    }
                    // toggle drawing faces
                    if (e.key.keysym.sym == SDLK_9) {
                        //for (auto& model : models) {
                        //    //model.draw_face = !model.draw_face;
                        //}
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
            cam_yaw = cam_yaw + middle_mouse_down_delta.data[0] * PI / (double)win_res_x;
            // calculate new pitch
            cam_pitch = cam_pitch + middle_mouse_down_delta.data[1] * PI / (double)win_res_y;
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

        // cam_pos and look_at_dir
        
        init_vertices();
        int num_verts = pos_in.size(); // TODO CHANGE
        int num_tris = indices.size()/12; // TODO CHANGE
        assert(indices.size() % 12 == 0);

        Vec3 up = {0.0, 1.0, 0.0};
        //Vec<double, 3> cam_pos = slice<double, 3, 4>(
        //        world_to_camera_translate.get_row(3));
        Vec3 look_dir = {0.0, 0.0, 1.0}; 
        look_dir = look_dir * rotate_x<double, 3, 3>(cam_pitch) * rotate_y<double, 3, 3>(cam_yaw);
        Vec3 target = cam_pos + look_dir;
        //look_at = look_at_matrix(cam_pos, target, up);
        world_to_camera_matrix = look_at_matrix(cam_pos, target, up);


        model_to_world(0, num_verts);
        model_to_world_normals(0, norms_in.size());
        model_to_world_part2(0, num_tris);  // does triangle level operations, obsolete?
        world_to_clip_matrix = world_to_camera_matrix * clip_matrix;
        world_to_clip(0, num_verts);

        clip(0, num_tris);
        // verts and tri counts have changed
        num_verts = pos_clip.size();
        num_tris = indices_clip.size() / NUM_ATTRB_PER_VERT;
        
        // 
        pos_ndc.resize( pos_clip.size() );
        pos_screen.resize( pos_clip.size() ); //crashes
        homogenous_divide(0, num_verts);

        map_to_window(0, num_verts);
        for (int i = 0; i < win_res_y; ++i) {
            for (int j = 0; j < win_res_x; ++j) {
                depth_buffer[i][j] = 2.0f;
            }
        }
        raster_triangles(0, num_tris);
        

        //
        // end rendering pipeline
        //

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
    IMG_Quit();
    SDL_Quit();
}

