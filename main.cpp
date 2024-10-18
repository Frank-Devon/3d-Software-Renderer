#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <exception>
#include "vec.hpp"
#include "test.hpp"
#include "sdl_init.hpp"

using std::cout, std::endl, std::vector, std::array;

// vertexes of object
#define NUM_VERTICES 9
#define NUM_MESHES 2
// clockwise winding means triangle visible, ccw means culled
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
};
// vertices drawn to screen
array<vector<Vec<double>>, NUM_MESHES> model_outputs 
= {
    vector<Vec<double>>{},
    vector<Vec<double>>{}
};


int win_res_x = 1000;
int win_res_y = 1000;
// meshes
//Control_Select control_select = Control_Select::model_0;
enum Control_Select : unsigned int { m_0 = 0, m_1 };
Control_Select control_select = m_0;
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
};

Mat<double> world_to_camera = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.0, 1.0}
};
Mat<double> translate_z_plus = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, 0.1, 1.0}
};
Mat<double> translate_z_minu = {  
    Vec{1.0, 0.0, 0.0, 0.0}, 
    Vec{0.0, 1.0, 0.0, 0.0},
    Vec{0.0, 0.0, 1.0, 0.0},
    Vec{0.0, 0.0, -0.1, 1.0}
};
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


Mat<double> clip_matrix(4, 4); // not used yet, to be used for cliping vertices outside view frustrum.
Mat<double> turn_left(4, 4);          
Mat<double> turn_right(4, 4);          

// clipped vertices go here
vector<Vec<double>> output_vertices;
vector<int> output_vertices_per_polygon;
vector<Vec<double>> temp_vertices;
vector<Vec<double>> clip_planes = {
    Vec{1.0, 0.0, 0.0, -1.0},
    Vec{1.0, 0.0, 0.0, 1.0},
    Vec{0.0, 1.0, 0.0, -1.0},
    Vec{0.0, 1.0, 0.0, 1.0},
    Vec{0.0, 0.0, 1.0, -1.0},
    Vec{0.0, 0.0, 1.0, 1.0}
};

int main() 
{
    int index = 0;
    for(auto& model : models) {
        model_outputs[index] = model; 
        index++;
    }

    // doesn't rotate by 15 degrees, but by 7.5 degrees
    double rad = (PI * 7.5)/360.0;
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

    // camera space to clip space
    double near = 1.00;  // 1.0
    double far = 200000.0;  // 20.0
    double projection_width = 0.6;
    double projection_height = 0.6;
    clip_matrix.set(0, 0, near / (projection_width * 2));
    clip_matrix.set(1, 1, near / (projection_height * 2));
    clip_matrix.set(2, 2, (near + far)/(far - near));
    clip_matrix.set(2, 3, 1.0 );
    clip_matrix.set(3, 2, (- 2 * near * far) / (far - near));
    cout << "here's the clip matrix" << endl;
    cout << clip_matrix << endl;


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
                    // rotate camera
                    if (e.key.keysym.sym == SDLK_u) {
                        Vec<double> translation = world_to_camera.get_row(3);
                        world_to_camera.zero_row(3);
                        world_to_camera = turn_left * world_to_camera;
                        world_to_camera.add_to_row(3, translation);
                    }
                    if (e.key.keysym.sym == SDLK_o) {
                        Vec<double> translation = world_to_camera.get_row(3);
                        world_to_camera.zero_row(3);
                        world_to_camera = turn_right * world_to_camera;
                        world_to_camera.add_to_row(3, translation);
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
                    if (e.key.keysym.sym == SDLK_j) {
                        world_to_camera.add_to_row( 
                                3, Vec{0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_l) {
                        world_to_camera.add_to_row( 
                                3, Vec{-0.08,0.0,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_i) {
                        world_to_camera.add_to_row( 
                                3, Vec{0.0,0.0,-0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_k) {
                        world_to_camera.add_to_row( 
                                3, Vec{0.0,0.0,0.16,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_n) {
                        world_to_camera.add_to_row( 
                                3, Vec{0.0,0.16,0.0,0.0});
                    }
                    if (e.key.keysym.sym == SDLK_m) {
                        world_to_camera.add_to_row( 
                                3, Vec{0.0,-0.16,0.0,0.0});
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
                        cout << "control_select " << control_select << endl;
                        model_transforms[control_select] = model_scale[control_select] 
                            * model_rotate[control_select] * model_translate[control_select];
                        cout << "model_transform" << endl;
                        cout << model_transforms[control_select] << endl;
                        cout << "model_scale" << endl;
                        cout << model_scale[control_select] << endl;
                        cout << "model_rotate" << endl;
                        cout << model_rotate[control_select] << endl;
                        cout << "model_translate" << endl;
                        cout << model_translate[control_select] << endl;
                        cout << "scale_shrink" << endl;
                        cout << scale_shrink<< endl;
                        //cout << "world_to_camera" << endl;
                        //cout << world_to_camera << endl;
                        //cout << "clip_matrix" << endl;
                        //cout << clip_matrix << endl;

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


        //print model vertices
        if (!posted_info) { 
            cout << "object space vertices" << endl;
        }
        //model to world transform
        if (!posted_info) { cout << "world space vertices" << endl; }
        int model_index = 0;
        for (auto& model : models) {
            // compute model_transform for given translation, rotation and scale
            model_transforms[model_index] = model_scale[model_index] 
                * model_rotate[model_index] * model_translate[model_index];
            for (int i = 0; i < model.size(); ++i) {
                model_outputs[model_index][i] = model[i] * model_transforms[model_index];
            }
            model_index++;
        }
        //world to camera
        if (!posted_info) { cout << "camera space vertices" << endl; }
        model_index = 0;
        for (auto& model : model_outputs) {
            for (int i = 0; i < model.size(); ++i) {
                model[i] = model[i] * world_to_camera;
            }
            model_index++;
        }
        //camera to clip
        model_index = 0;
        if (!posted_info) { 
            cout << "clip matrix " << endl << clip_matrix << endl;
            cout << "clip space vertices" << endl; }
        for (auto& model : model_outputs) {
            for (int i = 0; i < model.size(); ++i) {
                model[i] = model[i] * clip_matrix;
            }
            model_index++;
        }
        // homogenous divide
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

        if (!posted_info) { cout << "screen vertices" << endl; }

        for (auto& model : model_outputs) {
            for (int i = 0; i < model.size(); ++i) {
                model[i].data[0] = model[i].data[0] * win_res_x
                    / (model[i].data[2] * 2.0) + 1000.0/2.0;
                model[i].data[1] = -1.0 * model[i].data[1] * win_res_y
                    / (model[i].data[2] * 2.0) + 1000.0/2.0;
            }
        }

        // draw lines between vertices
        for (auto& model : model_outputs) {
            for (int i_tri = 0; i_tri < model.size(); i_tri += 3) {
                // i_tri is index to begining of triangle
                for (int i = i_tri; i < i_tri + 3; ++i) {
                    // i at vertex of triangle
                    int j = (i + 1 < i_tri + 3) ? i + 1 : i_tri;
                    //// give each edge a different color
                    //if (i % 3 == 0) { // red
                    //    SDL_SetRenderDrawColor(gsdl.renderer, 150, 0, 0, 20);
                    //} else if (i % 3 == 1) { // green
                    //    SDL_SetRenderDrawColor(gsdl.renderer, 0, 150, 0, 20);
                    //} else if (i % 3 == 2) { // blue
                    //    SDL_SetRenderDrawColor(gsdl.renderer, 0, 0, 150, 20);
                    //}
                    SDL_SetRenderDrawColor(gsdl.renderer, 250, 250, 250, 20);
                    SDL_RenderDrawLine(gsdl.renderer, 
                            model[i].data[0], model[i].data[1],
                            model[j].data[0], model[j].data[1]);
                }
            }
        }

        if (posted_info == false) {
            posted_info = true;
        }


        SDL_RenderPresent(gsdl.renderer);
        time_frame_duration = SDL_GetTicks64() - time_frame_start;
        if (time_frame_duration < 1000/30) {
            SDL_Delay(1000/30 - time_frame_duration);
        }
    }
    SDL_DestroyWindow(gsdl.window);
    SDL_Quit();
}

