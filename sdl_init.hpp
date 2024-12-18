#ifndef SDL_INIT_HPP
#define SDL_INIT_HPP

#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_ttf.h>
#include <iostream>
//#include "vector.hpp"

class SDL_Main {
public:
    SDL_Window *window;
    SDL_Renderer *renderer;
    TTF_Font *font;
    int window_width;
    int window_height;
    bool created_ok;
    SDL_Main(int screen_width, int screen_height);
};

extern SDL_Main gsdl;

#endif



