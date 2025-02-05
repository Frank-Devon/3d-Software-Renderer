# 3D Software Renderer

![graphics_demo](https://github.com/user-attachments/assets/25486ebc-72e5-484a-a86f-6cbeb0416f93)


**Overview**
The program displays objects in 3D perspective.  The teapot and spaceship can be translated, rotated, and scaled. The camera can be translated and rotated (yaw and pitch). Basic lighting is implemented (A simple directional lighting model is used). 

**Controls**
- camera: hold middle mouse and move mouse to pan camera, j left, k back, i forward, l right, n down, m up, u turn left, o turn right, ] turn down, [ turn up. Note while translating the camera, you can pan the camera with the middle mouse.
- objects: a left, s back, w forward, d right, z down, x up, q turn left, e turn right, r turn down, f turn up, c scale down, v scale up
- There are two models to move, to have the movement controls apply to the first model press 1 (default is 1). To manipulate model 2, press 2.
- toggle wireframes 8, toggle faces 9, toggle back face culling 0

**More details**
- Back face culling enabled (CCW triangles are culled)
- Clipping on the near plane is done in clip space.
- Row vector's and left hand coordinates are used.
- Painters algorithm used

**What's next**
- Vec type is very slow. All of it's memory is stored on the heap. Static allocation will be used. This should increase speed greatly.
- Texturing
- Depth buffer
- All rasterizing done with draw pixel function of SDL, no more draw line and draw triangle function calls.
- Per pixel lighting
- General refactoring
- More geometric primitives. Add more information to vertices, like normal's.
- Clip on Left, Right, Top, Bottom, and Far planes. Currently only clipping near plane.
- Improve rotations, add shear transform to objects. add banking to camera
- Multi-threading

**Build**
- builds on linux with g++ *.cpp -g -lSDL2 -lSDL2_image -lSDL2_ttf
- requires SDL2 for windowing

**Credits and Resources used**
- The spaceship, teapot and terrain models are taken from Javidx9 aka (
OneLoneCoder) https://github.com/OneLoneCoder/Javidx9/tree/master/ConsoleGameEngine/BiggerProjects/Engine3D
His youtube playlist for graphics is https://www.youtube.com/watch?v=ih20l3pJoeU . I learned a lot from his channel.
- 3D Math Primer for Graphics and Game Development by Fletcher Dunn and Ian Parbery
- Cem Yuksel https://www.youtube.com/@cem_yuksel
- Wolfgang Huerst https://www.youtube.com/@huerst
- Keenan Crane https://www.youtube.com/@keenancrane
- thebennybox https://www.youtube.com/@thebennybox
- Jorge Rodriguez https://www.youtube.com/@JorgeVinoRodriguez
