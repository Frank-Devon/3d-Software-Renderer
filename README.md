# 3D Software Renderer

![torus](https://github.com/user-attachments/assets/58ee3aa8-98a6-41e5-8eae-16e81fa5bf03)

**Overview**
The program displays textured objects with smooth diffuse lighting in 3D perspective.  The torus can be translated, rotated, and scaled. The camera can be translated and rotated (yaw and pitch). 

**Controls**
- camera: hold middle mouse and move mouse to pan camera, j left, k back, i forward, l right, n down, m up, u turn left, o turn right, ] turn down, [ turn up. Note while translating the camera, you can pan the camera with the middle mouse.
- object: a left, s back, w forward, d right, z down, x up, q turn left, e turn right, r turn down, f turn up, c scale down, v scale up

**More details**
- Indexed rendering (each vertex transformed only once through the pipeline)
- Texturing and depth buffer used.
- Per-pixel diffuse lighting accomplished by interpolated normals.
- Clipping is done in clip space.
- Row vector's and left hand coordinates are used.

**What's next**
- Multi-threading
- Specular reflections
- Fix minor glitch: perspective correct interpolation of vertex attributes.
- Geometry and texture editing
- General refactoring

**Build**
- builds on linux with g++ *.cpp -g -lSDL2 -lSDL2_image -lSDL2_ttf
- requires SDL2 for windowing

**Credits and Resources used**
- The spaceship, teapot and terrain models are taken from Javidx9 aka OneLoneCoder https://github.com/OneLoneCoder/Javidx9/tree/master/ConsoleGameEngine/BiggerProjects/Engine3D
His youtube playlist for graphics is https://www.youtube.com/watch?v=ih20l3pJoeU . I learned a lot from his channel. His interpolation code was used and edited by me. I removed his perspective correct interpolation fix, but added normal interpolation. This code will be improved later.
- Fundementals of Computer Graphics. Steve Marschner and Peter Shirley.
- 3D Math Primer for Graphics and Game Development by Fletcher Dunn and Ian Parbery
- Cem Yuksel https://www.youtube.com/@cem_yuksel
- Wolfgang Huerst https://www.youtube.com/@huerst
- Keenan Crane https://www.youtube.com/@keenancrane
- thebennybox https://www.youtube.com/@thebennybox
- Jorge Rodriguez https://www.youtube.com/@JorgeVinoRodriguez
