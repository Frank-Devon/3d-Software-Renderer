# 3d Software Renderer

![3D](https://github.com/user-attachments/assets/358386ec-6f33-4b92-9480-ba660e4510ed)

**Overview**
The program displays two wireframe models in 3d.  The two models can be translated, rotated, and resized. The camera can also be translated and rotated.

**Controls**
- move camera with: j left, k back, i forward, l right, n down, m up, u turn left, o turn right
- move object with: a left, s back, w forward, d right, z down, x up, q turn left, e turn right, c scale down, v scale up
- There are two models to move, to have the movement controls apply to the first model press 1 (default is 1). To manipulate model 2, press 2.

**What's next**
- load model/vertex data from a file.
- render sides of object and apply global lighting
- and culling and clipping

**Bugs**
- Moving objects outside the view frustrum can cause visual bugs and preformance issues. Proper culling and clipping will fix that.

**Build**
- builds on linux with g++ *.cpp -g -lSDL2 -lSDL2_image -lSDL2_ttf
- requires SDL2 for windowing
