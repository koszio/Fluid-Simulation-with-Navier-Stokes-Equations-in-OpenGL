#define WIN32_LEAN_AND_MEAN


#pragma comment(lib, "opengl32.lib")
#pragma comment(lib, "glu32.lib")
#pragma comment(linker, "/subsystem:console")

#include <iostream>
#include "include/gl/glut.h" // OpenGL utilties
#include "windows.h"
#include <vector>
#include <tuple>
#include <algorithm>
#include <math.h>

#include "stdlib.h"
//#include "FluidPlane.h"
#include "FluidCube.h"


// Camera Settings
int rotationY = 25; // press 'r' to rotate around y-axis
int zPos = -5;

// Fluid Container Params
//FluidPlane fluidPlane; 
FluidCube fluidCube;
int angleCounter = 0;
float TIME_STEP = 0.00000001;
float VISCOSITY = 10;
float DIFFUSION = 0.5;
int MAX_DENSITY = 10;
int M = 30;
int HALF_M = round(M / 2.0);
float BASE_VELOCITY = 50;

//prototypes for our callback functions
void DisplayScene(void);                  //our drawing routine
void idle(void);                          //what to do when nothing is happening
void key(unsigned char k, int x, int y);  //handle key presses
void reshape(int width, int height);      //when the window is resized
void init_drawing(void);                  //drawing intialisation

// Basic draw functions
void draw_plane(float x, float y, float z);
void draw_square(float x, float y, float z, float length);
void draw_line(float x0, float y0, float z0, float x1, float y1, float z1);
void draw_cube(float x, float y, float z);

void drawFluidCube(float x, float y, float z);
void drawFluidCubeBounds(void);

int main(int argc, char* argv[]) {
    //Initialise Glut and create a window
    glutInit(&argc, argv);

    //sets up our display mode
    //here we've selected an RGBA display with depth testing
    //and double buffering enabled
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("Fluid Simulation");

    // initialise fluid object instance
    // fluidPlane = FluidPlane(0.000001, 0.25, 0.05);
    fluidCube = FluidCube(TIME_STEP, DIFFUSION, VISCOSITY);
    fluidCube.addCollisionCube(HALF_M, HALF_M, HALF_M, 10);

    //run our own drawing initialisation routine
    init_drawing();

    glutReshapeFunc(reshape);
    glutKeyboardFunc(key);
    glutIdleFunc(idle);
    glutDisplayFunc(DisplayScene);

    //request a window size of 600 x 600
    glutInitWindowSize(600, 600);
    glutReshapeWindow(600, 600);

    glutMainLoop();
    return 0;
}

void DisplayScene(void) {
    //clear the current window
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //make changes to the modelview matridx
    glMatrixMode(GL_MODELVIEW);

    //initialise the modelview matridx to the identity matridx
    glLoadIdentity();

    // Apply camera translations
    glTranslatef(0, 0, zPos);
    glRotatef(rotationY, 0.0, 1.0, 0.0);
    glScalef(0.05, 0.05, 0.05);

    // Draw axis lines (temp)
    //glColor3f(1, 0, 0);
    //draw_line(-1000, 0, 0, 1000, 0, 0);
    //glColor3f(0, 1, 0);
    //draw_line(0, -1000, 0, 0, 1000, 0);
    // z-axis
    //glColor3f(0, 0, 1);
    //draw_line(0, 0, 0, 0, 0, 1000);
    //glColor3f(1, 0, 1);
    //draw_line(0, 0, -1000, 0, 0, 0);

    // Draw current state of fluid plane at pos
    // draw_plane(0, 0, 0);
    // draw_cube(0,0,0);
    drawFluidCube(0,0,0);

    //flush what we've drawn to the buffer
    glFlush();
    //swap the back buffer with the front buffer
    glutSwapBuffers();
}

/*void draw_plane(float x, float y, float z) {
    // origin (x,y,z); (0,0,0)
    // __________
    // | | | | | |
    // | | |O| | |
    // | | | | | |
    // -----------

    std::vector<float> sorted_densities = fluidPlane.density;
    std::sort(sorted_densities.begin(), sorted_densities.end());
    float max_density = sorted_densities.at(sorted_densities.size() -1);
    

    float d = (fluidPlane.plane_size / 2.0);
    for (int i = 0; i < fluidPlane.plane_size; i++) {
        for (int j = 0; j < fluidPlane.plane_size; j++) {

            float pixel_density = fluidPlane.density[IDX(i, j)];
            pixel_density = pixel_density / max_density;

            // Add bounds for pixel density
            if (pixel_density > 1.0) pixel_density = 1.0;
            if (pixel_density < 0.0) pixel_density = 0.0;

            float color = 1 - pixel_density;
            
           
            glColor3f(color, color, 1);
            draw_square(i - d , j - d, z, 2);
        }
    }
}*/

void drawFluidCube(float x, float y, float z) {

    glPushMatrix();
    float d = (fluidCube.plane_size / 2.0);
    glTranslatef(-d, -d, -d);

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    for (auto k = 0; k < fluidCube.plane_size; k++) { // increment closer to screen
        for (auto i = 0; i < fluidCube.plane_size; i++) {
            for (auto j = 0; j < fluidCube.plane_size; j++) {
 

                float pixel_density = fluidCube.getDensity(i,j,k);
                
                if (fluidCube.isCollision(i,j,k)) { 
                    glColor4f(1, 1, 0.0, 1);
                    draw_cube(i, j, k);
                } else if (pixel_density > 0.0) {
                    pixel_density /= MAX_DENSITY;

                    float alpha = pixel_density < 0.25 ? 0.25 : pixel_density;
                    //glColor4f(pixel_density/i, pixel_density/j, pixel_density/k, alpha);
                    glColor4f(1, 1, 1, 0.15 + pixel_density);
                    draw_cube(i, j, k);
                }
            }
        }
    }
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    drawFluidCubeBounds();

    glPopMatrix();
}

void drawFluidCubeBounds() {
    //glEnable(GL_DEPTH_TEST);
    //glDisable(GL_BLEND);

    glLineWidth(1.5);
    glColor3f(1, 1, 1);
    float offset = 0;
    float min = 0 - offset;
    float max = fluidCube.plane_size + offset;

    // closest plane
    draw_line(min, min, min, max, min, min);
    draw_line(min, min, min, min, max, min);
    draw_line(min, max, min, max, max, min);
    draw_line(max, min, min, max, max, min);

    // farthest plane
    draw_line(max, max, max, max, min, max);
    draw_line(max, max, max, min, max, max);
    draw_line(max, min, max, min, min, max);
    draw_line(min, max, max, min, min, max);

    // connectors
    draw_line(min, min, min, min, min, max);
    draw_line(min, max, min, min, max, max);
    draw_line(max, max, min, max, max, max);
    draw_line(max, min, min, max, min, max);
}

void draw_line(float x0, float y0, float z0, float x1, float y1, float z1) {
    glBegin(GL_LINES);
    glVertex3f(x0, y0, z0);
    glVertex3f(x1, y1, z1);
    glEnd();
}

void draw_cube(float x, float y, float z) {
    glPushMatrix();
    glTranslatef(x,y,z);

    //glEnable(GL_DEPTH_TEST);
    glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
    float l = 0.5; // each side of cube has length 1 
    //l -= 0.01; // Add small offset to avoid quads drawn on same coordinates
    //l += 0.01; // Add small offset to make quads intersect

    glVertex3f(l, l, -l);
    glVertex3f(-l, l, -l);
    glVertex3f(-l, l, l);
    glVertex3f(l, l, l);

    // Bottom face (y = -l)
    glVertex3f(l, -l, l);
    glVertex3f(-l, -l, l);
    glVertex3f(-l, -l, -l);
    glVertex3f(l, -l, -l);

    // Front face  (z = l)
    glVertex3f(l, l, l);
    glVertex3f(-l, l, l);
    glVertex3f(-l, -l, l);
    glVertex3f(l, -l, l);

    // Back face (z = -l)
    glVertex3f(l, -l, -l);
    glVertex3f(-l, -l, -l);
    glVertex3f(-l, l, -l);
    glVertex3f(l, l, -l);

    // Left face (x = -l)
    glVertex3f(-l, l, l);
    glVertex3f(-l, l, -l);
    glVertex3f(-l, -l, -l);
    glVertex3f(-l, -l, l);

    // Right face (x = l)
    glVertex3f(l, l, -l);
    glVertex3f(l, l, l);
    glVertex3f(l, -l, l);
    glVertex3f(l, -l, -l);

    glEnd();  // End of drawing color-cube
    glPopMatrix();
}

void draw_square(float x, float y, float z, float length) {
    glBegin(GL_QUADS);
    glVertex3f(x - 0.5f, y + 0.5f, z);
    glVertex3f(x + 0.5f, y + 0.5f, z);
    glVertex3f(x + 0.5f, y - 0.5f, z);
    glVertex3f(x - 0.5f, y - 0.5f, z);
    glEnd();
}

void idle(void) {
    /*angleCounter = (angleCounter + 1) % 360;
    float random = (float)rand() / (float)(RAND_MAX / 360);
    float xMult = cos(0.05 * random);
    float yMult = sin(0.05 * random);
    float baseVelocity = 10000;
    float baseDensity = 1000;

    int c = round(N / 2);
    std::vector<std::tuple<int, int>> emittors = { {c - 1,c - 1}, {c - 1,c}, {c - 1,c + 1}, {c, c - 1}, {c, c}, {c, c + 1}, {c + 1,c - 1}, {c + 1, c}, {c + 1,c + 1} };

    for (int i = 0; i < emittors.size(); i++) {
        int x = std::get<0>(emittors[i]);
        int y = std::get<1>(emittors[i]);
        fluidPlane.addDensity(x, y, baseDensity);
        fluidPlane.addVelocity(x, y, baseVelocity * xMult, baseVelocity * yMult);
    }*/

    //float random = (float)rand() / (float)(RAND_MAX / 1);

    angleCounter = (angleCounter + 1) % 360;
    float random = (float)rand() / (float)(RAND_MAX / 360);
    float xMult = cos(0.00005 * random);
    float yMult = sin(0.00005 * random);
    fluidCube.addDensity(HALF_M - 5, HALF_M, HALF_M, MAX_DENSITY);
    fluidCube.addVelocity(HALF_M - 5, HALF_M, HALF_M, xMult  * BASE_VELOCITY, yMult  * BASE_VELOCITY, 0);

    // eject from beneath cube at an angle
    //fluidCube.addDensity(HALF_M - 4, HALF_M - 4, HALF_M - 4, MAX_DENSITY);
    //fluidCube.addVelocity(HALF_M - 4, HALF_M - 4, HALF_M - 4, velocity, velocity, velocity);

    // eject from side of cube in positive x direction
    
    //fluidCube.addDensity(2, HALF_M, HALF_M, MAX_DENSITY);
    //fluidCube.addVelocity(2, HALF_M, HALF_M, BASE_VELOCITY, 0, 0);
    //fluidCube.addDensity(2, HALF_M+1, HALF_M, MAX_DENSITY);
    //fluidCube.addVelocity(2, HALF_M+1, HALF_M, BASE_VELOCITY, 0, 0);

    // eject from side of cube in negative x direction
    //fluidCube.addDensity(M - 3, HALF_M, HALF_M, MAX_DENSITY);
    //fluidCube.addVelocity(M - 3 , HALF_M, HALF_M, -BASE_VELOCITY, 0, 0);


    fluidCube.update();
    DisplayScene();
    
}

//key callback function - called whenever the user presses a
void key(unsigned char k, int x, int y) {

    switch (k) {
    case 27:  //27 is the ASCII code for the ESCAPE key
        exit(0);
        break;
    case 'f':
        fluidCube.update();
        break;
    case 'd':
        fluidCube.addDensity(2, HALF_M, HALF_M, MAX_DENSITY);
        fluidCube.addVelocity(2, HALF_M, HALF_M, BASE_VELOCITY, 0, 0);
        fluidCube.addDensity(2, HALF_M + 1, HALF_M, MAX_DENSITY);
        fluidCube.addVelocity(2, HALF_M + 1, HALF_M, BASE_VELOCITY, 0, 0);
        //fluidCube.update();
        break;
    case 'r':
        rotationY++;
        break;
   }
    

    DisplayScene();
}

//reshape callback function - called when the window size changed
void reshape(int width, int height) {
    //set the viewport to be the same width and height as the window
    glViewport(0, 0, width, height);

    //make changes to the projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    //set up our projection type
    gluPerspective(45.0, (float)width / (float)height, 1.0, 100.0);
    DisplayScene();
}

//set up OpenGL before we do any drawing
void init_drawing(void) {
    //blend colours across the surface of the polygons
    //another mode to try is GL_FLAT which is flat shading
    glShadeModel(GL_SMOOTH);

    //turn lighting off
    glDisable(GL_LIGHTING);

    // discard faces not visible
    glEnable(GL_CULL_FACE);

    //enable OpenGL hidden surface removal
    //glDepthFunc(GL_LEQUAL);
    //glEnable(GL_DEPTH_TEST);
    //glClearColor(1.0,0.0,0.0,0.0);
}


