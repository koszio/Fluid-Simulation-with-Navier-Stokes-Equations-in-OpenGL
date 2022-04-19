#pragma once
#include <vector>
#include <iostream>

extern int M;
extern int MAX_DENSITY;

struct CollisionObj {
    int x;
    int y;
    int z;
    int size;

};
class FluidCube {

public:

    float plane_size;
    float dt;
    float diffusion;
    float viscocity;

    std::vector<CollisionObj> collisionCubes = {};

    std::vector<float> temp; // temporary storage
    std::vector<float> density; // density of dye in fluid
    std::vector<float> Vx; // current velocities in x direction
    std::vector<float> Vy; // current velocities in y direction
    std::vector<float> Vz;
    std::vector<float> Vx0; // previous velocities in x direction
    std::vector<float> Vy0; // previous velocities in y direction
    std::vector<float> Vz0;
    FluidCube();
    FluidCube(float dt, float diffusion, float viscosity);
    ~FluidCube(void) { ; };

    void addDensity(int x, int y, int z, float amount);
    void addVelocity(int x, int y, int z, float amountX, float amountY, float amountZ);
    float FluidCube::getDensity(int x, int y, int z);
    void update();
    void randomInit();

    static std::vector<float> diffuse(CollisionObj obj, int b, std::vector<float> x, std::vector<float> x0, float diffusion, float dt, int iter);
    static std::vector<std::vector<float>> project(CollisionObj obj,std::vector<float> Vx, std::vector<float> Vy, std::vector<float> Vz, std::vector<float> p, std::vector<float> div, int iter);
    static std::vector<float> advect(CollisionObj obj,int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& Vx, std::vector<float>& Vy, std::vector<float>& Vz, float dt);
    static std::vector<float> lin_solve(CollisionObj obj,int b, std::vector<float> x, std::vector<float> x0, float a, float c, int iter);
    static std::vector<float> set_bounds(CollisionObj obj,int b, std::vector<float> x);
    void addCollisionCube(int x, int y, int z, int size);
    bool isCollision(int x, int y, int z);
   
};
