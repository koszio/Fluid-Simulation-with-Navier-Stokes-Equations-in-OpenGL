#include <vector>
#include <iostream>

// number of surfaces along each axis
const int N = 120;

class FluidPlane {

public:

    int plane_size = N;
    float dt;
    float diffusion;
    float viscocity;

    std::vector<float> temp; // temporary storage
    std::vector<float> density; // density of dye in fluid
    std::vector<float> Vx; // current velocities in x direction
    std::vector<float> Vy; // current velocities in y direction
    std::vector<float> Vx0; // previous velocities in x direction
    std::vector<float> Vy0; // previous velocities in y direction

    FluidPlane();
    FluidPlane(float dt, float diffusion, float viscosity);
    ~FluidPlane(void) { ; };

    void addDensity(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);
    void update();
    void randomInit();

    static std::vector<float> diffuse(int b, std::vector<float> x, std::vector<float> x0, float diffusion, float dt, int iter);
    static std::vector<std::vector<float>> project(std::vector<float> Vx, std::vector<float> Vy, std::vector<float> p, std::vector<float> div, int iter);
    static std::vector<float> advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& Vx, std::vector<float>& Vy, float dt);
    static std::vector<float> lin_solve(int b, std::vector<float> x, std::vector<float> x0, float a, float c, int iter);
    static std::vector<float> set_bounds(int b, std::vector<float> x);

};