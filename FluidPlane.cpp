#include "FluidPlane.h"

int IX(int x, int y) {
	return x + y * N;
}

void print(int x, int y, std::vector<float> arr, char* name) {
	bool print = false;
	if (print) std::cout << name << ": " << arr[IX(x, y)] << std::endl;
}

FluidPlane::FluidPlane() {
	this->plane_size = 10;
	this->dt = 0.001;
	this->diffusion = 0.01;
	this->viscocity = 0.01; 

	int n = this->plane_size * this->plane_size;
	this->temp = { 0 };
	this->density = { 0 };
	this->Vx = { 0 };
	this->Vy = { 0 };
	this->Vx0 = { 0 };
	this->Vy0 = { 0 };
}

FluidPlane::FluidPlane(float dt, float diffusion, float viscosity) {
	this->dt = dt;
	this->diffusion = diffusion;
	this->viscocity = viscosity;

	int n = this->plane_size * this->plane_size;
	this->temp.resize(n, 0);
	this->density.resize(n, 0);
	this->Vx.resize(n, 0);
	this->Vy.resize(n, 0);
	this->Vx0.resize(n, 0);
	this->Vy0.resize(n, 0);
};

void FluidPlane::randomInit() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			float random = (float)rand() / (float)(RAND_MAX / 10000);
			this->addDensity(i, j, random);
			this->addVelocity(i, j, random, random);
		}
	}
}

// Add density to the fluid, this represents the dye
// which is what we can observe
void FluidPlane::addDensity(int x, int y, float amount) {
	this->density[IX(x,y)] += amount;
}

// If there is no velocity the fluid will not move around
// much, which would be a dull simulation.
void FluidPlane::addVelocity(int x, int y, float amountX, float amountY) {
	this->Vx[IX(x,y)] += amountX;
	this->Vy[IX(x,y)] += amountY;
}

// An update step will do all the operations necessary to
// calculate the state of the fluid plane at the next timestep
void FluidPlane::update() {

	std::vector<std::vector<float>> projectResult;

	//print(50, 50, this->Vx, "Vx Before Diffuse");
	//print(50, 50, this->Vx0, "Vx0 Before Diffuse");

	// diffuse along x and y axes, velocity components
	this->Vx0 = diffuse(1, this->Vx0, this->Vx, this->viscocity, this->dt, 4);
	this->Vy0 = diffuse(2, this->Vy0, this->Vy, this->viscocity, this->dt, 4);

	//print(50, 50, this->Vx, "Vx After Diffuse");
	//print(50, 50, this->Vx0, "Vx0 After Diffuse");

	// project things, fidx up velocities which maintains incompressibility
	projectResult = project(this->Vx0, this->Vy0, this->Vx, this->Vy, 4);
	this->Vx0 = projectResult[0];
	this->Vy0 = projectResult[1];
	this->Vx = projectResult[2];
	this->Vy = projectResult[3];

	//print(50, 50, this->Vx, "Vx After Project");
	//print(50, 50, this->Vx0, "Vx0 After Project");

	// advect, move velocities around according to velocities of fluid
	this->Vx = advect(1, this->Vx, this->Vx0, this->Vx0, this->Vy0, this->dt);
	this->Vy = advect(2, this->Vy, this->Vy0, this->Vx0, this->Vy0, this->dt);

	//print(50, 50, this->Vx, "Vx After Advect");
	//print(50, 50, this->Vx0, "Vx0 After Advect");

	// project, fidx up velocities again
	projectResult = project(this->Vx, this->Vy, this->Vx0, this->Vy0, 4);
	this->Vx = projectResult[0];
	this->Vy = projectResult[1];
	this->Vx0 = projectResult[2];
	this->Vy0 = projectResult[3];

	//print(50, 50, this->Vx, "Vx After 2nd Project");
	//print(50, 50, this->Vx0, "Vx0 After 2nd Project");

	// diffuse density, diffuse the dye
	this->temp = diffuse(0, this->temp, this->density, this->diffusion, this->dt, 4);

	// advect density, move dye according to velocities
	this->density = advect(0, this->density, this->temp, this->Vx, this->Vy, this->dt);
}

std::vector<float> FluidPlane::advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& Vx, std::vector<float>& Vy, float dt) {
	float i0, i1, j0, j1, k0, k1;

	float dtx = dt * (N - 2);
	float dty = dt * (N - 2);

	float s0, s1, t0, t1, u0, u1;
	float tmp1, tmp2, tmp3, x, y;

	float Nfloat = N;
	float ifloat, jfloat;
	int i, j;

	for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
		for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
			tmp1 = dtx * Vx[IX(i, j)];
			tmp2 = dty * Vy[IX(i, j)];
			x = ifloat - tmp1;
			y = jfloat - tmp2;


			if (x < 0.5f) x = 0.5f;
			if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
			i0 = floorf(x);
			i1 = i0 + 1.0f;
			if (y < 0.5f) y = 0.5f;
			if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
			j0 = floorf(y);
			j1 = j0 + 1.0f;

			s1 = x - i0;
			s0 = 1.0f - s1;
			t1 = y - j0;
			t0 = 1.0f - t1;

			int i0i = i0;
			int i1i = i1;
			int j0i = j0;
			int j1i = j1;

			d[IX(i, j)] =
				s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)])
				+ s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
		}
	}

	return FluidPlane::set_bounds(b, d);
}

std::vector<float> FluidPlane::diffuse(int b, std::vector<float> x, std::vector<float> x0, float diffusion, float dt, int iter) {
	float a = dt * diffusion * (N - 2) * (N - 2); // 0.1 * 0.001 * 13 * 13 = 0.0169
	return FluidPlane::lin_solve(b, x, x0, a, 1 + 6 * a, iter);
}

std::vector<std::vector<float>> FluidPlane::project(std::vector<float> Vx, std::vector<float> Vy, std::vector<float> p, std::vector<float> div, int iter) {
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			div[IX(i, j)] = -0.5f * (Vx[IX(i + 1, j)] - Vx[IX(i - 1, j)] + Vy[IX(i, j + 1)] - Vy[IX(i, j - 1)]) / N;
			p[IX(i, j)] = 0;
		}
	}

	div = FluidPlane::set_bounds(0, div);
	p = FluidPlane::set_bounds(0, p);
	p = FluidPlane::lin_solve(0, p, div, 1, 6, iter);

	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			Vx[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
			Vy[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
		}
	}

	Vx = FluidPlane::set_bounds(1, Vx);
	Vy = FluidPlane::set_bounds(2, Vy);

	std::vector<std::vector<float>> result(4);
	result[0] = Vx;
	result[1] = Vy;
	result[2] = p;
	result[3] = div;

	return result;
}

std::vector<float> FluidPlane::lin_solve(int b, std::vector<float> x, std::vector<float> x0, float a, float c, int iter) {
	float cRecip = 1.0 / c;
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			float result = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)]));
			x[IX(i, j)] = result * cRecip;
		}
	}
	return FluidPlane::set_bounds(b, x);
}

std::vector<float> FluidPlane::set_bounds(int b, std::vector<float> x) {
	// top row and bottom row
	for (int i = 1; i < N - 1; i++) {
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
	}
	// left col and right col
	for (int j = 1; j < N - 1; j++) {
		x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
		x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
	}

	// corner cells are equal to the avg of surrounding cells
	// (0,0), (0,N-1), (N-1,0), (N-1,N-1)
	x[IX(0, 0)] = 0.33f * (x[IX(1, 0)] + x[IX(0, 1)] + x[IX(0, 0)]);
	x[IX(0, N - 1)] = 0.33f * (x[IX(1, N - 1)] + x[IX(0, N - 2)] + x[IX(0, N - 1)]);
	x[IX(N - 1, 0)] = 0.33f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)] + x[IX(N - 1, 0)]);
	x[IX(N - 1, N - 1)] = 0.33f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)] + x[IX(N - 1, N - 1)]);

	return x;
}
