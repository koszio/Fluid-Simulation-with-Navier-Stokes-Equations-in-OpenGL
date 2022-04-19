#include "FluidCube.h"

int IX(int x, int y,int z) {
	return x + y * M + z * M * M;
}


FluidCube::FluidCube() {
	this->plane_size = 10;
	this->dt = 0.001;
	this->diffusion = 0.01;
	this->viscocity = 0.01;
	int n = this->plane_size * this->plane_size * this->plane_size;
	this->temp = { 0 };
	this->density = { 0 };
	this->Vx = { 0 };
	this->Vy = { 0 };
	this->Vz = { 0 };
	this->Vx0 = { 0 };
	this->Vy0 = { 0 };
	this->Vz0 = { 0 };
}

FluidCube::FluidCube(float dt, float diffusion, float viscosity) {
	this->plane_size = M;
	this->dt = dt;
	this->diffusion = diffusion;
	this->viscocity = viscosity;

	int n = this->plane_size * this->plane_size* this->plane_size;
	this->temp.resize(n, 0);
	this->density.resize(n, 0);
	this->Vx.resize(n, 0);
	this->Vy.resize(n, 0);
	this->Vz.resize(n, 0);
	this->Vx0.resize(n, 0);
	this->Vy0.resize(n, 0);
	this->Vz0.resize(n, 0);
}


void FluidCube::addDensity(int x, int y, int z, float amount) {
	if (amount > MAX_DENSITY) amount = MAX_DENSITY;
	this->density[IX(x, y, z)] += amount;
}

// If there is no velocity the fluid will not move around
// much, which would be a dull simulation.
void FluidCube::addVelocity(int x, int y, int z,float amountX, float amountY, float amountZ) {
	this->Vx[IX(x, y,z)] += amountX;
	this->Vy[IX(x, y,z)] += amountY;
	this->Vz[IX(x, y,z)] += amountZ;
}

void FluidCube::update() {

	std::vector<std::vector<float>> projectResult;
	CollisionObj obj = this->collisionCubes[0];
	// diffuse along x and y axes, velocity components
	this->Vx0 = diffuse(obj, 1, this->Vx0, this->Vx, this->viscocity, this->dt, 4);
	this->Vy0 = diffuse(obj, 2, this->Vy0, this->Vy, this->viscocity, this->dt, 4);
	this->Vz0 = diffuse(obj, 3, this->Vz0, this->Vz, this->viscocity, this->dt, 4);

	// project things, fidx up velocities which maintains incompressibility
	projectResult = project(obj, this->Vx0, this->Vy0, this->Vz0, this->Vx, this->Vy, 4);
	this->Vx0 = projectResult[0];
	this->Vy0 = projectResult[1];
	this->Vz0 = projectResult[2];
	this->Vx = projectResult[3];
	this->Vy = projectResult[4];
	

	// advect, move velocities around according to velocities of fluid
	this->Vx = advect(obj,1, this->Vx, this->Vx0, this->Vx0, this->Vy0, this->Vz0,this->dt);
	this->Vy = advect(obj,2, this->Vy, this->Vy0, this->Vx0, this->Vy0, this->Vz0, this->dt);
	this->Vz = advect(obj,3, this->Vz, this->Vz0, this->Vx0, this->Vy0, this->Vz0, this->dt);

	// project, fidx up velocities again
	projectResult = project(obj,this->Vx, this->Vy,this->Vz, this->Vx0, this->Vy0, 4);
	this->Vx = projectResult[0];
	this->Vy = projectResult[1];
	this->Vz = projectResult[2];
	this->Vx0 = projectResult[3];
	this->Vy0 = projectResult[4];

	// diffuse density, diffuse the dye
	this->temp = diffuse(obj, 0, this->temp, this->density, this->diffusion, this->dt, 4);

	// advect density, move dye according to velocities
	this->density = advect(obj,0, this->density, this->temp, this->Vx, this->Vy, this->Vz, this->dt);
}



std::vector<float> FluidCube::diffuse(CollisionObj obj,int b, std::vector<float> x, std::vector<float> x0, float diffusion, float dt, int iter) {
	float a = dt * diffusion * (M - 2) * (M - 2); // 0.1 * 0.001 * 13 * 13 = 0.0169
	return FluidCube::lin_solve(obj,b, x, x0, a, 1 + 6 * a, iter);
}

std::vector<std::vector<float>> FluidCube::project(CollisionObj obj,std::vector<float> Vx, std::vector<float> Vy, std::vector<float> Vz, std::vector<float> p, std::vector<float> div, int iter) {

	for (int k = 1; k < M - 1; k++) {
		for (int j = 1; j < M - 1; j++) {
			for (int i = 1; i < M - 1; i++) {

				div[IX(i, j, k)] = -0.5f * (Vx[IX(i + 1, j, k)] - Vx[IX(i - 1, j, k)] + Vy[IX(i, j + 1, k)] - Vy[IX(i, j - 1, k)] + Vy[IX(i, j, k + 1)] - Vy[IX(i, j, k - 1)]) / M;
				p[IX(i, j, k)] = 0;
			}
		}
	}

	div = FluidCube::set_bounds(obj, 0, div);
	p = FluidCube::set_bounds(obj, 0, p);
	p = FluidCube::lin_solve(obj, 0, p, div, 1, 6, iter);
	for (int k = 1; k < M - 1; k++){
		for (int j = 1; j < M - 1; j++) {
			for (int i = 1; i < M - 1; i++) {
			 
				Vx[IX(i, j, k)] -= 0.5f * (p[IX(i + 1, j,k)] - p[IX(i - 1, j,k)]) * M;
				Vy[IX(i, j, k)] -= 0.5f * (p[IX(i, j + 1,k)] - p[IX(i, j - 1,k)]) * M;
				Vy[IX(i, j, k)] -= 0.5f * (p[IX(i, j,k + 1)] - p[IX(i, j,k - 1)]) * M;
			}
		}
	}
	Vx = FluidCube::set_bounds(obj, 1, Vx);
	Vy = FluidCube::set_bounds(obj, 2, Vy);
	Vz = FluidCube::set_bounds(obj, 3, Vz);

	std::vector<std::vector<float>> result(5);
	result[0] = Vx;
	result[1] = Vy;
	result[2] = Vz;
	result[3] = p;
	result[4] = div;

	return result;
}


std::vector<float> FluidCube::advect(CollisionObj obj, int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& Vx, std::vector<float>& Vy, std::vector<float>& Vz, float dt) {
	float i0, i1, j0, j1, k0, k1;

	float dtx = dt * (M - 2);
	float dty = dt * (M -2);
	float dtz = dt * (M - 2);

	float s0, s1, t0, t1, u0, u1;
	float tmp1, tmp2, tmp3, x, y, z;

	float Nfloat = M;
	float ifloat, jfloat, kfloat;
	int i, j,k;
	for (k = 1, kfloat = 1; k < M - 1; k++, kfloat++)
		for (j = 1, jfloat = 1; j < M - 1; j++, jfloat++) {
			for (i = 1, ifloat = 1; i < M - 1; i++, ifloat++) {
				tmp1 = dtx * Vx[IX(i, j, k)];
				tmp2 = dty * Vy[IX(i, j, k)];
				tmp3 = dtz * Vz[IX(i, j, k)];
				x = ifloat - tmp1;
				y = jfloat - tmp2;
				z = kfloat - tmp3;


				if (x < 0.5f) x = 0.5f;
				if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
				i0 = floorf(x);
				i1 = i0 + 1.0f;
				if (y < 0.5f) y = 0.5f;
				if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
				j0 = floorf(y);
				j1 = j0 + 1.0f;

				if (z < 0.5f) z = 0.5f;
				if (z > Nfloat + 0.5f) z = Nfloat + 0.5f;
				k0 = floorf(z);
				k1 = k0 + 1.0f;

				s1 = x - i0;
				s0 = 1.0f - s1;
				t1 = y - j0;
				t0 = 1.0f - t1;
				u1 = z - k0;
				u0 = 1.0f - u1;

				int i0i = i0;
				int i1i = i1;
				int j0i = j0;
				int j1i = j1;
				int k0i = k0;
				int k1i = k1;

				d[IX(i, j, k)] =

					s0 * (t0 * (u0 * d0[IX(i0i, j0i, k0i)]
						+ u1 * d0[IX(i0i, j0i, k1i)])
						+ (t1 * (u0 * d0[IX(i0i, j1i, k0i)]
							+ u1 * d0[IX(i0i, j1i, k1i)])))
					+ s1 * (t0 * (u0 * d0[IX(i1i, j0i, k0i)]
						+ u1 * d0[IX(i1i, j0i, k1i)])
						+ (t1 * (u0 * d0[IX(i1i, j1i, k0i)]
							+ u1 * d0[IX(i1i, j1i, k1i)])));
			}
		}
	

	return FluidCube::set_bounds(obj, b, d);
}


std::vector<float> FluidCube::lin_solve(CollisionObj obj, int b, std::vector<float> x, std::vector<float> x0, float a, float c, int iter) {
	float cRecip = 1.0 / c;
	for (int k = 0; k < iter; k++) {
		for (int m = 1; m < M - 1; m++) {
			for (int j = 1; j < M - 1; j++) {
				for (int i = 1; i < M - 1; i++) {
					x[IX(i, j, m)] =
						(x0[IX(i, j, m)]
							+ a * (x[IX(i + 1, j, m)]
								+ x[IX(i - 1, j, m)]
								+ x[IX(i, j + 1, m)]
								+ x[IX(i, j - 1, m)]
								+ x[IX(i, j, m + 1)]
								+ x[IX(i, j, m - 1)]
								)) * cRecip;
				}
			}
		}
	}
	return FluidCube::set_bounds(obj, b, x);
}

std::vector<float> FluidCube::set_bounds(CollisionObj obj, int b, std::vector<float> x) {

	// Handle inverted component along z-axis, front and back plane
	for (int j = 1; j < M - 1; j++) {
		for (int i = 1; i < M - 1; i++) {
			x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
			x[IX(i, j, M - 1)] = b == 3 ? -x[IX(i, j, M - 2)] : x[IX(i, j, M - 2)];
		}
	}

	// Handle inverted component along y-axis, top and bottom plane
	for (int k = 1; k < M - 1; k++) {
		for (int i = 1; i < M - 1; i++) {
			x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
			x[IX(i, M - 1, k)] = b == 2 ? -x[IX(i, M - 2, k)] : x[IX(i, M - 2, k)];
		}
	}
	// Handle inverted component along x-axis, left and right plane
	for (int k = 1; k < M - 1; k++) {
		for (int j = 1; j < M - 1; j++) {
			x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
			x[IX(M - 1, j, k)] = b == 1 ? -x[IX(M - 2, j, k)] : x[IX(M - 2, j, k)];
		}
	}

	// Handle inverted component along x-axis for collision object
	for (auto j = obj.y; j < obj.y + obj.size; j++) {
		for (auto k = obj.z; k < obj.z + obj.size; k++) {

			// invert velocity coming from positive x direction against collision cube
			// if b !== 1, the velocity component is not going towards collision cube (y || z)
			x[IX(obj.x - 1, j, k)] = b == 1 ? -x[IX(obj.x - 2, j, k)] : x[IX(obj.x - 2, j, k)];

			// invert velocity coming from negative x direction against collision cube
			// if b !== 1, the velocity component is not going towards collision cube (y || z)
			x[IX(obj.x + obj.size, j, k)] = b == 1 ? -x[IX(obj.x + obj.size + 1, j, k)] : x[IX(obj.x + obj.size + 2, j, k)];
		}
	}

	// Handle inverted component along y-axis for collision object
	for (auto i = obj.x; i < obj.x + obj.size; i++) {
		for (auto k = obj.z; k < obj.z + obj.size; k++) {

			// invert velocity coming from positive y direction against collision cube
			// if b !== 2, the velocity component is not going towards collision cube (x || z)
			x[IX(i, obj.y - 1, k)] = b == 2 ? -x[IX(i, obj.y - 2, k)] : x[IX(i, obj.y - 2, k)];

			// invert velocity coming from negative y direction against collision cube
			// if b !== 2, the velocity component is not going towards collision cube (x || z)
			x[IX(i, obj.y + obj.size, k)] = b == 2 ? -x[IX(i, obj.y + obj.size + 1, k)] : x[IX(i, obj.y + obj.size + 1, k)];
		}
	}

	// Handle inverted component along z-axis for collision object
	for (auto j = obj.y; j < obj.y + obj.size; j++) {
		for (auto i = obj.x; i < obj.x + obj.size; i++) {

			// invert velocity coming from positive z direction against collision cube
			// if b !== 3, the velocity component is not going towards collision cube (y || x)
			x[IX(i, j, obj.z - 1)] = b == 3 ? -x[IX(i, j, obj.z - 2)] : x[IX(i, j, obj.z - 2)];

			// invert velocity coming from negative z direction against collision cube
			// if b !== 3, the velocity component is not going towards collision cube (y || x)
			x[IX(i, j, obj.z + obj.size)] = b == 3 ? -x[IX(i, j, obj.z + obj.size + 1)] : x[IX(i, j, obj.z + obj.size + 1)];
		}
	}




	
	x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)]
		+ x[IX(0, 1, 0)]
		+ x[IX(0, 0, 1)]);
	x[IX(0, M - 1, 0)] = 0.33f * (x[IX(1, M - 1, 0)]
		+ x[IX(0, M - 2, 0)]
		+ x[IX(0, M - 1, 1)]);
	x[IX(0, 0, M - 1)] = 0.33f * (x[IX(1, 0, M - 1)]
		+ x[IX(0, 1, M - 1)]
		+ x[IX(0, 0, M-1)]);
	x[IX(0, M - 1, M - 1)] = 0.33f * (x[IX(1, M - 1, M - 1)]
		+ x[IX(0, M - 2, M - 1)]
		+ x[IX(0, M - 1, M - 2)]);
	x[IX(M - 1, 0, 0)] = 0.33f * (x[IX(M - 2, 0, 0)]
		+ x[IX(M - 1, 1, 0)]
		+ x[IX(M - 1, 0, 1)]);
	x[IX(M - 1, M - 1, 0)] = 0.33f * (x[IX(M - 2, M - 1, 0)]
		+ x[IX(M - 1, M - 2, 0)]
		+ x[IX(M - 1, M - 1, 1)]);
	x[IX(M - 1, 0, M - 1)] = 0.33f * (x[IX(M - 2, 0, M - 1)]
		+ x[IX(M - 1, 1, M - 1)]
		+ x[IX(M - 1, 0, M - 2)]);
	x[IX(M - 1, M - 1, M - 1)] = 0.33f * (x[IX(M - 2, M - 1, M - 1)]
		+ x[IX(M - 1, M - 2, M - 1)]
		+ x[IX(M - 1, M - 1, M - 2)]);
	return x;
}

void FluidCube::addCollisionCube(int x, int y, int z, int size) {
	CollisionObj c;
	c.x = x;
	c.y = y;
	c.z = z;
	c.size = size;

	this->collisionCubes.push_back(c);

}

bool FluidCube::isCollision(int x, int y, int z) {
	for (int k = 0; k < this->collisionCubes.size(); k++) {
		CollisionObj obj = this->collisionCubes[k];
		if (x >= obj.x && x < obj.x + obj.size &&
			y >= obj.y && y < obj.y + obj.size &&
			z >= obj.z && z < obj.z + obj.size) {
			return true;
		}
	}
		return false;
}

float FluidCube::getDensity(int x, int y, int z) {
	return this->density[IX(x, y, z)];
}
