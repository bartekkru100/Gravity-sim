#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "winbgi2.h"

struct Vector
{
	/*////////////////////////////////////////////////////////////////
	Not to be confused with C++ vectors
	/*////////////////////////////////////////////////////////////////

	double x;
	double y;
	double z;
};
struct Body
{
	int index;
	double mass;
	Vector pos;
	Vector vel;
	double rad;
	bool has_collided;
	int collision_count;
};

const double G = 6.674301515E-11;
const double PI = atan(1) * 4.;
const double SOL_MASS = 2E30;
const double SOL_RAD = 7E8;
const double EARTH_MASS = 6E24;
const double EARTH_RAD = 6E6;
const double AU = 1.496E11;
const Vector ZERO_V = { 0,0,0 };

double betterRand();
double distance(Vector point1, Vector point2);
Vector crossProduct(Vector vector1, Vector vector2);
Vector randUnitVector();
Vector findUnitVector(Vector vector);
bool initialiseFromFile(Body** body, int* n, double* dt, double* scale);
bool initialiseRandomStarSystem(Body** body, int* n, double* dt, double* scale);
void eulerForce(Body* body, int n, double dt);
void symplecticEulerForce(Body* body, int n, double dt);
void rk4Force(Body* body, int n, double dt);
void verletForce(Body* body, int n, double dt);
void collisionDetection(Body** body, int* n, double dt);
void collide(Body** body, int* n, double dt, int i, int j, double ct);

int main()
{
	graphics(1600, 900);
	int i, j;
	int n = 0;
	Body* body = (Body*)malloc(sizeof(Body));
	bool is_init = false;
	char init_choice[256];
	double scale;
	double dt;

	/*////////////////////////////////////////////////////////////////
	Here we choose the initialisation method
	/*////////////////////////////////////////////////////////////////

	while (is_init == false)
	{
		printf("Choose initialisation method:\n1 - from file\n2 - random star system\n3 - exit\n");
		scanf("%s", &init_choice);
		{
			switch (atoi(init_choice))
			{
			case 1:
				if (initialiseFromFile(&body, &n, &dt, &scale) == true)is_init = true;
				break;
			case 2:
				if (initialiseRandomStarSystem(&body, &n, &dt, &scale) == true)is_init = true;
				break;
			case 3:
				return 0;
			default:
				puts("Invalid input!\nPlease try again:");
				break;
			}
		}
	}

	/*////////////////////////////////////////////////////////////////
	Setting up the graphics window
	/*////////////////////////////////////////////////////////////////

	Vector* trail[50];
	for (i = 0; i < 50; i++)trail[i] = (Vector*)malloc(n * sizeof(Vector));
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < 50; j++)
		{
			trail[j][i].x = body[i].pos.x;
			trail[j][i].y = body[i].pos.y;
			trail[j][i].z = body[i].pos.z;
		}
	}
	line(0, 20, 1600, 20);
	outtextxy(5, 5, "Z-View");
	outtextxy(805, 5, "Y-View");
	int old_n = n;
	int k = 99;
	int l = 0;
	int end_l = 100000;
	FILE* result = fopen("result.txt", "w");
	double init_total_energy = 0;
	double total_energy[10];
	double avg_total_energy = 0;
	for (i = 0; i < 10; i++)
	{
		total_energy[i] = 0;
	}
	for (i = 0; i < n; i++)
	{
		init_total_energy += body[i].mass * pow(distance(body[i].vel, ZERO_V), 2.) / 2.; 
		double temp_dist = 0;
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				init_total_energy -= G * body[i].mass * body[j].mass / temp_dist / 2.;
			}
		}
	}
	double error;
	while (animate(3000))
	{
		for (i = 9; i > 0; i--)
		{
			total_energy[i] = total_energy[i - 1];
		}
		total_energy[0] = 0;
		for (i = 0; i < n; i++)
		{
			total_energy[0] += body[i].mass * pow(distance(body[i].vel, ZERO_V), 2.) / 2.;
			double temp_dist = 0;
			for (j = 0; j < n; j++)
			{
				if (j == i)continue;
				else
				{
					temp_dist = distance(body[i].pos, body[j].pos);
					total_energy[0] -= G * body[i].mass * body[j].mass / temp_dist / 2.;
				}
			}
		}
		avg_total_energy = 0;
		for (i = 0; i < 10; i++)
		{
			avg_total_energy += total_energy[i];
		}
		avg_total_energy /= 10.;
		error = fabs(avg_total_energy / init_total_energy - 1);
		printf("error\t%lf\n", error);
		fprintf(result, "%d\t%lf\n", l, error);
		if (l >= end_l)break;
		/*////////////////////////////////////////////////////////////////
		By commenting or uncommenting four functions below
		you can switch between different methods of integration.
		Verlet integration is enabled by default, since it gives a nice balance
		between accuracy and ability to conserve energy and momentum.
		/*////////////////////////////////////////////////////////////////

		symplecticEulerForce(body, n, dt);
		//eulerForce(body, n, dt);
		//rk4Force(body, n, dt);
		//verletForce(body, n, dt);

		/*////////////////////////////////////////////////////////////////
		This function takes an array of bodies
		and checks for collisions in a given timestep.
		More on it later...
		/*////////////////////////////////////////////////////////////////

		collisionDetection(&body, &n, dt);

		/*////////////////////////////////////////////////////////////////
		Clear() doesn't clear the entire screen,
		so I had to improvise with using black lines instead.
		/*////////////////////////////////////////////////////////////////

		k++;
		if (k > 99)
		{
			setcolor(BLACK);
			for (i = 21; i < 900; i++)line(0, i, 1600, i);
			if (old_n != n)
			{
				for (i = 0; i < 50; i++)trail[i] = (Vector*)realloc(trail[i], n * sizeof(Vector));
				for (i = old_n; i < n; i++)
				{
					for (j = 0; j < 50; j++)
					{
						trail[j][i].x = body[i].pos.x;
						trail[j][i].y = body[i].pos.y;
						trail[j][i].z = body[i].pos.z;
					}
				}
			}
			old_n = n;
			for (i = 0; i < n; i++)
			{
				setcolor(DARKGRAY);
				for (j = 49; j > 0; j--)
				{
					trail[j][i].x = trail[j - 1][i].x;
					trail[j][i].y = trail[j - 1][i].y;
					trail[j][i].z = trail[j - 1][i].z;
				}
				trail[0][i].x = body[i].pos.x;
				trail[0][i].y = body[i].pos.y;
				trail[0][i].z = body[i].pos.z;
				for (j = 0; j < 49; j++)
				{
					if (400 + trail[j][i].x / scale < 800 && 450 + trail[j][i].y / scale > 21)
					{
						line(400 + trail[j][i].x / scale, 450 + trail[j][i].y / scale, 400 + trail[j + 1][i].x / scale, 450 + trail[j + 1][i].y / scale);
					}
					if (1200 + trail[j][i].x / scale > 800 && 450 + trail[j][i].z / scale > 21)
					{
						line(1200 + trail[j][i].x / scale, 450 + trail[j][i].z / scale, 1200 + trail[j + 1][i].x / scale, 450 + trail[j + 1][i].z / scale);
					}
				}

				setcolor(WHITE);
				line(800, 20, 800, 900);
				if (400 + body[i].pos.x / scale < 800 && 450 + body[i].pos.y / scale > 21)
				{
					circle(400 + body[i].pos.x / scale, 450 + body[i].pos.y / scale, ceil(10. * body[i].rad / scale));
				}
				if (1200 + body[i].pos.x / scale > 800 && 450 + body[i].pos.z / scale > 21)
				{
					circle(1200 + body[i].pos.x / scale, 450 + body[i].pos.z / scale, ceil(10. * body[i].rad / scale));
				}
			}
			k = 0;
		}
		l++;
	}
	free(body);
	for (i = 0; i < 50; i++)
	{
		free(trail[i]);
	}
	fclose(result);
	return 0;
}

double betterRand()
{
	/*////////////////////////////////////////////////////////////////
	Just a small function made to save me time
	when I need to generate a random number from between 0 and 1.
	/*////////////////////////////////////////////////////////////////

	return double(rand()) / double(RAND_MAX);
}
double distance(Vector point1, Vector point2)
{
	/*////////////////////////////////////////////////////////////////
	Another small function to save me some typing.
	There are a lot of instances where I need to calculate
	distances between two bodies or their relative speeds.
	/*////////////////////////////////////////////////////////////////

	return sqrt(pow(point1.x - point2.x, 2.) + pow(point1.y - point2.y, 2.) + pow(point1.z - point2.z, 2.));
}
Vector crossProduct(Vector vector1, Vector vector2)
{
	/*////////////////////////////////////////////////////////////////
	This function calculates a cross product of two vectors
	/*////////////////////////////////////////////////////////////////
	Vector product;
	product.x = vector1.y * vector2.z - vector1.z * vector2.y;
	product.y = vector1.z * vector2.x - vector1.x * vector2.z;
	product.z = vector1.x * vector2.y - vector1.y * vector2.x;
	return product;
}
Vector randUnitVector()
{

	/*////////////////////////////////////////////////////////////////
	This function generates a random unit vector
	/*////////////////////////////////////////////////////////////////
	Vector vector;
	double alpha = acos(1. - 2. * betterRand());
	double beta = PI * betterRand();
	vector.x = sin(alpha) * cos(beta);
	vector.y = sin(alpha) * sin(beta);
	vector.z = cos(alpha);
	return vector;
}
Vector findUnitVector(Vector vector)
{
	/*////////////////////////////////////////////////////////////////
	This takes a vector as an imput and finds its init vector
	/*////////////////////////////////////////////////////////////////
	Vector unit_vector = vector;
	unit_vector.x = vector.x / sqrt(pow(vector.x, 2.) + pow(vector.y, 2.) + pow(vector.z, 2.));
	unit_vector.y = vector.y / sqrt(pow(vector.x, 2.) + pow(vector.y, 2.) + pow(vector.z, 2.));
	unit_vector.z = vector.z / sqrt(pow(vector.x, 2.) + pow(vector.y, 2.) + pow(vector.z, 2.));
	return unit_vector;
}
bool initialiseFromFile(Body** body, int* n, double* dt, double* scale)
{

	/*////////////////////////////////////////////////////////////////
	This initialises a simulation from a file specified by the user
	/*////////////////////////////////////////////////////////////////
	char file_name[255];
	FILE* file = NULL;
	do
	{
		bool want_try = true;
		printf("Enter the file name:\n");
		scanf("%s", file_name);
		file = fopen(file_name, "r");
		if (file == NULL)
		{
			char response[255];
			printf("Couldn't find %s. Do you wish to try again?\nyes / no?\n");
			do
			{
				scanf("%s", response);
				if (strcmp(response, "no") == 0 || strcmp(response, "No") == 0)
				{
					want_try = false;
					break;
				}
				else if (strcmp(response, "yes") == 0 || strcmp(response, "Yes") == 0)
				{
					want_try = true;
					break;
				}
				else
				{
					printf("Invalid response. Type yes or no:\n");
				}
			} while (1);
		}
		if (want_try == false)return false;
	} while (file == NULL);
	char file_reader;
	*body = (Body*)realloc(*body, *n * sizeof(Body));
	int i = 0;
	fscanf(file, "step size:\t%lf\n", &(*dt));
	fscanf(file, "scale:\t\t%lf\n", &(*scale));
	while (!feof(file))
	{
		(*n)++;
		*body = (Body*)realloc(*body, *n * sizeof(Body));
		fscanf(file, "%c\n", &file_reader);
		fscanf(file, "mass\t%lf\n", &(*body)[i].mass);
		fscanf(file, "position:\nx\t%lf\n", &(*body)[i].pos.x);
		fscanf(file, "y\t%lf\n", &(*body)[i].pos.y);
		fscanf(file, "z\t%lf\n", &(*body)[i].pos.z);
		fscanf(file, "velocity:\nx\t%lf\n", &(*body)[i].vel .x);
		fscanf(file, "y\t%lf\n", &(*body)[i].vel.y);
		fscanf(file, "z\t%lf\n", &(*body)[i].vel.z);
		fscanf(file, "radius\t%lf\n", &(*body)[i].rad);
		(*body)[i].has_collided = false;
		(*body)[i].collision_count = 0;
		i++;
	}
	return true;
}
bool initialiseRandomStarSystem(Body** body, int* n, double* dt, double* scale)
{
	/*////////////////////////////////////////////////////////////////
	This function initialises the simulation with
	a procedurally generated star system
	/*////////////////////////////////////////////////////////////////

	srand(time(NULL));
	rand();
	int is_binary = rand() % 2;
	if (is_binary)*n = 3 + rand() % 7;
	else*n = 2 + rand() % 12;
	*body = (Body*)realloc(*body, *n * sizeof(Body));
	Vector pos_axis, vel_axis;
	double star_mass, star_rad, star_dist;
	double temp_rand;
	Vector temp_axis = randUnitVector();;
	temp_axis.z *= 100;
	temp_axis = findUnitVector(temp_axis);
	if (is_binary)
	{
		pos_axis = randUnitVector();
		pos_axis.z /= 10;
		pos_axis = findUnitVector(pos_axis);
		(*body)[0].index = 0;
		(*body)[1].index = 1;
		(*body)[0].mass = (pow(betterRand(), 4.) + 0.005) * 4E31;
		(*body)[1].mass = (pow(betterRand(), 4.) + 0.005) * 4E31;
		(*body)[0].rad = SOL_RAD* ((*body)[0].mass / SOL_MASS - pow((*body)[0].mass / SOL_MASS, 2.) / 80.)* (0.9 + 0.2 * betterRand());
		(*body)[1].rad = SOL_RAD * ((*body)[1].mass / SOL_MASS - pow((*body)[1].mass / SOL_MASS, 2.) / 80.) * (0.9 + 0.2 * betterRand());
		temp_rand = betterRand();
		(*body)[0].pos.x = (10. + 100. * temp_rand) * pos_axis.x * ((*body)[0].rad + (*body)[1].rad);
		(*body)[0].pos.y = (10. + 100. * temp_rand) * pos_axis.y * ((*body)[0].rad + (*body)[1].rad);
		(*body)[0].pos.z = (10. + 100. * temp_rand) * pos_axis.z * ((*body)[0].rad + (*body)[1].rad);
		(*body)[1].pos.x = -((*body)[0].mass * (*body)[0].pos.x) / (*body)[1].mass;
		(*body)[1].pos.y = -((*body)[0].mass * (*body)[0].pos.y) / (*body)[1].mass;
		(*body)[1].pos.z = -((*body)[0].mass * (*body)[0].pos.z) / (*body)[1].mass;
		star_mass = (*body)[0].mass + (*body)[1].mass;
		star_rad = (*body)[0].rad + (*body)[1].rad;
		star_dist = distance((*body)[0].pos, (*body)[1].pos);
		vel_axis = findUnitVector(crossProduct((*body)[0].pos, temp_axis));
		(*body)[0].vel.x = (0.9 + 0.2 * betterRand()) * vel_axis.x * sqrt(G * (*body)[1].mass * distance((*body)[0].pos, ZERO_V) / pow(star_dist, 2.));
		(*body)[0].vel.y = (0.9 + 0.2 * betterRand()) * vel_axis.y * sqrt(G * (*body)[1].mass * distance((*body)[0].pos, ZERO_V) / pow(star_dist, 2.));
		(*body)[0].vel.z = (0.9 + 0.2 * betterRand()) * vel_axis.z * sqrt(G * (*body)[1].mass * distance((*body)[0].pos, ZERO_V) / pow(star_dist, 2.));
		(*body)[1].vel.x = -(*body)[0].vel.x * (*body)[0].mass / (*body)[1].mass;
		(*body)[1].vel.y = -(*body)[0].vel.y * (*body)[0].mass / (*body)[1].mass;
		(*body)[1].vel.z = -(*body)[0].vel.z * (*body)[0].mass / (*body)[1].mass;
		(*body)[0].has_collided = false;
		(*body)[0].collision_count = 0;
		(*body)[1].has_collided = false;
		(*body)[1].collision_count = 0;
	}
	else
	{
		(*body)[0].index = 0;
		(*body)[0].mass = (pow(betterRand(), 4.) + 0.005) * 4E31;
		(*body)[0].rad = SOL_RAD * ((*body)[0].mass / SOL_MASS - pow((*body)[0].mass / SOL_MASS, 2.) / 80.) * (0.9 + 0.2 * betterRand());
		(*body)[0].pos.x = 0;
		(*body)[0].pos.y = 0;
		(*body)[0].pos.z = 0;
		star_mass = (*body)[0].mass;
		star_rad = (*body)[0].rad;
		star_dist = 0;
		(*body)[0].vel.x = 0;
		(*body)[0].vel.y = 0;
		(*body)[0].vel.z = 0;
		(*body)[0].has_collided = false;
		(*body)[0].collision_count = 0;
	}
	for (int i = 1 + is_binary; i < *n; i++)
	{
		(*body)[i].index = i;
		(*body)[i].mass = EARTH_MASS * (pow(4 * betterRand(), 5.) + 0.02);
		(*body)[i].rad = EARTH_RAD * cbrt((*body)[i].mass / EARTH_MASS / (0.7 + 0.4 * betterRand()) / (0.5 / (0.1 * (*body)[i].mass / EARTH_MASS + 0.5) + 0.2));
		pos_axis = randUnitVector();
		pos_axis.z /= (100-i*5);
		pos_axis = findUnitVector(pos_axis);
		temp_rand = betterRand();
		(*body)[i].pos.x = pos_axis.x * (0.95 + 0.1 * temp_rand) * pow(2.5 * (10. / *n), 2. / 3. * i) * (star_dist + star_rad / SOL_RAD * 0.2 * AU);
		(*body)[i].pos.y = pos_axis.y * (0.95 + 0.1 * temp_rand) * pow(2.5 * (10. / *n), 2. / 3. * i) * (star_dist + star_rad / SOL_RAD * 0.2 * AU);
		(*body)[i].pos.z = pos_axis.z * (0.95 + 0.1 * temp_rand) * pow(2.5 * (10. / *n), 2. / 3. * i) * (star_dist + star_rad / SOL_RAD * 0.2 * AU);
		temp_rand = betterRand();
		vel_axis = findUnitVector(crossProduct((*body)[i].pos, temp_axis));
		(*body)[i].vel.x = vel_axis.x * (0.9 + 0.2 * temp_rand) * sqrt(G * star_mass / distance((*body)[i].pos, ZERO_V));
		(*body)[i].vel.y = vel_axis.y * (0.9 + 0.2 * temp_rand) * sqrt(G * star_mass / distance((*body)[i].pos, ZERO_V));
		(*body)[i].vel.z = vel_axis.z * (0.9 + 0.2 * temp_rand) * sqrt(G * star_mass / distance((*body)[i].pos, ZERO_V));
		(*body)[i].has_collided = false;
		(*body)[i].collision_count = 0;
	}
	if (distance((*body)[0].vel, ZERO_V) > distance((*body)[1].vel, ZERO_V))
	{
		*dt = 0.1 * distance((*body)[0].pos, ZERO_V) / distance((*body)[0].vel, ZERO_V);
	}
	else
	{
		*dt = 0.1 * distance((*body)[1].pos, ZERO_V) / distance((*body)[1].vel, ZERO_V);
	}
	*scale = distance((*body)[*n-1].pos, ZERO_V)/500.;
	return true;
}
void eulerForce(Body* body, int n, double dt)
{
	/*////////////////////////////////////////////////////////////////
	This function integrates path using Euler method,
	for some weird reason it's actually harder to implement than
	semi-implicit Eurler method
	/*////////////////////////////////////////////////////////////////

	int i, j;
	double temp_dist;
	Vector* acc = (Vector*)malloc(n * sizeof(Vector));
	for (i = 0; i < n; i++)
	{
		acc[i].x = 0;
		acc[i].y = 0;
		acc[i].z = 0;
		temp_dist = 0;
		for (j = 0; j < n; j++)
		{
			temp_dist = distance(body[i].pos, body[j].pos);
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				acc[i].x += G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				acc[i].y += G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				acc[i].z += G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
	for (i = 0; i < n; i++)
	{
		body[i].pos.x += dt * body[i].vel.x;
		body[i].pos.y += dt * body[i].vel.y;
		body[i].pos.z += dt * body[i].vel.z;
	}
	for (i = 0; i < n; i++)
	{
		body[i].vel.x += dt * acc[i].x;
		body[i].vel.y += dt * acc[i].y;
		body[i].vel.z += dt * acc[i].z;
	}
	free(acc);
}
void symplecticEulerForce(Body* body, int n, double dt)
{
	/*////////////////////////////////////////////////////////////////
	This function integrates path using semi-implicit Euler method
	/*////////////////////////////////////////////////////////////////

	int i, j;
	double temp_dist;
	for (i = 0; i < n; i++)
	{
		temp_dist = 0;
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				body[i].vel.x += dt * G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				body[i].vel.y += dt * G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				body[i].vel.z += dt * G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
	for (i = 0; i < n; i++)
	{
		body[i].pos.x += dt * body[i].vel.x;
		body[i].pos.y += dt * body[i].vel.y;
		body[i].pos.z += dt * body[i].vel.z;
	}
}
void rk4Force(Body* body, int n, double dt)
{
	/*////////////////////////////////////////////////////////////////
	RK4 caused me a lot of headaches.
	The method itself isn't that hard to understand,
	but having navigate through a bunch of half a dozen for loops
	which all look almost exactly the same was painful,
	especially since I'm working in 3D here
	/*////////////////////////////////////////////////////////////////

	Vector* old_pos = (Vector*)malloc(n * sizeof(Vector));
	Vector* k1_vel = (Vector*)malloc(n * sizeof(Vector));
	Vector* k2_vel = (Vector*)malloc(n * sizeof(Vector));
	Vector* k3_vel = (Vector*)malloc(n * sizeof(Vector));
	Vector* k4_vel = (Vector*)malloc(n * sizeof(Vector));
	Vector* k1_acc = (Vector*)malloc(n * sizeof(Vector));
	Vector* k2_acc = (Vector*)malloc(n * sizeof(Vector));
	Vector* k3_acc = (Vector*)malloc(n * sizeof(Vector));
	Vector* k4_acc = (Vector*)malloc(n * sizeof(Vector));
	double h_dt = dt / 2.;
	double temp_dist;
	int i, j;
	//k1
	for (i = 0; i < n; i++)
	{
		old_pos[i].x = body[i].pos.x;
		old_pos[i].y = body[i].pos.y;
		old_pos[i].z = body[i].pos.z;
		k1_acc[i].x = 0;
		k1_acc[i].y = 0;
		k1_acc[i].z = 0;
		//vel
		k1_vel[i].x = body[i].vel.x;
		k1_vel[i].y = body[i].vel.y;
		k1_vel[i].z = body[i].vel.z;
		//acc
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				k1_acc[i].x += G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				k1_acc[i].y += G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				k1_acc[i].z += G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
	//k2
	for (i = 0; i < n; i++)
	{
		body[i].pos.x = old_pos[i].x + h_dt * k1_vel[i].x;
		body[i].pos.y = old_pos[i].y + h_dt * k1_vel[i].y;
		body[i].pos.z = old_pos[i].z + h_dt * k1_vel[i].z;
	}
	for (i = 0; i < n; i++)
	{
		k2_acc[i].x = 0;
		k2_acc[i].y = 0;
		k2_acc[i].z = 0;
		//vel
		k2_vel[i].x = k1_vel[i].x + h_dt * k1_acc[i].x;
		k2_vel[i].y = k1_vel[i].y + h_dt * k1_acc[i].y;
		k2_vel[i].z = k1_vel[i].z + h_dt * k1_acc[i].z;
		//acc
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				k2_acc[i].x += G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				k2_acc[i].y += G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				k2_acc[i].z += G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
	//k3
	for (i = 0; i < n; i++)
	{
		body[i].pos.x = old_pos[i].x + h_dt * k2_vel[i].x;
		body[i].pos.y = old_pos[i].y + h_dt * k2_vel[i].y;
		body[i].pos.z = old_pos[i].z + h_dt * k2_vel[i].z;
	}
	for (i = 0; i < n; i++)
	{
		k3_acc[i].x = 0;
		k3_acc[i].y = 0;
		k3_acc[i].z = 0;
		//vel
		k3_vel[i].x = k1_vel[i].x + h_dt * k2_acc[i].x;
		k3_vel[i].y = k1_vel[i].y + h_dt * k2_acc[i].y;
		k3_vel[i].z = k1_vel[i].z + h_dt * k2_acc[i].z;
		//acc
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				k3_acc[i].x += G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				k3_acc[i].y += G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				k3_acc[i].z += G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
	//k4
	for (i = 0; i < n; i++)
	{
		body[i].pos.x = old_pos[i].x + dt * k3_vel[i].x;
		body[i].pos.y = old_pos[i].y + dt * k3_vel[i].y;
		body[i].pos.z = old_pos[i].z + dt * k3_vel[i].z;
	}
	for (i = 0; i < n; i++)
	{
		k4_acc[i].x = 0;
		k4_acc[i].y = 0;
		k4_acc[i].z = 0;
		//vel
		k4_vel[i].x = k1_vel[i].x + dt * k3_acc[i].x;
		k4_vel[i].y = k1_vel[i].y + dt * k3_acc[i].y;
		k4_vel[i].z = k1_vel[i].z + dt * k3_acc[i].z;
		//acc
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				k4_acc[i].x += G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				k4_acc[i].y += G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				k4_acc[i].z += G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
	for (i = 0; i < n; i++)
	{
		body[i].pos.x = old_pos[i].x + dt * (k1_vel[i].x + 2. * k2_vel[i].x + 2. * k3_vel[i].x + k4_vel[i].x) / 6.;
		body[i].pos.y = old_pos[i].y + dt * (k1_vel[i].y + 2. * k2_vel[i].y + 2. * k3_vel[i].y + k4_vel[i].y) / 6.;
		body[i].pos.z = old_pos[i].z + dt * (k1_vel[i].z + 2. * k2_vel[i].z + 2. * k3_vel[i].z + k4_vel[i].z) / 6.;
		body[i].vel.x = k1_vel[i].x + dt * (k1_acc[i].x + 2. * k2_acc[i].x + 2. * k3_acc[i].x + k4_acc[i].x) / 6.;
		body[i].vel.y = k1_vel[i].y + dt * (k1_acc[i].y + 2. * k2_acc[i].y + 2. * k3_acc[i].y + k4_acc[i].y) / 6.;
		body[i].vel.z = k1_vel[i].z + dt * (k1_acc[i].z + 2. * k2_acc[i].z + 2. * k3_acc[i].z + k4_acc[i].z) / 6.;
	}
	free(old_pos);
	free(k1_vel);
	free(k2_vel);
	free(k3_vel);
	free(k4_vel);
	free(k1_acc);
	free(k2_acc);
	free(k3_acc);
	free(k4_acc);
}
void verletForce(Body* body, int n, double dt)
{
	/*////////////////////////////////////////////////////////////////
	After RK4 turned out to be mediocre I asked for help on Reddit on r/Physics
	where somebody pointed me towards velocity Verlet.
	This method is much simpler and causes almost no losses
	/*////////////////////////////////////////////////////////////////

	int i, j;
	double h_dt = dt / 2.;
	double temp_dist;
	for (i = 0; i < n; i++)
	{
		temp_dist = 0;
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				body[i].vel.x += h_dt * G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				body[i].vel.y += h_dt * G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				body[i].vel.z += h_dt * G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
	for (i = 0; i < n; i++)
	{
		body[i].pos.x += dt * body[i].vel.x;
		body[i].pos.y += dt * body[i].vel.y;
		body[i].pos.z += dt * body[i].vel.z;
	}
	for (i = 0; i < n; i++)
	{
		temp_dist = 0;
		for (j = 0; j < n; j++)
		{
			if (j == i)continue;
			else
			{
				temp_dist = distance(body[i].pos, body[j].pos);
				body[i].vel.x += h_dt * G * body[j].mass * (body[j].pos.x - body[i].pos.x) / pow(temp_dist, 3.);
				body[i].vel.y += h_dt * G * body[j].mass * (body[j].pos.y - body[i].pos.y) / pow(temp_dist, 3.);
				body[i].vel.z += h_dt * G * body[j].mass * (body[j].pos.z - body[i].pos.z) / pow(temp_dist, 3.);
			}
		}
	}
}
void collisionDetection(Body** body, int* n, double dt)
{
	/*////////////////////////////////////////////////////////////////
	This checks for collisions using two methods:
	first and the most straight forward is   
	just checking whether on not two bodies overlap,
	that's how most video game engines do it if I'm not mistaken,
	it's still not perfect though, at orbital speeds
	bodies will often move by several times their radius
	in each timestep, causing most would-be impacts to miss.
	For the second method I first came up with a formula
	for distance between bodies as a function of time.
	To make my life easier, I just assumed that the movement of bodies
	in between timesteps is roughly linear.
	After that I calculated its derivative to find
	the time of closest approach (local extremum of the function).
	At the end I made a copy of each body's position vector,
	moved it back in time by the time we got from looking for zeros
	in our derivative and checked whether or not the bodies overlap
	at that time.
	/*////////////////////////////////////////////////////////////////

	Vector* temp_pos = (Vector*)malloc(*n * sizeof(Vector));
	int i, j;
	int old_n = *n;
	for (i = 0; i < *n; i++)
	{
		temp_pos[i].x = (*body)[i].pos.x;
		temp_pos[i].y = (*body)[i].pos.y;
		temp_pos[i].z = (*body)[i].pos.z;
	}
	for (i = 0; i < *n; i++)
	{
		for (j = i+1; j < *n; j++)
		{
			double ct = ((*body)[j].pos.x + (*body)[j].pos.y + (*body)[j].pos.z - (*body)[i].pos.x - (*body)[i].pos.y - (*body)[i].pos.z) /
			((*body)[i].vel.x + (*body)[i].vel.y + (*body)[i].vel.z - (*body)[j].vel.x - (*body)[j].vel.y - (*body)[j].vel.z);
			if ((ct < 0. && -ct < dt) || distance((*body)[i].pos, (*body)[j].pos) <= (*body)[i].rad + (*body)[j].rad)
			{
				temp_pos[i].x += (*body)[i].vel.x * ct;
				temp_pos[i].y += (*body)[i].vel.y * ct;
				temp_pos[i].z += (*body)[i].vel.z * ct;
				temp_pos[j].x += (*body)[j].vel.x * ct;
				temp_pos[j].y += (*body)[j].vel.y * ct;
				temp_pos[j].z += (*body)[j].vel.z * ct;
				if ((distance(temp_pos[i], temp_pos[j]) <= (*body)[i].rad + (*body)[j].rad || distance((*body)[i].pos, (*body)[j].pos) <= (*body)[i].rad + (*body)[j].rad)
				&& (*body)[i].has_collided == false && (*body)[j].has_collided == false)
				{
					collide(&(*body), n, dt, i, j, ct);

					/*////////////////////////////////////////////////////////////////
					For some reason reallocating memory for
					temp_pos would give me random errors whenever
					the memory reallocated was smaller than
					initial memory from malloc().
					It took me hours to find the cause of the problem,
					but I've finally noticed that the crashes happen
					only after minor collisions that don't generate any debris.
					and therefore decrease the total number of bodies.
					All it took to fix it was adding an if statement to 
					realloc() to make sure it wouldn't decrease the allocated memory.
					/*////////////////////////////////////////////////////////////////

						temp_pos = (Vector*)realloc(temp_pos, *n * sizeof(Vector));
				}
				else
				{
					temp_pos[i].x = (*body)[i].pos.x;
					temp_pos[i].y = (*body)[i].pos.y;
					temp_pos[i].z = (*body)[i].pos.z;
					temp_pos[j].x = (*body)[j].pos.x;
					temp_pos[j].y = (*body)[j].pos.y;
					temp_pos[j].z = (*body)[j].pos.z;
				}
			}
		}
		(*body)[i].has_collided = false;
	}
	free(temp_pos);
}
void collide(Body** body, int* n, double dt, int i, int j, double ct)
{
	/*////////////////////////////////////////////////////////////////
	I consider this function to be the meat of this project.
	It takes everything I've learned during CS 1 classes
	and merges it into one. It gave me a better insight into
	how pointers work as they were to me (and most other students)
	basically like black magic. It causes bodies to marge on contact
	and spawn debris.
	/*////////////////////////////////////////////////////////////////

	/*////////////////////////////////////////////////////////////////
	First the bodies are moved to the point of collision.
	/*////////////////////////////////////////////////////////////////

	(*body)[i].pos.x += (*body)[i].vel.x * ct;
	(*body)[i].pos.y += (*body)[i].vel.y * ct;
	(*body)[i].pos.z += (*body)[i].vel.z * ct;
	(*body)[j].pos.x += (*body)[j].vel.x * ct;
	(*body)[j].pos.y += (*body)[j].vel.y * ct;
	(*body)[j].pos.z += (*body)[j].vel.z * ct;

	/*////////////////////////////////////////////////////////////////
	Momentum is calculated.
	/*////////////////////////////////////////////////////////////////

	Vector mom1, mom2, total_mom;
	double total_mass = (*body)[i].mass + (*body)[j].mass;
	double impact_vel = distance((*body)[i].vel, (*body)[j].vel);
	mom1.x = (*body)[i].vel.x * (*body)[i].mass;
	mom1.y = (*body)[i].vel.y * (*body)[i].mass;
	mom1.z = (*body)[i].vel.z * (*body)[i].mass;
	mom2.x = (*body)[j].vel.x * (*body)[j].mass;
	mom2.y = (*body)[j].vel.y * (*body)[j].mass;
	mom2.z = (*body)[j].vel.z * (*body)[j].mass;
	total_mom.x = mom1.x + mom2.x;
	total_mom.y = mom1.y + mom2.y;
	total_mom.z = mom1.z + mom2.z;
	double mass_ratio;
	if ((*body)[i].mass > (*body)[j].mass) mass_ratio = (*body)[j].mass / (*body)[i].mass;
	else mass_ratio = (*body)[i].mass / (*body)[j].mass;

	/*////////////////////////////////////////////////////////////////
	To make sure we don't end up with virtual Kessler syndrome,
	I've set up two little safeguards limiting debris creation.
	First one gives a limit on how small the smaller of the colliding bodies
	can be in order to be able to create debris on impact.
	Second one is an internal counter in each body that counts
	which generation of debris this body belongs to.
	If a body comes from a secondary collision (collision between debris)
	it won't generate any new bodies.
	/*////////////////////////////////////////////////////////////////

	bool is_major = false;
	if ((*body)[i].mass < 100 * (*body)[j].mass && (*body)[j].mass < 100 * (*body)[i].mass && (*body)[i].collision_count <= 0 && (*body)[j].collision_count <= 0)is_major = true;

	/*////////////////////////////////////////////////////////////////
	Colliding bodies merge according to laws of conservation of mass and momentum.
	/*////////////////////////////////////////////////////////////////

	(*body)[i].pos.x = ((*body)[i].mass * (*body)[i].pos.x + (*body)[j].mass * (*body)[j].pos.x) / total_mass;
	(*body)[i].pos.y = ((*body)[i].mass * (*body)[i].pos.y + (*body)[j].mass * (*body)[j].pos.y) / total_mass;
	(*body)[i].pos.z = ((*body)[i].mass * (*body)[i].pos.z + (*body)[j].mass * (*body)[j].pos.z) / total_mass;
	(*body)[i].mass = total_mass;
	(*body)[i].vel.x = total_mom.x / (*body)[i].mass;
	(*body)[i].vel.y = total_mom.y / (*body)[i].mass;
	(*body)[i].vel.z = total_mom.z / (*body)[i].mass;
	(*body)[i].rad = cbrt(pow((*body)[i].rad, 3.) + pow((*body)[j].rad, 3.));

	/*////////////////////////////////////////////////////////////////
	The hole after one of the consumed bodies is filled by the next body in an array.
	/*////////////////////////////////////////////////////////////////

	int k;
	for (k = j; k < *n - 1; k++)
	{
		(*body)[k].mass = (*body)[k + 1].mass;
		(*body)[k].pos.x = (*body)[k + 1].pos.x;
		(*body)[k].pos.y = (*body)[k + 1].pos.y;
		(*body)[k].pos.z = (*body)[k + 1].pos.z;
		(*body)[k].vel.x = (*body)[k + 1].vel.x;
		(*body)[k].vel.y = (*body)[k + 1].vel.y;
		(*body)[k].vel.z = (*body)[k + 1].vel.z;
		(*body)[k].rad = (*body)[k + 1].rad;
	}

	/*////////////////////////////////////////////////////////////////
	Realloc is used to free the last space in the array.
	/*////////////////////////////////////////////////////////////////

	(*n)--;
	*body = (Body*)realloc(*body, *n * sizeof(Body));

	/*////////////////////////////////////////////////////////////////
	Debris is generated below.
	/*////////////////////////////////////////////////////////////////

	if (is_major)
	{
		srand(time(NULL));
		Vector ejected_mom = ZERO_V;
		double alpha, beta;

		/*////////////////////////////////////////////////////////////////
		The number and mass of new bodies generated is based the impact velocity.
		/*////////////////////////////////////////////////////////////////

		int m = floor((10 + rand() % 20) * 2. * atan(impact_vel) / PI);
		*n += m;
		*body = (Body*)realloc(*body, *n * sizeof(Body));
		for (k = *n - m; k < *n; k++)
		{

			/*////////////////////////////////////////////////////////////////
			We randomly chose direction in which the fragments will fly away.
			/*////////////////////////////////////////////////////////////////

			alpha = acos(1. - 2. * betterRand());
			beta = 2. * PI * betterRand();
			(*body)[k].index = k;
			(*body)[k].mass = sqrt(betterRand()) * total_mass * atan(impact_vel/1000.) / PI / m * mass_ratio;
			(*body)[i].mass -= (*body)[k].mass;
			(*body)[k].vel.x = (*body)[i].vel.x + sin(alpha) * cos(beta) * impact_vel * sqrt(betterRand()) / cbrt(2.);
			(*body)[k].vel.y = (*body)[i].vel.y + sin(alpha) * sin(beta) * impact_vel * sqrt(betterRand()) / cbrt(2.);
			(*body)[k].vel.z = (*body)[i].vel.z + cos(alpha) * impact_vel * sqrt(betterRand()) / cbrt(2.);
			ejected_mom.x += (*body)[k].vel.x * (*body)[k].mass;
			ejected_mom.y += (*body)[k].vel.y * (*body)[k].mass;
			ejected_mom.z += (*body)[k].vel.z * (*body)[k].mass;
			(*body)[k].rad = cbrt((*body)[k].mass * pow((*body)[i].rad, 3.) / total_mass);

			/*////////////////////////////////////////////////////////////////
			Debris is generated at a small distance from the main body
			so they won't smash into it immidiately after being spawned.
			/*////////////////////////////////////////////////////////////////

			(*body)[k].pos.x = (*body)[i].pos.x + (*body)[k].vel.x * (dt + ct) + sin(alpha) * cos(beta) * ((*body)[k].rad + (*body)[i].rad);
			(*body)[k].pos.y = (*body)[i].pos.y + (*body)[k].vel.y * (dt + ct) + sin(alpha) * sin(beta) * ((*body)[k].rad + (*body)[i].rad);
			(*body)[k].pos.z = (*body)[i].pos.z + (*body)[k].vel.z * (dt + ct) + cos(alpha) * ((*body)[k].rad + (*body)[i].rad);
			(*body)[k].collision_count = (*body)[i].collision_count + 1;
			(*body)[k].has_collided = true;
		}

		/*////////////////////////////////////////////////////////////////
		We correct main body's velocity to satisfy law of conservation of momentum.
		/*////////////////////////////////////////////////////////////////

		(*body)[i].vel.x = (total_mom.x - ejected_mom.x) / (*body)[i].mass;
		(*body)[i].vel.y = (total_mom.y - ejected_mom.y) / (*body)[i].mass;
		(*body)[i].vel.z = (total_mom.z - ejected_mom.z) / (*body)[i].mass;
		(*body)[i].rad = cbrt(pow((*body)[i].rad, 3.) * (*body)[i].mass / total_mass);
	}

	/*////////////////////////////////////////////////////////////////
	We move the body forward in time because we moved it back at the beginning.
	/*////////////////////////////////////////////////////////////////

	(*body)[i].pos.x -= (*body)[i].vel.x * ct;
	(*body)[i].pos.y -= (*body)[i].vel.y * ct;
	(*body)[i].pos.z -= (*body)[i].vel.z * ct;
	(*body)[i].has_collided = true;
}