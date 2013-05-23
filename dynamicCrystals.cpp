#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <time.h>
#include <stdlib.h>

#define PI 3.14159265

using namespace std;

//struct for colloid
struct Colloid{
	double 	x,y;		//Dimensions in pixels
	int 	r;
	double 	theta;
	double	v;
	bool 	ctc;		//in contact of not
};

//Parameters
int 	SC	= 20;			//Surface Coverage
int 	L	= 500;			//500x500 px box
int 	radius  = 10;			//radius in px
int 	tp 	= 15;			//Persistant Time
int 	xtp	= 50;			//Persistant Time Cycles
double 	steps	= 200.;			//Steps per persistant time
double 	deltaT	= 1/steps;		
double  v0	= 1.0;

int duration	= (int)(xtp*tp*steps);	//Number of iterations in simulation
double threshold= 3.0;			//Distance for attraction
double	coeff	= 0.87/2.0;


//Declare Prototypes
void initialize(Colloid* p, int N, int L);
void propagate(Colloid* p, int N);
void correct_overlap(Colloid &A, Colloid &B, double dist, double ang);
void pbc(Colloid* p, int N, int L);
double distance(int x1, int x2, int y1, int y2);
bool check_collision(Colloid A, Colloid B);
int numColloids(double rho, int L, int r);

int main(){
	/*
	 *	Set up simulation
	 */
	int N = numColloids(10.,10,1);
	Colloid *particles = new Colloid[N];
	//Initialize random starting points for particles
	initialize(particles, N, L);
	/*
	 *	Begin the simulation
	 */
	for(int i=0; i<duration; ++i){
		//Move particles in some direction
		propagate(particles, N);
		//Come up with a smarter way to do this, so it uses less 
		//memory...too lazy now o.O
		pbc(particles, N, L);
		//Record every quarter of a tp
		if(i%(int)(0.25*steps)==0){
			cout << "Print Data" << endl;	
		}
	}

	//Free up memory because we are nice people
	delete [] particles;
	return 0;
}

/*
 *	Initiate starting positions of particles
 */
void initialize(Colloid* p, int N, int L){
	srand (time(NULL));	//initialize random seed
	for(int i=0; i<N; ++i){
		p[i].x 		= rand()%L+1.0;
		p[i].y 		= rand()%L+1.0;
		p[i].theta	= rand()%2*PI+0.5;
		p[i].v	 	= v0;
		p[i].r 		= radius;
		p[i].ctc	= false;
	}
}

/*
 *	Move particles one step - includes P.B.C.s
 */
void propagate(Colloid* p, int N){
	srand(time(NULL));	//initialize random seed
	for(int i=0; i<N; ++i){
		//Move particle
		p[i].x = p[i].x+deltaT*cos(p[i].theta)*p[i].v;
		p[i].y = p[i].y+deltaT*sin(p[i].theta)*p[i].v;
		for(int j=0; j<N; ++j){
			if(i==j)	//skip redundancy
				continue;
			double dist = distance(p[i].x,p[j].x,p[i].y,p[j].y);
			//Apply attractive force
			if(dist<threshold){
				double ang	= 0.0;
				double tan_ang 	= (p[j].y-p[i].y)/(p[j].x-p[i].x);
				if(p[j].x-p[i].x>0)
					ang	= atan(tan_ang);
				else
					ang	= atan(tan_ang)+PI;
				//Move i and j closer together
				p[j].x = p[j].x-deltaT*coeff/(pow(dist,2))*cos(ang);
				p[j].y = p[j].y-deltaT*coeff/(pow(dist,2))*sin(ang);
				p[i].x = p[i].x+deltaT*coeff/(pow(dist,2))*cos(ang);
				p[i].y = p[i].x+deltaT*coeff/(pow(dist,2))*sin(ang);
			//}	//CHECK THIS TO MAKE SURE IT IS SMART
			//Update the distance between particles
			dist = distance(p[i].x,p[j].x,p[i].y,p[j].y);
			//Correct overlap
			if(check_collision(p[i],p[j]))
				correct_overlap(p[i],p[j],dist,ang);
			}	
		}
		//Apply smearing to the direction
		p[i].theta=p[i].theta+sqrt(2/((steps)*tp))*rand();

	}
}
/*
 *	Correct for overlaps
 */
void correct_overlap(Colloid &A, Colloid &B, double dist, double ang){
	A.x = A.x-(1-dist)*cos(ang)/2.;
	A.y = A.y-(1-dist)*sin(ang)/2.;
	B.x = B.x+(1-dist)*cos(ang)/2.;
	B.y = B.y+(1-dist)*sin(ang)/2.;
}
/*
 *	Periodic Boundary Conditions
 */
void pbc(Colloid* p, int N, int L){
	for(int i=0; i<N; ++i){
		if(p[i].x>=L)
			p[i].x=p[i].x-L;
		if(p[i].y>=L)
			p[i].y=p[i].y-L;
		if(p[i].x<=0)
			p[i].x=p[i].x+L;
		if(p[i].y<=0)
			p[i].y=p[i].y+L;
	}
}

/*
 *	Calculate distance between two points
 */
double distance(int x1, int x2, int y1, int y2){
	return sqrt(pow(x2-x1,2)+pow(y2-y1,2));
}

/*
 *	Check for collisions
 */
bool check_collision(Colloid A, Colloid B){
	if(distance(A.x,B.x,A.y,B.y)<(A.r+B.r)){
		return true;
	}else
		return false;
}
/*
 *	Determine number of colloids needed
 */
int numColloids(double rho, int L, int r){
	return (int)(1.0*rho*L*L)/(PI*pow(r,2));
}
