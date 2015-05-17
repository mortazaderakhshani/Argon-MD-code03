/* INCLUDE ESSENTIAL LIBRARIES */

#include <stdio.h>
#include <math.h>                                                     
#include <time.h>
#include <stdlib.h>
#include "stringlib.h"  

///////////////////////////////////////////////////////////////////////////////////////////////////////

/* DECLARATION OF VARIABLES */

void readConfigFile(char *, char *, char *, char *, double *, int *, int *, int *, double*, double *);

void initPositions(double **, int, double *);
void write_xyz(double **, int, int, double, FILE *);
void PrevCoord(double, double, double **, double **, double **);
void initVelocities(int, int, double, double, double **);
void write_Vel(int, int, double **, FILE *);
void computeForce_Energy(int, int, int, double, double, double, double, double *, double **, double **); 
void write_Force(int, int, double **, FILE *);
void CoordUpdate(int, double, double, double, double, double **, double **, double **, double **);
void Velupdate(int, int, double , double **, double **, double **, double *);

///////////////////////////////////////////////////////////////////////////////////////////////////////

/*STATIC VARIABLES */

const double R = 8.314;			 // Boltzman Constat times Avogadro # (Gass constant) J/(mol.k)
const double mass = 0.03995;		 // Mass of Ar or any Lennard-Jones Liquid,Units in kg/mol
const double eps = 1.6567;		 // Units of J/mol-------- or 0.210849 Units of kcal/mol, it is a meassure of strength
const double sigma = 3.345;		 // Sigma value for Ar or any Lennard-Jones Liquid, it is a measure of range of potential, Units of Angstrums

///////////////////////////////////////////////////////////////////////////////////////////////////////

/*MAIN CODE */

int main(){
	int n_atoms;				// Number of atoms
	double box;				// Size of the cubic box
	double length;				// Length of cubic box
	int n_iter;				// Number of MD iterations
	double temp;				// Temperature
	int delta_write;			// How often to write coordinates and log file in MD iterations
	double cutoff;				// cutoff distance (Angstrom)
	double cutoff_squ;			// nonbonding interaction cutoff distance squared
	double RT;				// Boltzman constant*Temperature Units of energy (J)
	double **atom_vel;			// Velocity array	
	double **coord;				// Particle coordinates array
	double **O_coord;			// Old coordinate array of Argon Particles
	int i,j,k;				// Generic indeces
	int seed=1;				// Random seed for velocity initialization
	double sigma6;				// LJL sigma^6 value
	double **forces_on_atom;		// Force acting on atom array fx, fy, fz
	//double Tff;				// add LJ potential energy 
	double Tpe;                  		// Total LJ Energy
	double dt;				// delta t value	
	double dt2;				// delta t value squared
	double im;				// Argon's inverse mass
	double Tke;				// Total Kinetic Energy
	double TE;				// Total Energy
	int iter=0;
	char log_FileName[1024];         	// Log file name
	char vel_FileName[1024];         	// Velocity file name
	char traj_FileName[1024];         	// Trajectoruy file name
	char force_FileName[1024];         	// Force file name



	FILE *logOut;
	FILE *xyzOut;
	FILE *velOut;
	FILE *forceOut;
	
	//printf("Tff is=%.3f\n", Tff);

/* READ CONFIG FILE FROM STANDARD IN */

        readConfigFile(log_FileName, vel_FileName, traj_FileName, force_FileName, &temp, &n_atoms, &n_iter, &delta_write, &cutoff, &dt);
	
	sigma6 = sigma*sigma*sigma*sigma*sigma*sigma;
	
	dt2=dt*dt;
	
	printf("dt2 is=%.2f\n\n", dt2);

	printf("sigma6 is=%.2f\n\n", sigma6);

	cutoff_squ = cutoff*cutoff;
	RT = R*temp;
	im = 1/mass;
	printf("R and T are =%.3f %.3f\n\n",R,temp);
	
	printf("RT and im are =%.3f %.3f\n\n",RT, im);

	printf("Cutoff is equal=%.3f\n\n",cutoff);

	// allocate coordinate  and velocity arrays
	coord = (double**) calloc(n_atoms,sizeof(double*));
	O_coord = (double**) calloc(n_atoms,sizeof(double*));
	// allocate velocity array memory
	atom_vel = (double**) calloc(n_atoms,sizeof(double*));   
	
	// allocate force array memory
	forces_on_atom = (double**) calloc(n_atoms,sizeof(double*));        
	  
	for (i=0;i<n_atoms;i++) {
		coord[i] = (double*) calloc(3,sizeof(double));
		O_coord[i] = (double*) calloc(3,sizeof(double));	
		atom_vel[i] = (double*) calloc(3,sizeof(double));	// Does not complie wiht calloc???
		forces_on_atom[i] = (double*) calloc(3,sizeof(double));
	}
	// open log files
	logOut = fopen(log_FileName,"w");
	xyzOut = fopen(traj_FileName,"w");
	velOut = fopen(vel_FileName,"w");
	forceOut = fopen(force_FileName, "w");

	// inintialize positions and velocities + compute forces//

	           			
	initPositions(coord, n_atoms, &box);
	write_xyz(coord, n_atoms, n_iter, box, xyzOut);
	

	initVelocities(n_atoms, seed, mass, RT, atom_vel);
	write_Vel(n_atoms, n_iter, atom_vel, velOut);
	
	Tke=0.5*mass*R*temp;
	printf("KE=%16.6E\n", Tke);
	computeForce_Energy(n_atoms, iter, delta_write, box, cutoff_squ, eps, sigma6, &Tpe, coord, forces_on_atom);	
	write_Force(n_atoms, n_iter, forces_on_atom, forceOut);
	
	PrevCoord(n_atoms, dt, coord, atom_vel, O_coord);

	fflush(xyzOut);
	fflush(velOut);
	fflush(forceOut);
	

////////////////////////////////////////////////////

/* Run MD itterations  */
// Use Verlet integration
	
/*	Tke=0.0;*/
/*	Tpe=0.0;*/
	
	for(iter=0;iter<n_iter;iter++) {
		
				

		CoordUpdate(n_atoms, im, dt, dt2, box, coord, atom_vel, forces_on_atom, O_coord);
		

		computeForce_Energy(n_atoms, iter, delta_write, box, cutoff_squ, eps, sigma6, &Tpe, coord, forces_on_atom);
 
		
		Velupdate(iter, delta_write, dt, coord, atom_vel, O_coord, &Tke);

		if(iter%delta_write==0) {

		printf("Tke=%12.6E\n", Tke);
/*TE=Tpe+Tke;*/
		
		//printf("Iteration=%d\t", iter);
		//printf("cutoff_squ is=%.3f\n", cutoff_squ);
			
			write_xyz(coord, n_atoms, iter, box, xyzOut);
			
			
			write_Vel(n_atoms, n_iter, atom_vel, velOut);

			write_Force(n_atoms, iter, forces_on_atom, forceOut);

			fflush(xyzOut);
			
			fflush(forceOut);
			
			fflush(velOut);

		}	

	}



	fclose(xyzOut);
	fclose(forceOut);
	fclose(velOut );
}

	



///////////////////////////////////////////////////////////////////////////////////////////////////////

/* SUBROUTINS AND FUNCTIONS */

/* READ CONFIG FILE */
	
void readConfigFile(char *log_FileName, char *vel_FileName, char *traj_FileName, char *force_FileName, double *temp, int *n_atoms, int *n_iter, int *delta_write, double *cutoff, double *dt)
	{

	char buffer[1024];
        char tempBuffer[1024];
        char check[15];
        char *firstWord;	
	
	while (fgets(buffer,1024,stdin) != NULL) {
		strncpy(tempBuffer,buffer,1024);
		firstWord=string_firstword(tempBuffer);
                if (strncmp(firstWord,"log_File",8)==0) {
		       strcpy(log_FileName,string_secondword(buffer));
               } else if (strncmp(firstWord,"vel_File",8)==0) {
		       strcpy(vel_FileName,string_secondword(buffer));
	       } else if (strncmp(firstWord,"traj_File",9)==0) {
                       strcpy(traj_FileName,string_secondword(buffer));
	       } else if (strncmp(firstWord,"force_File",10)==0) {
		       strcpy(force_FileName,string_secondword(buffer));
	       } else if (strncmp(firstWord,"temperature",11)==0) {
	               *temp = atof(string_secondword(buffer));
	       } else if (strncmp(firstWord,"n_atoms",7)==0) {
	               *n_atoms = atoi(string_secondword(buffer));
	       } else if (strncmp(firstWord,"n_iter",6)==0) {
	               *n_iter = atoi(string_secondword(buffer));
	       } else if (strncmp(firstWord,"delta_write",11)==0) {
	               *delta_write = atoi(string_secondword(buffer));
	       } else if (strncmp(firstWord,"cutoff",6)==0) {
		       *cutoff = atof(string_secondword(buffer));
	       } else if (strncmp(firstWord,"dt",2)==0) {
		       *dt = atof(string_secondword(buffer));
	       }
	}
}

		

		

////////////////////////////////////////////////////


/* Initilazie Positions */

void initPositions(double **coord, int n_atoms, double *box) {

	double cbrt(double x);          // cube root function
	int iBoxD;                      // integer box dimension
	double fBoxD;                   // float box dimension
	
	int x, y, z;
	double xPos, yPos, zPos;
	int atomCount;

	// determine how many bins to divide the box into
	iBoxD = (int) cbrt((double) n_atoms);	// number of atoms should be an integer to the third power
	if (iBoxD*iBoxD*iBoxD < n_atoms) {
		iBoxD++;
	}

	// determine the size of the bins
	fBoxD = 3.55;        	// seems unitless since multipleied by (x+0.5) gives xPos but I think it gives dimention to particles 
	*box = iBoxD*fBoxD;     // calculate dimension of the box

	// add an atom to each created bin 
	atomCount = 0;
	for(x=0;x<iBoxD;x++) {
		xPos = (x+0.5)*fBoxD;
		for(y=0;y<iBoxD;y++) {
			yPos = (y+0.5)*fBoxD;
			for(z=0; z<iBoxD; z++) {
				if (atomCount < n_atoms) {
					zPos = (z+0.5)*fBoxD;
					coord[atomCount][0]=xPos;
					coord[atomCount][1]=yPos;
					coord[atomCount][2]=zPos;
					atomCount++;
				} else {
					break;
				}
			}

			if (atomCount>=n_atoms) {
				break;
			}
		}
	}
}

/////////////////////////////////////////////////////////

// Compute pervious step coordinates
	
void PrevCoord(double n_atoms, double dt, double **coord, double **atom_vel, double **O_coord) {

	int i, j;
	
	for(i=0; i<n_atoms; i++) {
		for(j=0; j<3; j++) {
		O_coord[i][j] = coord[i][j] - atom_vel[i][j]*dt;	

		}
	}
}	
//////////////////////////////////////////////////////////////

/* Initialize Velocities */

void initVelocities(int n_atoms, int seed, double mass, double RT, double **atom_vel) {
	 
	double sumv[3] = {0.0, 0.0, 0.0};
	double msv = 0;
	int i, j;
	double scale_factor;
	double total_number_of_element = n_atoms*3;

	double sqrt(double x);

	srand((unsigned) seed);

	for(i=0; i<n_atoms; i++) {
		for(j=0; j<3; j++) {		
			// randomly assigned velocities in th range of -0.5 - 0.5
			atom_vel[i][j] = (rand()/((double) RAND_MAX))-0.5;
			//printf("%f  ", atom_vel[i][j]);
			// sum in each coordinate
			sumv[j] = sumv[j] + atom_vel[i][j]/n_atoms;
			// mean square velocity
			msv = msv + (atom_vel[i][j]*atom_vel[i][j])/total_number_of_element; //unit of m2/S2
		}
                //printf("\n");
	}
	scale_factor = (1E-2)*sqrt(3.0*RT/(mass*msv)); 		//Scale factor is unitless -- 1E-2 is to convert m2/S2 to A2/fs2
	// printf("scale factor %f mass %f", RT, mass);
	for(i=0; i<n_atoms; i++) {
		for(j=0; j<3; j++) {
			atom_vel[i][j] = (atom_vel[i][j] - sumv[j])*scale_factor;
			//printf("%f  ", atom_vel[i][j]);
		}

	}

}
////////////////////////////////////////////////////

// Sub: Compute Energy and forces

void computeForce_Energy(int n_atoms, int iter, int delta_write, double box, double cutoff_squ, double eps, double sigma6, double *Tpe, double **coord, double **forces_on_atom) {
	
	
	int atom1;
        int atom2;
        double r2;
        double r2rev;
        double r6rev;
        double rev6sig6;

        int j;
	

        double element[3] = {0.0, 0.0, 0.0};

	double ff;
	//double Tff;

	//Tff = 0.0;
	*Tpe = 0.0;

	for(atom1=0; atom1<n_atoms; atom1++) {
		for(j=0;j<3;j++) {
			forces_on_atom[atom1][j] = 0.0;
			forces_on_atom[atom2][j] = 0.0;
		}
	}
	

	for(atom1=0; atom1<n_atoms; atom1++){

		for(atom2=atom1+1;atom2<n_atoms;atom2++){
			r2=0.0;
			
			for(j=0;j<3;j++) {

				element[j]= coord[atom2][j]-coord[atom1][j];
				if(element[j]<-box/2.0){
					element[j]=element[j]+box;
				} else if (element[j]>box/2.0) {
					element[j]=element[j]-box;
				}
				r2+=element[j]*element[j];

			}
				//printf("r2=%f.3\n", r2);
		
			if(r2<cutoff_squ) {
				r2rev=1.0/r2;
				//printf("r2rev=%f.3\n", r2rev);
				r6rev=r2rev*r2rev*r2rev;
				//printf("r6rev=%f.3\n", r6rev);
				rev6sig6=sigma6*r6rev;
				//printf("rev6sig6=%f.3\n", rev6sig6);

				ff=48*eps*r2rev*rev6sig6*(rev6sig6-0.5);
				//printf("ff=%f.3\n", ff);
				for(j=0; j<3; j++) {
					forces_on_atom[atom1][j]-= ff*element[j];	  // forces acting on atoms
					forces_on_atom[atom2][j]+= ff*element[j];
				}
				if(iter%delta_write==0) {
					*Tpe += 4.0*eps*rev6sig6*(rev6sig6 - 1.0);	// Total potential Energy in each deltawrite
					printf("Tpe=%16.6E\n", *Tpe);
				}


			} 

		}
	}
	

}



////////////////////////////////////////////////////

//MD LOOP SUBROUTINS


// compute new coordinates

void CoordUpdate(int n_atoms, double im, double dt, double dt2, double box, double **coord, double **atom_vel, double **forces_on_atom, double **O_coord) {
	
	int i, j;
	double element;
	
	for (i=0; i<n_atoms; i++) {
		for (j=0; j<3; j++) {
			element = 2*coord[i][j] - O_coord[i][j] + im*dt2*forces_on_atom[i][j];
			// wrapping molecules within box
			if(element<0) {
				element += box;
			} else if(element > box) {
				element -= box;
			}
			
			O_coord[i][j] = coord[i][j];			
			coord[i][j] = element;
			
			//

			element = 0.0;
		}
	}
}
////////////////////////////////////////////////////




// compute new velocities


void Velupdate(int iter, int delta_write, double dt, double **coord, double **atom_vel, double **O_coord, double *Tke) {
	
	int i, j;
	int n_atoms;
	double mvs=0.0;
	*Tke=0.0;

	for (i=0; i<n_atoms; i++) {
		for (j=0; j<3; j++) {
			atom_vel[i][j] = (coord[i][j] - O_coord[i][j])/dt;
			}
			if(iter%delta_write==0) {
			mvs +=(atom_vel[i][j]*atom_vel[i][j]);
			
		}

	}
			
/*			printf("mvs=%f\n", mvs);*/
/*			*Tke=0.5*mass*mvs;*/
/*			printf("Tke=%12.6E\n", *Tke);*/
}

////////////////////////////////////////////////////

// Sub: Write Output data

// WRITE 3D COORDINATES

void write_xyz(double **coord, int n_atoms, int n_iter, double box, FILE *xyzOut) {

	int i,j,k;
	int atom;
	
	fprintf(xyzOut, "%d \n", n_atoms);
	fprintf(xyzOut, "Step %d , box: %9.6f  %9.6f  %9.6f\n", n_iter, box, box, box);
	for (atom=0;atom<n_atoms;atom++) {
		fprintf(xyzOut, "Ar%12.6f%12.6f%12.6f\n",coord[atom][0],coord[atom][1],coord[atom][2]);

	}
}



// Write Velocities 


void write_Vel(int n_atoms, int n_iter, double **atom_vel, FILE *velOut) {
	int atom;

	fprintf(velOut, "Step %d\n", n_iter);
	for(atom=0; atom<n_atoms; atom++) {
		fprintf(velOut, "Ar%12.6f%12.6f%12.6f\n", atom_vel[atom][0], atom_vel[atom][1], atom_vel[atom][2]);
	}
}

// Write Forces 


void write_Force(int n_atoms, int n_iter, double **forces_on_atom, FILE *forceOut) {

	int atom;

	fprintf(forceOut, "Step %d\n", n_iter);
	for(atom=0; atom<n_atoms; atom++) {
		fprintf(forceOut, "Ar%16.6E%16.6E%16.6E\n", forces_on_atom[atom][0], forces_on_atom[atom][1], forces_on_atom[atom][2]);
	}
}

////////////////////////////////////////////////////
//
//
//
//
//
//
