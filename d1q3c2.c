/****************************************************************/
/*        Isothermal D1Q3 Lattice Boltzmann Simulation          */
/*                    Two Component System                      */
/*            Kyle Strand:  kyle.t.strand@ndsu.edu              */
/*               North Dakota State University                  */
/*                    25 November 2015                          */
/****************************************************************/

/* There are currently two initialization routines. init() is the 
   default init. This definies sine density profiles with null 
   initial velocities of componenets. init2() is a constant density 
   model which requires initial velcoities. All values can be edited 
   through the interface. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mygraph.h>
#include <unistd.h>
#include <time.h>

#define xdim 5000                                   // Number of x lattice points
#define C 2

double f0[C][xdim], f1[C][xdim], f2[C][xdim], n[C][xdim], F0[C][xdim], F1[C][xdim], F2[C][xdim], u[C][xdim];
double u0[C]={0.1,-0.1}, th[C], uprev[C];
double omega = 1, n0[C]={0.6,0.4}, Amp=0.001, T0 = 1./3., Fric=0.001;                  //Fric is the scalar quantity attached to the F_x variable
//double ntot=n0[0]+n0[1], ptot=n0[0]*u0[0]+n0[1]*u0[1];
double ug[C][xdim];
int ugreq[C] = {0,0};
int next = 0, Pause = 1, done = 0, repeat = 1, iterations, writedata = 0, interval=100,printtheory=0;
char filename[500]="finite.txt";

void DataToFile() {                                 // Routine to output data
  char path[500]="Data/";
  strcat(path,filename);
  FILE *data = fopen(path,"a");
  if (printtheory == 1) fprintf(data,"%f %f %f %f %i\n",u[0][0],u[1][0],th[0],th[1],iterations);
  if (printtheory == 0) fprintf(data,"%f %f %i\n",u[0][0],u[1][0],iterations);
  fclose(data);
}

void theory() {                                      // This will output theoretical values if desired.
  double ntot=n[0][0]+n[1][0], ptot=n[0][0]*u[0][0]+n[1][0]*u[1][0];
  double ubar=(n[0][0]*u0[0]+n[1][0]*u0[1])/(n[0][0]+n[1][0]);
  //for (int i=0; i<C; i++) th[i]=ubar+(u0[i]-ubar)*exp(-Fric*(n[0][0]+n[1][0])*iterations); //Exact solution
  /*if (iterations==0) {
    th[0]=uprev[0];
    th[1]=uprev[1];
  } else {
    th[0]=Fric*n[1][0]*(u[1][0]-u[0][0])+uprev[0];  //Finite difference solutions
    th[1]=Fric*n[0][0]*(u[0][0]-u[1][0])+uprev[1];
  }*/
  th[0]=(n[1][0]*u[1][0])/n[0][0]+(ptot*(n[0][0]-n[1][0])-ptot*(n[0][0]-n[1][0])*pow(1-Fric*ntot,iterations)+(n0[0]*u0[0]-n0[1]*u0[1])*pow(1-Fric*ntot,iterations))/(ntot*n[0][0]);
  th[1]=(n[0][0]*u[0][0])/n[1][0]+(ptot*(n[1][0]-n[0][0])-ptot*(n[1][0]-n[0][0])*pow(1-Fric*ntot,iterations)+(n0[1]*u0[1]-n0[0]*u0[0])*pow(1-Fric*ntot,iterations))/(ntot*n[1][0]);
}

void init2() {                                       // Initializing Eq. Dists

  iterations = 0;

  for (int c = 0; c < C; c++) {
    for (int i = 0; i < xdim; i++) {                  // not sure if these are correct
      n[c][i] = n0[c];     //try equal densities out of phase
      u[c][i] = u0[c];
      //n[i] = (i<xdim/3)? n0+Amp:n0;                 // step function
      f0[c][i] = n[c][i] * (1 - T0);
      f1[c][i] = n[c][i]/2 * T0 - n[c][i]*u[c][i]/2;
      f2[c][i] = n[c][i]/2 * T0 + n[c][i]*u[c][i]/2;
      //u[c][i] = (f2[c][i] - f1[c][i])/n[c][i];
    }
  }
}


void init() {                                       // Initializing Eq. Dists

  iterations = 0;

  for (int c = 0; c < C; c++) {
    for (int i = 0; i < xdim; i++) {                  
      n[c][i] = n0[c] + Amp * sin(2*M_PI*i/xdim);     //try equal densities out of phase
      //n[i] = (i<xdim/3)? n0+Amp:n0;                 // step function
      f0[c][i] = n[c][i] * (1 - T0);
      f1[c][i] = n[c][i]/2 * T0;
      f2[c][i] = n[c][i]/2 * T0;
      u[c][i] = (f2[c][i] - f1[c][i])/n[c][i];
    }
  }
}

void iteration() {                                  // Iteration step
  
  double tmp, ax, F;

  uprev[0]=u[0][0];
  uprev[1]=u[1][0];

  //iterations++;

  for (int c = 0; c < C; c++) {                                     //all u's must be calculated first
      for (int i = 0; i < xdim; i++) {
        n[c][i] = f0[c][i] + f1[c][i] + f2[c][i];                   // Sum of Eq dists = density
        u[c][i] = (f2[c][i] - f1[c][i])/n[c][i];
      }
  }
  for (int c=0; c<C; c++) {

    for (int i=0; i<xdim; i++) {
        ax = 0;
        for (int d = 0; d < C; d++) ax+= n[d][i]*(u[d][i]-u[c][i]);
        F = Fric*ax;
        F0[c][i] = -n[c][i]*F*(2*u[c][i]+F);
        F1[c][i] = 0.5*n[c][i]*F*(F+2*u[c][i]-1);                         // F_-1
        F2[c][i] = 0.5*n[c][i]*F*(F+2*u[c][i]+1);                         // F_1
        f0[c][i] += omega*(n[c][i]*(1-T0-u[c][i]*u[c][i])-f0[c][i]) + F0[c][i];                 // f_0
        f1[c][i] += omega*(n[c][i]/2*(-u[c][i]+T0+u[c][i]*u[c][i])-f1[c][i]) + F1[c][i];      // f_-1
        f2[c][i] += omega*(n[c][i]/2*(u[c][i]+T0+u[c][i]*u[c][i])-f2[c][i]) + F2[c][i];       // f_1
    }
  }

  for (int c = 0; c < C; c++) {
    tmp = f1[c][0];
    memmove(&f1[c][0], &f1[c][1], (xdim-1)*sizeof(double));
    f1[c][xdim-1] = tmp;
    tmp=f2[c][xdim-1];
    memmove(&f2[c][1], &f2[c][0], (xdim-1)*sizeof(double));
    f2[c][0] = tmp;
  }

  if (printtheory == 1) theory();

  if (iterations % interval == 0) DataToFile();

  iterations++;
}

void GUI() {
  static int Xdim = xdim;

  DefineGraphN_R("n0",&n[0][0],&Xdim,NULL);
  DefineGraphN_R("n1",&n[1][0],&Xdim,NULL);
  DefineGraphN_R("u0",&ug[0][0],&Xdim,&ugreq[0]);
  DefineGraphN_R("u1",&ug[1][0],&Xdim,&ugreq[1]);
  StartMenu("D1Q3",1);
  DefineInt("Iterations",&iterations);
  DefineDouble("Omega",&omega);
  DefineDouble("Amp",&Amp);
  DefineDouble("Fric",&Fric);
  DefineDouble("n0",&n0[0]);
  DefineDouble("n1",&n0[1]);
  DefineDouble("u0",&u0[0]);
  DefineDouble("u1",&u0[1]);
  DefineDouble("T0",&T0);
  DefineFunction("init",&init);
  DefineFunction("init2",&init2);
  DefineGraph(curve2d_,"Graphs");
  DefineInt("Repeat",&repeat);
  DefineBool("Next",&next);
  DefineBool("Pause",&Pause);
  StartMenu("Data Output",0);
  DefineString("Name: ",filename,200);
  DefineInt("Output Interval",&interval);
  DefineBool("Theory Data",&printtheory);
  DefineBool("Write Data",&writedata);
  EndMenu();
  DefineBool("Done",&done);
  EndMenu();
}

void GetData() {
  int i;

  for (int c = 0; c < C; c++) {
    if (ugreq[c]) {
      for (i = 0; i < xdim; i++) ug[c][i] = (f2[c][i] - f1[c][i])/n[c][i];
      ugreq[c] = 0;
    }
  }
}

int main(int argc, char *argv[]) {
  int newdata = 1;
  int i;

  init();
  GUI();

  while (done == 0) {
    Events(newdata);
    GetData();
    DrawGraphs();
    if (next || !Pause) {
      newdata = 1;
      next = 0;
      for (i = 0; i < repeat; i++) {
        iteration();
      }
    }
    else sleep(1);
  }

  return 0;
}
