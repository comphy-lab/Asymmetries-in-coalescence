/* Title: Coalescence of bubbles

# Version 1.0
# Last modified: Oct 15, 2024

# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

// 1 is bubble
#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

#if !MPI
#include "distance.h"
#endif

int MAXlevel; // command line input

#define tsnap (1e-2)
#define tsnap2 (1e-4)
// Error tolerances
#define fErr (1e-3)                                 // error tolerance in VOF
#define VelErr (1e-2)                            // error tolerances in velocity
#define TOL (1e-2)                                 // error tolerance in position

// boundary conditions
f[left] = dirichlet(0.0);
u.t[left] = dirichlet(0.0);

// Other command - line inputs.
double tmax, MuRin, OhOut, RhoIn;
double Rr, zWall;
double Ldomain; // Dimension of the bugger drop and the domain
char nameOut[80], dumpFile[80];

int main(int argc, char const *argv[]) {
  if (argc < 7){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 7-argc);
    return 1;
  }

  // Values taken from the terminal
  MuRin = 1e-2;

  OhOut = atof(argv[1]);
  RhoIn = atof(argv[2]);
  Rr = atof(argv[3]);
  MAXlevel = atoi(argv[4]);
  tmax = atof(argv[5]);
  zWall = atof(argv[6]);

  Ldomain = zWall+2.+2.*Rr+4.0;

  fprintf(ferr, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, OhOut %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);

  L0=Ldomain;
  origin(-2.0-zWall, 0.0);
  init_grid (1 << (6));

  rho1 = RhoIn; mu1 = MuRin*OhOut;
  rho2 = 1e0; mu2 = OhOut;
  f.sigma = 1.0;

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "dump");

  run();
}

event init(t = 0){
#if _MPI
  if (!restore(file = dumpFile)) {
    fprintf(ferr, "Cannot restored from a dump file!\n");
  }
#else
  if (!restore (file = dumpFile)){
    char filename[60];
    sprintf(filename,"InitialConditionRr-%3.2f.dat", Rr);

    char comm[160];
    sprintf (comm, "scp -r DataFiles/%s .", filename);
    system(comm);

    FILE * fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-8, 1e-8}, MAXlevel).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fraction of the domain. */
    fractions (phi, f);
    foreach(){
      u.x[] = 0.0;
      u.y[] = 0.0;
      // p[] = 2*(1.-f[]);
    }
    dump (file = dumpFile);
    static FILE * fp2;
    fp2 = fopen("InitialConditionStatus.dat","w");
    fprintf(fp2, "Initial condition is written to %s\n", dumpFile);
    fclose(fp2);
    // return 1;
  }
#endif
}

event adapt(i++){
  adapt_wavelet ((scalar *){f, u.x, u.y},
     (double[]){fErr, VelErr, VelErr},
      MAXlevel, MAXlevel-6);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

event end (t = end) {
  fprintf(ferr, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, Oh2 %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);
}

scalar posEq[], posPoles[];
event logWriting (t = 0; t += tsnap2; t <= tmax+tsnap) {

  double ke = 0., wt = 0., xCOM = 0., Vcm = 0.;

  foreach (reduction(+:ke) reduction(+:wt) reduction(+:xCOM) reduction(+:Vcm)){
    ke += 2*pi*y*(0.5*clamp(f[], 0., 1.)*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    xCOM += 2*pi*y*x*clamp(f[], 0.0, 1.0)*sq(Delta);
    Vcm += 2*pi*y*u.x[]*clamp(f[], 0.0, 1.0)*sq(Delta);
    wt += 2*pi*y*clamp(f[], 0.0, 1.0)*sq(Delta);
  }
  xCOM /= wt;
  
  static FILE * fp;

  if (pid() == 0) {
    if (i == 0) {
      fprintf (ferr, "i dt t ke Xc Vcm Re\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, OhOut %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);
      fprintf (fp, "i dt t ke Xc Vcm Re\n");
      fprintf (fp, "%d %g %g %g %g %g\n", i, dt, t, ke, xCOM, Vcm/wt);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g %g %g\n", i, dt, t, ke, xCOM, Vcm/wt);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g %g %g\n", i, dt, t, ke, xCOM, Vcm/wt);
  }

  assert(ke > -1e-10);
  // dump (file = "dump");
  if (ke < 1e-5 && i > 1000){
    if (pid() == 0){
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      fp = fopen ("log", "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
    }
    return 1;
  }
}
