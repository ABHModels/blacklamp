// Code to create jobs for running RELXILL_NK ray-tracing code on cluster
// Creates 30x30 (spin x deformation parameter) jobs in folder "jobs" - must create folder before running this code
// Uses isco_*.dat file to generate jobs. Create isco_*.dat file is isco.cpp code.
// Modify isco_*.dat filename here is necessary.
// Modify fprintf inputs and text after "./photons4trf" as necessary.

#include <stdio.h>
#include <iostream>
#include <stdlib.h>

int main()
{
  double spin, defpar, garb;
  int s, a;
  FILE *prog, *in;
  char filename[128];
  
  in = fopen("isco_a13_shafqat.dat", "r");
  
  s = 0;
  a = 0;
  
	while(fscanf(in, "%lf %lf %lf", &spin, &defpar, &garb) > 0)
  {
		sprintf(filename, "jobs/ring_%d_%d", s, a);
    prog=fopen(filename, "w");

	  fprintf(prog, "#!/bin/bash\n\n#PBS -N ring_%.2f_%.2f\n#PBS -j oe\n#PBS -l nodes=1:ppn=1,walltime=96:00:00\n#PBS -o output/output_%d_%d\n#PBS -e error/error_%d_%d\n\ncd $PBS_O_WORKDIR\ntime ./mainloop_ring %.8f 0.0 0.0 %.8f 0.0 0.0\n", spin, defpar, s, a, s, a, spin, defpar);

	  fclose(prog);
    
    a++;
    if(a > 19)
    {
      a = 0;
      s++;
    }
  }
  
  fclose(in);
    
  return 0;
}
