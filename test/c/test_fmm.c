#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "fmm.h"
#include <assert.h>

int FC_MAIN (int argc, char **argv) {
  FILE *input_file;
  char *input_file_name;
  long int total = 114537;
  double *pos;
  double *charges;
  long int dim = 3;
  int i = 0, j= 0, k=0;
  double *forces;
  double coulomb_energy;
  char *parameterstring;

  /* fmm parameter */
  FCSInput input_handle;
  FCSOutput output_handle;
  FCSError err;

  if(argc == 1) {
    printf("No input file was specified!\n");
    exit(-1);
  }
  
  if(argc == 2) {
    input_file_name = argv[1];
    if((input_file = fopen(input_file_name, "r")) == NULL) {
      printf("No input file was found!\n");
      exit(-1);
    }
  }

  pos = (double *)malloc(total*dim*sizeof(double));
  assert(pos);
  charges = (double *)malloc(total*sizeof(double));
  assert(charges);
  for(i=0; i<total; ++i) {
    fscanf(input_file, " %lf", charges+i);
    for(j = 0; j < dim; ++j)
      fscanf(input_file, " %lf", (pos + i*dim + j));
  }
  
  printf("i: %d\n", i);
  printf("done reading data\n");
 
  
  err = fcsInput_create(&input_handle,
			total,
			pos,
			dim,
			charges);

  printf(fcsError_getMessage(err));
  if(fcsError_getErrorCode(err) != FCS_SUCCESS)
    {
      fclose(input_file);
      exit(-1);
    }


  err = fcsOutput_create(&output_handle,
			 total,
			 dim);

  printf(fcsError_getMessage(err));

 
  if(fcsError_getErrorCode(err) != FCS_SUCCESS)
    {
      fclose(input_file);
      exit(-1);
    }
  
  printf("Run fmm...\n");
  /* allocate parameter string for storing 
   * optional parameters in form "param1,param2,param3"
   */
  parameterstring = (char *)calloc(strlen(FMM_COULOMB_ENERGY)+1+
				   strlen(FMM_ERR_TYPE)+1+
				   strlen(FMM_ERR_BOUND)+1, sizeof(char));
  strcat(parameterstring, FMM_COULOMB_ENERGY);
  strcat(parameterstring, ",");
  strcat(parameterstring, FMM_ERR_TYPE);
  strcat(parameterstring, ",");
  strcat(parameterstring, FMM_ERR_BOUND);

  err = run_fmm_with_opt_param_(input_handle, output_handle,
				parameterstring, &coulomb_energy, 1, 1.e-5);

  forces = fcsOutput_getForces(output_handle);
  printf("%36.31e, %36.31e, %36.31e, %36.31e\n",
	coulomb_energy, forces[0], forces[1],forces[2]);
  printf(fcsError_getMessage(err));
 
  if(fcsError_getErrorCode(err) != FMM_OK)
    {
      fclose(input_file);
      exit(-1);
    }

  fclose(input_file);

  fcsInput_destroy(&input_handle);
  fcsOutput_destroy(&output_handle);
  fcsError_destroy(&err);

  exit(0);

}
