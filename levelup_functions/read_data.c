#include "read_data.h"

//function reads an .mca file into a double array and returns the number of spectra read in
int readMCA(FILE * inp, const char * filename, double outHist[NSPECT][S32K])
{
	int i,j;
	int tmpHist[S32K];
	
	//get the number of spectra in the .mca file
  int numSpec=S32K;
  for (i=0;i<numSpec;i++)
    if(fread(tmpHist,S32K*sizeof(int),1,inp)!=1)
      {
        numSpec=i;
        break;
      }
  fclose(inp);
  //printf("number of spectra in file '%s': %i\n",filename,numSpec);
	if((inp=fopen(filename,"r"))==NULL) //reopen the file
    {
      printf("ERROR: Cannot open the input file: %s\n",filename);
      printf("Check that the file exists.\n");
      exit(-1);
    }

	for (i=0;i<numSpec;i++)
		{
			if(fread(tmpHist,S32K*sizeof(int),1,inp)!=1)
				{
					printf("ERROR: Cannot read spectrum %i from the .mca file: %s\n",i,filename);
					printf("Verify that the format and number of spectra in the file are correct.\n");
					exit(-1);
				}
			else
				for(j=0;j<S32K;j++)
					outHist[i][j]=(double)tmpHist[j];
		}
	
	return numSpec;
  
}

//function reads an .fmca file into a double array and returns the number of spectra read in
int readFMCA(FILE * inp, const char * filename, double outHist[NSPECT][S32K])
{
	int i,j;
	float tmpHist[S32K];
	
	//get the number of spectra in the .fmca file
  int numSpec=S32K;
  for (i=0;i<numSpec;i++)
    if(fread(tmpHist,S32K*sizeof(float),1,inp)!=1)
      {
        numSpec=i;
        break;
      }
  fclose(inp);
  //printf("number of spectra in file '%s': %i\n",filename,numSpec);
	if((inp=fopen(filename,"r"))==NULL) //reopen the file
    {
      printf("ERROR: Cannot open the input file: %s\n",filename);
      printf("Check that the file exists.\n");
      exit(-1);
    }

	for (i=0;i<numSpec;i++)
		{
			if(fread(tmpHist,S32K*sizeof(float),1,inp)!=1)
				{
					printf("ERROR: Cannot read spectrum %i from the .fmca file: %s\n",i,filename);
					printf("Verify that the format and number of spectra in the file are correct.\n");
					exit(-1);
				}
			else
				for(j=0;j<S32K;j++)
					outHist[i][j]=(double)tmpHist[j];
		}
	
	return numSpec;
  
}

//function reads an .spe file into a double array and returns the array
int readSPE(FILE * inp, const char * filename, double outHist[NSPECT][S32K])
{
	int i;
  char header[36];
  float inpHist[4096];
  memset(outHist,0,sizeof(*outHist));

	if(fread(header,36,1,inp)!=1)
    {
      printf("ERROR: Cannot read header from the .spe file: %s\n",filename);
      printf("Verify that the format of the file is correct.\n");
      exit(-1);
    }
  if(fread(inpHist,4096*sizeof(float),1,inp)!=1)
    {
      printf("ERROR: Cannot read spectrum from the .spe file: %s\n",filename);
      printf("Verify that the format of the file is correct.\n");
      exit(-1);
    }
  
  //convert input data to double
  for(i=0;i<4096;i++)
  	outHist[0][i]=(double)inpHist[i];
  for(i=4096;i<S32K;i++)
  	outHist[0][i]=0.;
  	
  return 1;
  
}

//reads a file containing spectrum data into an array
//returns the number of spectra read (0 if reading fails)
int readDataFile(const char * filename, double outHist[NSPECT][S32K])
{
	int numSpec=0;
	FILE *inp;
	if((inp=fopen(filename,"r"))==NULL)
    {
      printf("Cannot open the input file: %s\n",filename);
      printf("Check that the file exists.\n");
      //exit(-1);
    }
	else
		{
			const char *dot = strrchr(filename, '.');//get the file extension
			if(strcmp(dot + 1,"mca")==0)
				numSpec=readMCA(inp, filename, outHist);
			else if(strcmp(dot + 1,"fmca")==0)
				numSpec=readFMCA(inp, filename, outHist);
			else if(strcmp(dot + 1,"spe")==0)
				numSpec=readSPE(inp, filename, outHist);
			else
				{
					printf("Improper type of input file: %s\n",filename);
				  printf("Integer array (.mca), float array (.fmca), or radware (.spe) files are supported.\n");
				  //exit(-1);
				}
			fclose(inp);
		}
	
	return numSpec;
}

int readNuclideRankParameters(const char * filename, nuclide_rank_par * nrpar)
{
	// Read parameters from text file
	char *tok;
	char str[256],fullLine[256],parameter[256],value[256];

	FILE *parfile;

	if((parfile=fopen(filename,"r"))==NULL)
		{
			printf("Cannot open the input file: %s\n",filename);
			printf("Check that the file exists.\n");
			return 0;
			//exit(-1);
		}
	else
		{
			while(!(feof(parfile)))//go until the end of file is reached
				{
					if(fgets(str,256,parfile)!=NULL) //get an entire line
						{
							strcpy(fullLine,str);
								tok=strtok(str,"[");
							if(tok!=NULL){
								tok[strcspn(tok, "\r\n")] = 0;//strips newline characters from the string
								strcpy(parameter,tok);
								tok = strtok (NULL, "]");
								if(tok!=NULL){
									tok[strcspn(tok, "\r\n")] = 0;//strips newline characters from the string
									strcpy(value,tok);
									if(strcmp(parameter,"on_shell_closure")==0){
										nrpar->onShellClosure = atof(value);
									}else if(strcmp(parameter,"on_subshell_closure")==0){
										nrpar->onSubshellClosure = atof(value);
									}else if(strcmp(parameter,"radioactive")==0){
										nrpar->radioactive = atof(value);
									}else if(strcmp(parameter,"low_relative_level_density")==0){
										nrpar->lowRelativeLevelDensity = atof(value);
									}else if(strcmp(parameter,"high_relative_level_density")==0){
										nrpar->highRelativeLevelDensity = atof(value);
									}else if(strcmp(parameter,"max_distance_from_stability")==0){
										nrpar->maxDistFromStability = atoi(value);
									}else if(strcmp(parameter,"exclude_stable_nuclides")==0){
										if(strcmp(value,"yes")==0){
											nrpar->excludeStable = 1;
										}else{
											nrpar->excludeStable = 0;
										}
									}
								}
							}
						}
				}

			printf("\nFile '%s' read sucessfully!\n",filename);
			/*printf("Listing of weights:\n");
			printf("on_shell_closure = %.2f \n",nrpar->onShellClosure);
			printf("on_subshell_closure = %.2f \n",nrpar->onSubshellClosure);
			printf("radioactive = %.2f \n",nrpar->radioactive);
			printf("low_relative_level_density = %.2f \n",nrpar->lowRelativeLevelDensity);
			printf("high_relative_level_density = %.2f \n",nrpar->highRelativeLevelDensity);
			printf("Other parameters:\n");
			printf("max_distance_from_stability = %i \n",nrpar->maxDistFromStability);
			if(nrpar->excludeStable == 1){
				printf("Excluding stable nucluides.\n");
			}else{
				printf("Not excluding stable nucluides.\n");
			}*/
		}

	fclose(parfile);
	return 1;
}


