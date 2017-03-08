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

/*//function reads an .spe file into a double array and returns the array
void readSPE(FILE * inp, const char * filename, const int numSpec, double outHist[NSPECT][S32K])
{
	int i;
  char header[36];
  float inpHist[S32K];
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
  for(i=0;i<S32K;i++)
  	outHist[0][i]=(double)inpHist[i];
  
  //copy data to each spectrum
  for(i=1;i<numSpec;i++)
  	memcpy(outHist[i],outHist[0],S32K*sizeof(double));
  
}*/

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
			/*else if(strcmp(dot + 1,"spe")==0)
				readSPE(inp, filename, numSpec, outHist);*/
			else
				{
					printf("Improper type of input file: %s\n",filename);
				  printf("Integer array (.mca) files are supported.\n");
				  //exit(-1);
				}
			fclose(inp);
		}
	
	return numSpec;
}
