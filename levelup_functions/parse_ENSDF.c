#include "parse_ENSDF.h"

int fudgeNumbers(double num1, double num2)
{
  double fudgeFactor = 2.;
  //printf("val: %f\n",fabs(num1-num2));
  if(fabs(num1-num2)<=fudgeFactor)
    return 1;
  else
    return 0;
}

void generateCascadeData(gdata *gd)
{
  int i,j,k,l;
  
  //generate gamma cascade lists
	if(gd->numNucl>=0) //check that indices are valid
	  for(i=0;i<gd->numNucl;i++)
		  if(gd->nuclData[i].numLevels>=0) //check that indices are valid
			  for(j=0;j<gd->nuclData[i].numLevels;j++)
			    for(k=0;k<gd->nuclData[i].levels[j].numGammas;k++)
			      for(l=0;l<j;l++)
				      {
				        //check whether the level decays to another level
				        double trialLevelE = gd->nuclData[i].levels[j].energy - gd->nuclData[i].levels[j].gamma_energies[k];
				        if(fudgeNumbers(trialLevelE, gd->nuclData[i].levels[l].energy)==1)
				          {
				            printf("-->%s: Cascade detected between levels at energies %f and %f keV.\n",gd->nuclData[i].nuclName,gd->nuclData[i].levels[j].energy,gd->nuclData[i].levels[l].energy);
				          }
				      }
}

//set initial databae values prior to importing data
void initialize_database(gdata * gd) 
{
	int i,j;
	gd->numNucl = -1;
	for(i=0;i<MAXNUMNUCL;i++)
		{
			gd->nuclData[i].numLevels = -1;
			for(j=0;j<MAXLEVELSPERNUCL;j++)
			  gd->nuclData[i].levels[j].numGammas = -1;
		}
}


//function reads parameter files for the topspek code
void readENSDFFile(const char * fileName, gdata * gd) 
{

  FILE *config;
  char *tok;
  char str[256];//string to be read from file (will be tokenized)
  char line[256],val[MAXNUMVALS][256];
  int tokPos;//position when tokenizing
  
  //subsection of the entry for a particular nucleus that the parser is at
  //each nucleus has multiple entries, including adopted gammas, and gammas 
  //associated with a particlular reaction mechanism
  int subSec; 
  
  //open the file and read all parameters
  if((config=fopen(fileName,"r"))==NULL)
    {
      printf("WARNING: Cannot open the ENSDF file %s\n",fileName);
      return;
    }
  while(!(feof(config)))//go until the end of file is reached
    {
			if(fgets(str,256,config)!=NULL) //get an entire line
				{
					//sscanf(str,"%s",str1);
					strcpy(line,str); //store the entire line
					tok=strtok (str," ");
					tokPos=0;
					strcpy(val[tokPos],tok);
					while (tok != NULL)
					{
						tok = strtok (NULL, " ");
						if(tok!=NULL)
							{
								tokPos++;
								if(tokPos<MAXNUMVALS)
									strcpy(val[tokPos],tok);
								else
									break;
							}
					}
					
					//increment the nucleus if a new nucleus is found
					if(strcmp(val[1],"ADOPTED")==0)
						if(strcmp(val[3],"GAMMAS")==0)
							if(gd->numNucl<MAXNUMNUCL)
								{
									gd->numNucl++;
									subSec=0; //we're at the beginning of the entry for this nucleus 
									printf("Adding gamma data for nucleus %s\n",val[0]);
									strcpy(gd->nuclData[gd->numNucl].nuclName,val[0]);
								}
							
					//add gamma levels
					if(gd->numNucl>=0) //check that indices are valid
						if(strcmp(val[0],gd->nuclData[gd->numNucl].nuclName)==0)
							if(subSec==1)//adopted gamma levels subsection
								if(gd->nuclData[gd->numNucl].numLevels<MAXLEVELSPERNUCL)
									if(strcmp(val[1],"L")==0)
										{
										  printf("Found gamma level at %f keV.\n",atof(val[2]));
											gd->nuclData[gd->numNucl].numLevels++;
											gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].energy=atof(val[2]);
											
											
											
											//add gamma to cascade if neccessary
											/*for(i=0;i<gd->nuclData[gd->numNucl].numLevels-1;i++)
											  {
											    //check whether the level decays to another level
											    if(fudgeNumbers(gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].energy, gd->nuclData[gd->numNucl].levels[i].energy)==1)
											  }*/
										}
					
					//add gamma rays
					if(gd->numNucl>=0) //check that indices are valid
					  if(gd->nuclData[gd->numNucl].numLevels>=0) //check that indices are valid
						  if(strcmp(val[0],gd->nuclData[gd->numNucl].nuclName)==0)
							  if(subSec==1)//adopted gamma levels subsection
								  if(gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas<MAXGAMMASPERLEVEL)
									  if(strcmp(val[1],"G")==0)
										  {
										    printf("-> Found gamma ray with energy %f keV.\n",atof(val[2]));
											  gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas++;
											  gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].gamma_energies[gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas]=atof(val[2]);
											  
										  }
										  
					
					
					
					//if neccesary, increment which subsection we're on
					if(gd->numNucl>=0)
						if(strcmp(val[0],gd->nuclData[gd->numNucl].nuclName)==0)
							if(strcmp(val[1],"H")==0)
								{
									subSec++;
								}
					
				}
		}
	fclose(config);
	
	if(gd->numNucl>=MAXNUMNUCL)
		{
			printf("ERROR: Attempted to import data for too many nuclei.  Increase the value of MAXNUMNUCL in levelup.h\n");
			exit(-1);
		}
	
	printf("Finished reading ENSDF file: %s\n",fileName);
  
}
