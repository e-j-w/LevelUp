#include "parse_ENSDF.h"

int strTokCmp(const char * str, const char * cmp, int pos)
{
	if(strcmp(str,"")==0)
		return(strcmp(str,cmp));
	char line[256];
	int i=0;
	char *tok;
	strcpy(line,str); //store the entire line
	tok=strtok (line," ");
	while(i<pos)
		if (tok != NULL)
			{
				tok = strtok (NULL, " ");
				if(tok!=NULL)
					{
						i++;
					}
			}
	return strcmp(tok,cmp);
}

//returns true if two numbers are equal to within some fudge factor
int fudgeNumbers(double num1, double num2, double fudgeFactor)
{
  //printf("val: %f\n",fabs(num1-num2));
  if(fabs(num1-num2)<=fudgeFactor)
    return 1;
  else
    return 0;
}

//if the specified gamma cascade contains the specified level energy,
//returns the step of the cascade containing that energy
int levelInCascade(gamma_cascade * c, double energy)
{
	int i;
	
	for(i=0;i<c->numLevels;i++)
		if(i<MAXCASCDELENGTH)
			if(energy==c->energies[i])
				return i;
				
	return -1;
}

//returns true if the specified gamma cascade contains the specified level energy as its highest level
int lastLevelInCascade(gamma_cascade * c, double energy)
{
	if(c->numLevels>0)
		if(energy==c->energies[c->numLevels-1])
			return 1;

	return 0;
}

void generateCascadeData(gdata *gd)
{
	int i,j,k,l,m,n;// loop indices (!)
	int append=0;
	int ind;
	
  //generate gamma cascade lists
	if(gd->numNucl>=0) //check that indices are valid
	  for(i=0;i<gd->numNucl;i++)
	  	{
	  		gd->nuclData[i].numCascades=0;//initialize properly
				if(gd->nuclData[i].numLevels>=0) //check that indices are valid
					for(j=0;j<gd->nuclData[i].numLevels;j++)
					  for(k=0;k<gd->nuclData[i].levels[j].numGammas;k++)
					  	{
					  		//check whether the level decays to another level
					  		double trialLevelE = gd->nuclData[i].levels[j].energy - gd->nuclData[i].levels[j].gamma_energies[k];
					    	for(l=0;l<j;l++)
						      if(fudgeNumbers(trialLevelE, gd->nuclData[i].levels[l].energy, 2.)==1)
						        {
						        
						        	//printf("-->%s: Cascade detected between levels at energies %f and %f keV.\n",gd->nuclData[i].nuclName,gd->nuclData[i].levels[j].energy,gd->nuclData[i].levels[l].energy);
						        	//printf("Number of cascades: %i\n",gd->nuclData[i].numCascades);
						        
						          append=0;
						          //see if we can add the level to an existing cascade
						          for(m=0;m<gd->nuclData[i].numCascades;m++)
						          	if(m<MAXCASCDESPERNUCL)
						          		if(gd->nuclData[i].cascades[m].numLevels<MAXGAMMASPERLEVEL)
								        		if(lastLevelInCascade(&gd->nuclData[i].cascades[m], gd->nuclData[i].levels[l].energy))
													  	{
													  		//printf("Adding to cascade %i.\n",m+1);
													  		//add level to existing cascade
													  		gd->nuclData[i].cascades[m].energies[gd->nuclData[i].cascades[m].numLevels] = gd->nuclData[i].levels[j].energy;
													  		gd->nuclData[i].cascades[m].numLevels++;
													  		append=1;
													  		break;
													  	}
											if(append==0)
												{
													//see if we can make a new cascade by copying part of an older one
													if((gd->nuclData[i].numCascades+1)<MAXCASCDESPERNUCL)//verify that there is room for a new cascade
														for(m=0;m<gd->nuclData[i].numCascades;m++)
															{
																ind=levelInCascade(&gd->nuclData[i].cascades[m], gd->nuclData[i].levels[l].energy);
								        				if(ind>=0)
								        					{
								        						//printf("Copying existing data from cascade %i.  ind=%i\n",m+1,ind);
								        						append=1;
										      					if((ind+2)<MAXCASCDELENGTH)//verify that there is room for a new level
																			{
																				//copy the cascade
																				gd->nuclData[i].cascades[gd->nuclData[i].numCascades].numLevels=ind+2;
																				for(n=0;n<=ind;n++)
																					gd->nuclData[i].cascades[gd->nuclData[i].numCascades].energies[n]=gd->nuclData[i].cascades[m].energies[n];
																				//add the new level
																				gd->nuclData[i].cascades[gd->nuclData[i].numCascades].energies[ind+1]=gd->nuclData[i].levels[j].energy;
																				//increment cascade counter
																				gd->nuclData[i].numCascades++;
																				break;
																			}
																	}
													  	}
												}
											if(append==0)
												{
													//make a brand new cascade
													if((gd->nuclData[i].numCascades+1)<MAXCASCDESPERNUCL)//verify that there is room for a new cascade
														{
															//printf("Creating new cascade.\n");
															gd->nuclData[i].cascades[gd->nuclData[i].numCascades].numLevels=2;
															gd->nuclData[i].cascades[gd->nuclData[i].numCascades].energies[0]=gd->nuclData[i].levels[l].energy;
															gd->nuclData[i].cascades[gd->nuclData[i].numCascades].energies[1]=gd->nuclData[i].levels[j].energy;
															gd->nuclData[i].numCascades++;
														}
												}
											
						        }
						      /*else
						      	{
						      		if(strcmp(gd->nuclData[i].nuclName,"68SE")==0)
						      			{
										  		printf("j=%i, k=%i, l=%i, numLevels=%i\n",j,k,l,gd->nuclData[i].numLevels);
										  		printf("-->%s: No cascade detected between levels at energies %f and %f keV.\n",gd->nuclData[i].nuclName,gd->nuclData[i].levels[j].energy,gd->nuclData[i].levels[l].energy);
						      			}
						      	}*/
						    }
						    
				//dump cascade data
				/*if(gd->nuclData[i].numCascades>0)
					{
						printf("Cascades generated for nucleus: %s\n",gd->nuclData[i].nuclName);
						for(m=0;m<gd->nuclData[i].numCascades;m++)
							{
								printf("CASCADE %i:\nStep   Level Energy (keV)   Gamma energy (keV)\n",m+1);
								for(n=0;n<gd->nuclData[i].cascades[m].numLevels;n++)
									printf("%i      %f\n",n+1,gd->nuclData[i].cascades[m].energies[n]);
							}
						//if(strcmp(gd->nuclData[i].nuclName,"68SE")==0)
							//getc(stdin);
			  	}*/
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
			  gd->nuclData[i].levels[j].numGammas = 0;
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
						if((strcmp(val[2],"LEVELS")==0)||(strcmp(val[2],"LEVELS,")==0))
							if(gd->numNucl<MAXNUMNUCL)
								{
									gd->numNucl++;
									subSec=0; //we're at the beginning of the entry for this nucleus 
									//printf("Adding gamma data for nucleus %s\n",val[0]);
									strcpy(gd->nuclData[gd->numNucl].nuclName,val[0]);
								}
							
					//add gamma levels
					if(gd->numNucl>=0) //check that indices are valid
						if(strcmp(val[0],gd->nuclData[gd->numNucl].nuclName)==0)
							if(subSec==1)//adopted gamma levels subsection
								if(gd->nuclData[gd->numNucl].numLevels<MAXLEVELSPERNUCL)
									if(strcmp(val[1],"L")==0)
										{
										  //printf("Found gamma level at %f keV.\n",atof(val[2]));
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
										    //printf("-> Found gamma ray with energy %f keV.\n",atof(val[2]));
										    gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].gamma_energies[gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas]=atof(val[2]);
											  gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas++;
											  
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
