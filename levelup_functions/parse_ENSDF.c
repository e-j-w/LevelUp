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
										//see if we can add the level to existing cascade(s)
										for(m=0;m<gd->nuclData[i].numCascades;m++)
									  	if(m<MAXCASCDESPERNUCL)
									  		if(gd->nuclData[i].cascades[m].numLevels<MAXGAMMASPERLEVEL)
													if(lastLevelInCascade(&gd->nuclData[i].cascades[m], gd->nuclData[i].levels[l].energy))
														if(gd->nuclData[i].cascades[m].numLevels+1<=MAXCASCDELENGTH)
													  	{
													  		//printf("Adding to cascade %i.\n",m+1);
													  		//add level to existing cascade
													  		gd->nuclData[i].cascades[m].energies[gd->nuclData[i].cascades[m].numLevels] = gd->nuclData[i].levels[j].energy;
													  		gd->nuclData[i].cascades[m].numLevels++;
													  		append=1;
													  	}
										if(append==0)
											{
												//see if we can make a new cascade by copying part of an older one
												if((gd->nuclData[i].numCascades+1)<=MAXCASCDESPERNUCL)//verify that there is room for a new cascade
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
												if((gd->nuclData[i].numCascades+1)<=MAXCASCDESPERNUCL)//verify that there is room for a new cascade
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
	
	
	//generate cascade gamma ray energies
	if(gd->numNucl>=0) //check that indices are valid
	  for(i=0;i<gd->numNucl;i++)
	  	for(j=0;j<gd->nuclData[i].numCascades;j++)
	  		for(k=0;k<gd->nuclData[i].cascades[j].numLevels;k++)
	  			{
	  				if(k==0)
	  					gd->nuclData[i].cascades[j].gammaEnergies[k]=0.;
	  				else
	  					gd->nuclData[i].cascades[j].gammaEnergies[k]=gd->nuclData[i].cascades[j].energies[k] - gd->nuclData[i].cascades[j].energies[k-1];
	  			}
	
	
}


void getNuclNZ(nucl * nuc)
{
	char str[256];
	int A,Z,N;
	char *tok;
	
	//get mass number
	strcpy(str,nuc->nuclName); //copy the nucleus name
	tok=strtok (str,"ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	A=atoi(tok);
	
	//get proton number
	strcpy(str,nuc->nuclName); //copy the nucleus name
	tok=strtok (str,"0123456789");
	if(strcmp(tok,"H")==0)
		Z=1;
	else if(strcmp(tok,"HE")==0)
		Z=2;
	else if(strcmp(tok,"LI")==0)
		Z=3;
	else if(strcmp(tok,"BE")==0)
		Z=4;
	else if(strcmp(tok,"B")==0)
		Z=5;
	else if(strcmp(tok,"C")==0)
		Z=6;
	else if(strcmp(tok,"N")==0)
		Z=7;
	else if(strcmp(tok,"O")==0)
		Z=8;
	else if(strcmp(tok,"F")==0)
		Z=9;
	else if(strcmp(tok,"NE")==0)
		Z=10;
	else if(strcmp(tok,"NA")==0)
		Z=11;
	else if(strcmp(tok,"MG")==0)
		Z=12;
	else if(strcmp(tok,"AL")==0)
		Z=13;
	else if(strcmp(tok,"SI")==0)
		Z=14;
	else if(strcmp(tok,"P")==0)
		Z=15;
	else if(strcmp(tok,"S")==0)
		Z=16;
	else if(strcmp(tok,"CL")==0)
		Z=17;
	else if(strcmp(tok,"AR")==0)
		Z=18;
	else if(strcmp(tok,"K")==0)
		Z=19;
	else if(strcmp(tok,"CA")==0)
		Z=20;
	else if(strcmp(tok,"SC")==0)
		Z=21;
	else if(strcmp(tok,"TI")==0)
		Z=22;
	else if(strcmp(tok,"V")==0)
		Z=23;
	else if(strcmp(tok,"CR")==0)
		Z=24;
	else if(strcmp(tok,"MN")==0)
		Z=25;
	else if(strcmp(tok,"FE")==0)
		Z=26;
	else if(strcmp(tok,"CO")==0)
		Z=27;
	else if(strcmp(tok,"NI")==0)
		Z=28;
	else if(strcmp(tok,"CU")==0)
		Z=29;
	else if(strcmp(tok,"ZN")==0)
		Z=30;
	else if(strcmp(tok,"GA")==0)
		Z=31;
	else if(strcmp(tok,"GE")==0)
		Z=32;
	else if(strcmp(tok,"AS")==0)
		Z=33;
	else if(strcmp(tok,"SE")==0)
		Z=34;
	else if(strcmp(tok,"BR")==0)
		Z=35;
	else if(strcmp(tok,"KR")==0)
		Z=36;
	else if(strcmp(tok,"RB")==0)
		Z=37;
	else if(strcmp(tok,"SR")==0)
		Z=38;
	else if(strcmp(tok,"Y")==0)
		Z=39;
	else if(strcmp(tok,"ZR")==0)
		Z=40;
	else if(strcmp(tok,"NB")==0)
		Z=41;
	else if(strcmp(tok,"MO")==0)
		Z=42;
	else if(strcmp(tok,"TC")==0)
		Z=43;
	else if(strcmp(tok,"RU")==0)
		Z=44;
	else if(strcmp(tok,"RH")==0)
		Z=45;
	else if(strcmp(tok,"PD")==0)
		Z=46;
	else if(strcmp(tok,"AG")==0)
		Z=47;
	else if(strcmp(tok,"CD")==0)
		Z=48;
	else if(strcmp(tok,"IN")==0)
		Z=49;
	else if(strcmp(tok,"SN")==0)
		Z=50;
	else if(strcmp(tok,"SB")==0)
		Z=51;
	else if(strcmp(tok,"TE")==0)
		Z=52;
	else if(strcmp(tok,"I")==0)
		Z=53;
	else if(strcmp(tok,"XE")==0)
		Z=54;
	else if(strcmp(tok,"CS")==0)
		Z=55;
	else if(strcmp(tok,"BA")==0)
		Z=56;
	else if(strcmp(tok,"LA")==0)
		Z=57;
	else if(strcmp(tok,"CE")==0)
		Z=58;
	else if(strcmp(tok,"PR")==0)
		Z=59;
	else if(strcmp(tok,"ND")==0)
		Z=60;
	else if(strcmp(tok,"PM")==0)
		Z=61;
	else if(strcmp(tok,"SM")==0)
		Z=62;
	else if(strcmp(tok,"EU")==0)
		Z=63;
	else if(strcmp(tok,"GD")==0)
		Z=64;
	else if(strcmp(tok,"TB")==0)
		Z=65;
	else if(strcmp(tok,"DY")==0)
		Z=66;
	else if(strcmp(tok,"HO")==0)
		Z=67;
	else if(strcmp(tok,"ER")==0)
		Z=68;
	else if(strcmp(tok,"TM")==0)
		Z=69;
	else if(strcmp(tok,"YB")==0)
		Z=70;
	else if(strcmp(tok,"LU")==0)
		Z=71;
	else if(strcmp(tok,"HF")==0)
		Z=72;
	else if(strcmp(tok,"TA")==0)
		Z=73;
	else if(strcmp(tok,"W")==0)
		Z=74;
	else if(strcmp(tok,"RE")==0)
		Z=75;
	else if(strcmp(tok,"OS")==0)
		Z=76;
	else if(strcmp(tok,"IR")==0)
		Z=77;
	else if(strcmp(tok,"PT")==0)
		Z=78;
	else if(strcmp(tok,"AU")==0)
		Z=79;
	else if(strcmp(tok,"HG")==0)
		Z=80;
	else if(strcmp(tok,"TL")==0)
		Z=81;
	else if(strcmp(tok,"PB")==0)
		Z=82;
	else if(strcmp(tok,"BI")==0)
		Z=83;
	else if(strcmp(tok,"PO")==0)
		Z=84;
	else if(strcmp(tok,"AT")==0)
		Z=85;
	else if(strcmp(tok,"RN")==0)
		Z=86;
	else if(strcmp(tok,"FR")==0)
		Z=87;
	else if(strcmp(tok,"RA")==0)
		Z=88;
	else if(strcmp(tok,"AC")==0)
		Z=89;
	else if(strcmp(tok,"TH")==0)
		Z=90;
	else if(strcmp(tok,"PA")==0)
		Z=91;
	else if(strcmp(tok,"U")==0)
		Z=92;
	else if(strcmp(tok,"NP")==0)
		Z=93;
	else if(strcmp(tok,"PU")==0)
		Z=94;
	else if(strcmp(tok,"AM")==0)
		Z=95;
	else if(strcmp(tok,"CM")==0)
		Z=96;
	else if(strcmp(tok,"BK")==0)
		Z=97;
	else if(strcmp(tok,"CF")==0)
		Z=98;
	else if(strcmp(tok,"ES")==0)
		Z=99;
	else if(strcmp(tok,"FM")==0)
		Z=100;
	else if(strcmp(tok,"MD")==0)
		Z=101;
	else if(strcmp(tok,"NO")==0)
		Z=102;
	else if(strcmp(tok,"LR")==0)
		Z=103;
	else if(strcmp(tok,"RF")==0)
		Z=104;
	else if(strcmp(tok,"DB")==0)
		Z=105;
	else if(strcmp(tok,"SG")==0)
		Z=106;
	else if(strcmp(tok,"BH")==0)
		Z=107;
	else if(strcmp(tok,"HS")==0)
		Z=108;
	else if(strcmp(tok,"MT")==0)
		Z=109;
	else if(strcmp(tok,"DS")==0)
		Z=110;
	else if(strcmp(tok,"RG")==0)
		Z=111;
	else
		Z=-1;
	
	//get neutron number
	if(strcmp(tok,"NN")==0)
		{
			N=2;
			Z=A-N;
		}
	else
		{
			N=A-Z;
		}
	
	nuc->N=N;
	nuc->Z=Z;
	
}





//set initial databae values prior to importing data
void initialize_database(gdata * gd) 
{
	int i;
	
	memset(gd,0,sizeof(gdata));
	
	gd->numNucl = -1;
	for(i=0;i<MAXNUMNUCL;i++)
		{
			gd->nuclData[i].numLevels = -1;
		}
}


//function to parse ENSDF data files
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
  int subSec=0; 
  
  //open the file and read all parameters
  if((config=fopen(fileName,"r"))==NULL)
    {
      //printf("WARNING: Cannot open the ENSDF file %s\n",fileName);
      return;
    }
  while(!(feof(config)))//go until the end of file is reached
    {
			if(fgets(str,256,config)!=NULL) //get an entire line
				{
					//sscanf(str,"%s",str1);
					strcpy(line,str); //store the entire line
					if(strcmp(str,"")==0)
						{
							subSec++; //empty line, increment which subsection we're on
						}
					else
						{
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
									getNuclNZ(&gd->nuclData[gd->numNucl]); //get N and Z
								}
							
					//add gamma levels
					if(gd->numNucl>=0) //check that indices are valid
						if(strcmp(val[0],gd->nuclData[gd->numNucl].nuclName)==0)
							if(subSec==0)//adopted gamma levels subsection
								if(gd->nuclData[gd->numNucl].numLevels<MAXLEVELSPERNUCL)
									if(strcmp(val[1],"L")==0)
										{
											
										  //printf("Found gamma level at %f keV.\n",atof(val[2]));
											gd->nuclData[gd->numNucl].numLevels++;
											double levelE = atof(val[2]);
											if(gd->nuclData[gd->numNucl].numLevels==0)
												{
													//the level energy represents a new level
													gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].energy=levelE;
												}
											else if(levelE>gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels-1].energy)
												{
													//the level energy represents a new level
													gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].energy=levelE;
												}
											else
												{
													//the level energy is a repeat
													gd->nuclData[gd->numNucl].numLevels--;
												}
											
											
											
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
							  if(subSec==0)//adopted gamma levels subsection
								  if(gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas<MAXGAMMASPERLEVEL)
									  if(strcmp(val[1],"G")==0)
										  {
										    //printf("-> Found gamma ray with energy %f keV.\n",atof(val[2]));
										    gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].gamma_energies[gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas]=atof(val[2]);
										    gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].gamma_intensities[gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas]=atof(val[4]);
											  gd->nuclData[gd->numNucl].levels[gd->nuclData[gd->numNucl].numLevels].numGammas++;
											  
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
