//definitions
#include "levelup.h"
//functions and logic
#include "read_data.c"
#include "peak_find.c"
#include "parse_ENSDF.c"
#include "db_ops.c"



int main(int argc, char *argv[])
{

	printf("\nLevelUp - There Is No Hope\n--------------------------\n\n");

	if(!getenv("ENSDF")) 
		{
			printf("\nThe environment variable ENSDF is not set.\nWithout it, this program cannot find the ENSDF data files.  Set the envirnoment variable to point to the directory containing ENSDF data files using the following command (replace the path with the appropriate one):\n\n");
			printf ("export ENSDF=/path/to/ensdf/\n\n");
			exit(-1);
		}
  
	FILE *db;
	//allocate structures
	ndata *nd=(ndata*)malloc(sizeof(ndata));
	initialize_database(nd);
	
	
	//process the ENSDF database file
	char fileName[256];
	int i=0;
	int j=0;
	strcpy(fileName,"");
	strcat(fileName,getenv("ENSDF"));
	strcat(fileName,"ensdf_db");
	if((db=fopen(fileName,"r"))==NULL){
		printf("ENSDF database file not found, rebuilding database...\n");
		rebuildDatabase(nd,db);
		db=fopen(fileName,"r");//reopen file for reading
	}
	while(i!=1){	
		i=fread(nd,sizeof(ndata),1,db);
		if(i==1)
			{
				printf("ENSDF database retrieved from disk.  Data for %i nuclei found.\n",nd->numNucl);
			}
		else
			{
				j++;
				printf("ERROR: ENSDF database could not be read from file '%s'.  Rebuilding (attempt %i)...\n",fileName,j);
				rebuildDatabase(nd,db);
				db=fopen(fileName,"r");//reopen file for reading
				
			}
		if(j>3)
			{
				printf("ERROR: ENSDF database could not be built.  Exiting...\n");
				exit(-1);
			}

	}
	fclose(db);
	
  rl_bind_key('\t', rl_complete); // Configure readline to auto-complete on tab key

	//process user commands
	printf("\nEnter command.  Enter 'help' for a list of commands.\n");
	char cmd2[256];
	char *tok;
	while(1){
		//setup readline prompt and add commands to history
		char* cmd = readline("LevelUp > ");
		add_history(cmd);

		//cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
		if(strcmp(cmd,"help")==0){
			printf("\nCommand list:\n\n");
			printf("  casc NUCL - prints cascade data for the specified nucleus\n");
			printf("              (NUCL is the nucleus name, eg. '68SE')\n\n");
			printf("  id SP - finds nuclei with cascades which match peaks in the\n");
			printf("          specified spectrum (SP is the spectrum filename, eg.\n");
			printf("          'spectrum.mca')\n\n");
			printf("  qsp SP - finds nuclei corresponding to a specific gamma ray energy\n");
			printf("           in the specified spectrum, by finding cascades in the\n");
			printf("           spectrum which contain a gamma ray at or near that energy\n\n");
			printf("  lev NUCL - prints level data for the specified nucleus (can\n"); 
			printf("             optionally specify the number of levels to show as a\n");
			printf("             second parameter)\n\n");
			printf("  ol NUCL1 NUCL2 - finds overlapping gamma rays in the two\n"); 
			printf("                   specified nuclei\n\n");
			printf("  olr NUCL1 NUCL2 - finds overlapping gamma rays in the in the\n");
			printf("                    nucleus NUCL1 and nuclei in the region of\n"); 
			printf("                    the nucleus NUCL2 (a search width will be\n"); 
			printf("                    requested)\n\n");
			printf("  pfr NUCL - finds gamma rays in the in the region of the nucleus\n"); 
			printf("             NUCL which match an energy specified by the user\n\n"); 
			printf("  nz NUCL - show N, Z numbers for the specified nucleus\n\n");
			printf("  listnuc, ln - list names of nuclei in the ENSDF database\n\n");
			printf("  findcasc, fc - find nuclei which match a cascade that you enter\n\n");
			printf("  rebuild - rebuild the ENSDF database from ENSDF files on disk\n\n");
			printf("  help - list commands\n\n");
			printf("  exit, quit, ^C - exit the program\n\n");
		}else if((strcmp(cmd,"exit")==0)||(strcmp(cmd,"quit")==0)){
			printf("Exiting...\n");
			exit(0);
		}else if(strTokCmp(cmd,"casc",0)==0){
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No nucleus specified.\n");
			else
				{	
					int nucl=nameToNuclIndex(tok,nd);
					if(nucl<0)
						printf("Unknown nucleus: %s\n",tok);
					else
						showCascadeData(nucl,nd);
				}
		}else if(strTokCmp(cmd,"lev",0)==0){
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No nucleus specified.\n");
			else
				{
					int nucl=nameToNuclIndex(tok,nd);
					if(nucl<0)
						printf("Unknown nucleus: %s\n",tok);
					else if((tok = strtok (NULL, " "))==NULL)//read the 3rd entry in the command
						{
							showLevelData(nucl,nd,0);
						}
					else
						{
							if(atoi(tok)<=0)
								printf("Invalid number of levels specified: %s\n",tok);
							else
								showLevelData(nucl,nd,atoi(tok));
						}
				}
		}else if(strTokCmp(cmd,"ol",0)==0){
			int findOL=1;
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No nuclei specified.\n");
			else{
				int nucl1=nameToNuclIndex(tok,nd);
				if(nucl1<0){
					printf("Unknown nucleus: %s\n",tok);
					findOL=0;
				}else if((tok = strtok (NULL, " "))==NULL){//read the 3rd entry in the command
					printf("Two nuclei must be specified.\n");
				}else{
					int nucl2=nameToNuclIndex(tok,nd);
					if(nucl2<0){
						printf("Unknown nucleus: %s\n",tok);
						findOL=0;
					}
					if(findOL){
						findOverlappingLevels(nucl1,nucl2,nd);
					}
				}
			}
		}else if(strTokCmp(cmd,"olr",0)==0){
			int findOL=1;
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No nuclei specified.\n");
			else
				{
					int nucl1=nameToNuclIndex(tok,nd);
					if(nucl1<0)
						{
							printf("Unknown nucleus: %s\n",tok);
							findOL=0;
						}
					else if((tok = strtok (NULL, " "))==NULL)//read the 3rd entry in the command
						{
							printf("Two nuclei must be specified.\n");
						}
					else
						{
							int nucl2=nameToNuclIndex(tok,nd);
							if(nucl2<0)
								{
									printf("Unknown nucleus: %s\n",tok);
									findOL=0;
								}
							else if((tok = strtok (NULL, " "))==NULL)//read the 4th entry in the command
								{
									printf("Enter the width/height (in nucleon number) of the region in the nuclear chart to search: ");
									fgets(cmd,256,stdin);
									cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
								}
							else
								{
									strcpy(cmd,tok);
								}
							if (findOL)
								{
									int dim=atoi(cmd);
									if(dim<=0)
										printf("Invalid search width/height: %i\n",dim);
									else
										findOverlappingLevelsInRegion(nucl1,nucl2,dim,nd);
								}
						}
				}
		}else if(strTokCmp(cmd,"pfr",0)==0){
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No nucleus specified.\n");
			else{
				int nucl=nameToNuclIndex(tok,nd);
				if(nucl<0){
						printf("Unknown nucleus: %s\n",tok);
				}else{
					if((tok = strtok (NULL, " "))==NULL){ //read the 3rd entry in the command
						printf("Enter the width/height (in nucleon number) of the region in the nuclear chart to search: ");
						fgets(cmd,256,stdin);
						cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
					}else{
						strcpy(cmd,tok);
					}
					int dim=atoi(cmd);
					if(dim<=0){
						printf("Invalid search width/height: %i\n",dim);
					}else{
						if((tok = strtok (NULL, " "))==NULL){//read the 4th entry in the command
							printf("Enter the energy of the gamma ray to search for: ");
							fgets(cmd,256,stdin);
							cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
						}else{
							strcpy(cmd,tok);
						}
						double energy=atof(cmd);
						if(energy<=0){
							printf("Invalid gamma ray energy: %5.1f\n",energy);
						}else{
							findLevelInRegion(energy,nucl,dim,nd);
						}
					}
				}
			}
		}else if(strTokCmp(cmd,"nz",0)==0){
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No nucleus specified.\n");
			else{
				int nucl=nameToNuclIndex(tok,nd);
				if(nucl<0)
					printf("Unknown nucleus: %s\n",tok);
				else
					showNZ(nucl,nd);
			}
		}else if((strcmp(cmd,"listnuc")==0)||(strcmp(cmd,"ln")==0)){
			showNuclNames(nd);
		}else if((strcmp(cmd,"findcasc")==0)||(strcmp(cmd,"fc")==0)){
			int validE=1;
			gamma_cascade *c=(gamma_cascade*)malloc(sizeof(gamma_cascade));
			printf("Enter the number of gammas in the cascade: ");
			fgets(cmd,256,stdin);
			cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
			if(atoi(cmd)<=0)
				{
					validE=0;
					printf("Invalid number of gammas.  Returning...\n");
				}
			else
				c->numLevels=atoi(cmd);
			if(c->numLevels>MAXCASCDELENGTH)
				{
					printf("Maximum cascade length (%i) exceeded.  Setting the number of gammas to %i.\n",MAXCASCDELENGTH,MAXCASCDELENGTH);
					c->numLevels=MAXCASCDELENGTH;
				}
			for(i=0;i<c->numLevels;i++)
				{
					printf("Enter the energy of gamma ray %i [keV]: ",i+1);
					fgets(cmd,256,stdin);
					cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
					if(atof(cmd)<=0.)
						{
							validE=0;
							printf("Invalid energy.  Returning...\n");
							break;
						}
					else
						c->gammaEnergies[i]=atof(cmd);
				}
			if(validE==1)
				findCascade(c,c->numLevels,nd);
			free(c);
		}else if(strcmp(cmd,"rebuild")==0){
			printf("Rebuilding ENSDF database...\n");
			strcpy(fileName,"");
			strcat(fileName,getenv("ENSDF"));
			strcat(fileName,"ensdf_db");
			if((db=fopen(fileName,"r"))!=NULL){
				printf("ENSDF database file not found, rebuilding database...\n");
				rebuildDatabase(nd,db);
				db=fopen(fileName,"r");//reopen file for reading
				i=fread(nd,sizeof(ndata),1,db);
				while(i!=1)
					{	
						i=fread(nd,sizeof(ndata),1,db);
						if(i==1)
							{
								printf("ENSDF database retrieved from disk.  Data for %i nuclei found.\n",nd->numNucl);
							}
						else
							{
								j++;
								printf("ERROR: ENSDF database could not be read from file '%s'.  Rebuilding (attempt %i)...\n",fileName,j);
								rebuildDatabase(nd,db);
								db=fopen(fileName,"r");//reopen file for reading
				
							}
						if(j>3)
							{
								printf("ERROR: ENSDF database could not be built.  Exiting...\n");
								exit(-1);
							}
					}
				fclose(db);
			}else{
				printf("WARNING: could not rebuild database, database file %s is not readable!\n",fileName);
			}
		}else if(strTokCmp(cmd,"id",0)==0){
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No filename specified.\n");
			else
				{
					printf("Identifying species in spectrum data file: %s.\n",tok);
					inp_sp *sp=(inp_sp*)malloc(sizeof(inp_sp));
					sp->numSpectra=readDataFile(tok, sp->hist);
					if(sp->numSpectra>0)
						{
							//sum mca data
							for(i=0;i<S32K;i++)
								sp->sumHist[i]=0.;
							for(i=0;i<sp->numSpectra;i++)
								for(j=0;j<S32K;j++)
									sp->sumHist[j]+=sp->hist[i][j];
							
							if((tok = strtok (NULL, " "))==NULL)//read the 3rd entry in the command
								{
									printf("Enter the energy contraction [keV / channel]: ");
									fgets(cmd,256,stdin);
									cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
								}
							else
								{
									strcpy(cmd,tok);//3rd entry can be used as energy contraction
								}
							double cntr=atof(cmd);
							if(cntr<=0)
								{
									printf("Invalid energy contraction.  Returning...\n");
								}
							else
								{
									//find peak positions
									peak_fit_par pfpar = findPeak(sp->sumHist,cntr,5.,20);
									reportPeakPositions(&pfpar);
									//find matching cascades
									findCascadeInSpec(&pfpar,nd);
								}
						}
					free(sp);
				}
		}else if(strTokCmp(cmd,"qsp",0)==0){
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No filename specified.\n");
			else{
				printf("Finding gamma rays in spectrum data file: %s.\n",tok);
				inp_sp *sp=(inp_sp*)malloc(sizeof(inp_sp));
				sp->numSpectra=readDataFile(tok, sp->hist);
				if(sp->numSpectra>0)
					{
						//sum mca data
						for(i=0;i<S32K;i++)
							sp->sumHist[i]=0.;
						for(i=0;i<sp->numSpectra;i++)
							for(j=0;j<S32K;j++)
								sp->sumHist[j]+=sp->hist[i][j];
						
						if((tok = strtok (NULL, " "))==NULL)//read the 3rd entry in the command
							{
								printf("Enter the energy contraction [keV / channel]: ");
								fgets(cmd,256,stdin);
								cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
							}
						else
							{
								strcpy(cmd,tok);//3rd entry can be used as energy contraction
							}
						double cntr=atof(cmd);
						if(cntr<=0)
							{
								printf("Invalid energy contraction.  Returning...\n");
							}
						else
							{
								
								if((tok = strtok (NULL, " "))==NULL)//read the 4th entry in the command
									{
										printf("Enter the gamma ray energy to check [keV]: ");
										fgets(cmd,256,stdin);
										cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
									}
								else
									{
										strcpy(cmd,tok);//3rd entry can be used as energy contraction
									}
								double ge=atof(cmd);
								if(ge<=0)
									{
										printf("Invalid gamma ray energy.  Returning...\n");
									}
								else
									{
										//find peak positions (will use looser criteria here than in id)
										peak_fit_par pfpar = findPeak(sp->sumHist,cntr,5.,20);
										//reportPeakPositions(&pfpar);
										//find matching cascades
										findCascadeFromGammaEInSpec(&pfpar,nd,ge);
									}
								
							}
					}
				free(sp);
			}
		}else if(strTokCmp(cmd,"rank",0)==0){
			strcpy(cmd2,cmd);
			tok=strtok (cmd2," ");
			if((tok = strtok (NULL, " "))==NULL)//read the 2nd entry in the command
				printf("No filename specified.\n");
			else{
				printf("Ranking nuclides using parameters from file: %s.\n",tok);
				nuclide_rank_par *nrp=(nuclide_rank_par*)malloc(sizeof(nuclide_rank_par));
				if(readNuclideRankParameters(tok,nrp)==1){
					rankNuclides(nrp,nd);
				}
				free(nrp);
			}
		}else if(strTokCmp(cmd,"thereisnohope",0)==0){
				printf("Did I tell you the story about the man with a watch in a cave?\n");
		}else if(strcmp(cmd,"")==0){
				//do nothing
		}else{
			printf("Unknown command: %s\nEnter 'help' for a list of commands.\n",cmd);
		}
	}
	

  return 0; //great success
}
