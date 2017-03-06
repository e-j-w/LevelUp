//definitions
#include "levelup.h"
//functions and logic
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
  gdata *gd=(gdata*)malloc(sizeof(gdata));
  initialize_database(gd);
	printf("Data structures allocated.\n");
	
	
	//process the ENSDF database file
	char fileName[256];
	int i=0;
	int j=0;
	strcpy(fileName,"");
	strcat(fileName,getenv("ENSDF"));
	strcat(fileName,"ensdf_db");
	if((db=fopen(fileName,"r"))==NULL)
		{
			printf("ENSDF database file not found, rebuilding database...\n");
			rebuildDatabase(gd,db);
			db=fopen(fileName,"r");//reopen file for reading
		}
	while(i!=1)
		{	
			i=fread(gd,sizeof(gdata),1,db);
			if(i==1)
				{
					printf("ENSDF database retrieved from disk.  Data for %i nuclei found.\n",gd->numNucl);
				}
			else
				{
					j++;
					printf("ERROR: ENSDF database could not be read from file '%s'.  Rebuilding (attempt %i)...\n",fileName,j);
					rebuildDatabase(gd,db);
					db=fopen(fileName,"r");//reopen file for reading
					
				}
			if(j>3)
				{
					printf("ERROR: ENSDF database could not be built.  Exiting...\n");
					exit(-1);
				}
	
		}
	fclose(db);
	
	
	//process user commands
	printf("\nEnter command.  Enter 'help' for a list of commands.\n");
	char cmd[256],cmd2[256];
	char *tok;
	while(1)
		{
			printf("LevelUp > ");
  		fgets(cmd,256,stdin);
  		cmd[strcspn(cmd, "\r\n")] = 0;//strips newline characters from the string read by fgets
  		if(strcmp(cmd,"help")==0)
  			{
  				printf("Command list:\n");
  				printf("  casc NUCL - prints cascade data for the specified nucleus\n");
  				printf("              (NUCL is the nucleus name, eg. '68SE')\n");
  				printf("  lev NUCL - prints level data for the specified nucleus\n");
  				printf("  ol NUCL1 NUCL2 - finds overlapping gamma rays in the two\n"); 
  				printf("                   specified nuclei\n");
  				printf("  nz NUCL - show N, Z numbers for the specified nucleus\n");
  				printf("  listnuc, ln - list names of nuclei in the ENSDF database\n");
  				printf("  findcasc, fc - find nuclei which match a cascade that you enter\n");
  				printf("  rebuild - rebuild the ENSDF database from ENSDF files on disk\n");
  				printf("  help - list commands\n");
  				printf("  exit, quit - exit the program\n");
  			}
  		else if((strcmp(cmd,"exit")==0)||(strcmp(cmd,"quit")==0))
  			{
  				printf("Exiting...\n");
  				exit(0);
  			}
  		else if(strTokCmp(cmd,"casc",0)==0)
  			{
  				strcpy(cmd2,cmd);
  				tok=strtok (cmd2," ");
					tok = strtok (NULL, " ");//read the 2nd entry in the command
					int nucl=nameToNuclIndex(tok,gd);
					if ((tok != NULL) && (nucl>=0))
						{
  						showCascadeData(nucl,gd);
  					}
  				else
  					printf("Unknown nucleus: %s\n",tok);
  			}
  		else if(strTokCmp(cmd,"lev",0)==0)
  			{
  				strcpy(cmd2,cmd);
  				tok=strtok (cmd2," ");
					tok = strtok (NULL, " ");//read the 2nd entry in the command
					int nucl=nameToNuclIndex(tok,gd);
					if ((tok != NULL) && (nucl>=0))
						{
  						showLevelData(nucl,gd);
  					}
  				else
  					printf("Unknown nucleus: %s\n",tok);
  			}
  		else if((strTokCmp(cmd,"overlap",0)==0)||(strTokCmp(cmd,"ol",0)==0))
  			{
  				int findOL=1;
  				strcpy(cmd2,cmd);
  				tok=strtok (cmd2," ");
					tok = strtok (NULL, " ");//read the 2nd entry in the command
					int nucl1=nameToNuclIndex(tok,gd);
					if ((tok == NULL) || (nucl1<0))
						{
							printf("Unknown nucleus: %s\n",tok);
							findOL=0;
						}
					if (findOL)
						{
							tok = strtok (NULL, " ");//read the 3rd entry in the command
							int nucl2=nameToNuclIndex(tok,gd);
							if ((tok == NULL) || (nucl2<0))
								{
									printf("Unknown nucleus: %s\n",tok);
									findOL=0;
								}
							if (findOL)
								{
									findOverlappingLevels(nucl1,nucl2,gd);
								}
						}
					
  			}
  		else if(strTokCmp(cmd,"nz",0)==0)
  			{
  				strcpy(cmd2,cmd);
  				tok=strtok (cmd2," ");
					tok = strtok (NULL, " ");//read the 2nd entry in the command
					int nucl=nameToNuclIndex(tok,gd);
					if ((tok != NULL) && (nucl>=0))
						{
  						showNZ(nucl,gd);
  					}
  				else
  					printf("Unknown nucleus: %s\n",tok);
  			}
  		else if((strcmp(cmd,"listnuc")==0)||(strcmp(cmd,"ln")==0))
  			{
  				showNuclNames(gd);
  			}
  		else if((strcmp(cmd,"findcasc")==0)||(strcmp(cmd,"fc")==0))
  			{
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
  					findCascade(c,c->numLevels,gd);
  				free(c);
  			}
  		else if(strcmp(cmd,"rebuild")==0)
  			{
  				printf("Rebuilding ENSDF database...\n");
					rebuildDatabase(gd,db);
					db=fopen(fileName,"r");//reopen file for reading
					i=fread(gd,sizeof(gdata),1,db);
					while(i!=1)
						{	
							i=fread(gd,sizeof(gdata),1,db);
							if(i==1)
								{
									printf("ENSDF database retrieved from disk.  Data for %i nuclei found.\n",gd->numNucl);
								}
							else
								{
									j++;
									printf("ERROR: ENSDF database could not be read from file '%s'.  Rebuilding (attempt %i)...\n",fileName,j);
									rebuildDatabase(gd,db);
									db=fopen(fileName,"r");//reopen file for reading
					
								}
							if(j>3)
								{
									printf("ERROR: ENSDF database could not be built.  Exiting...\n");
									exit(-1);
								}
						}
					fclose(db);
  			}
  		else if(strcmp(cmd,"")==0)
  			{
  				//do nothing
  			}
  		else
  			{
  				printf("Unknown command: %s\nEnter 'help' for a list of commands.\n",cmd);
  			}
		}
	

  return 0; //great success
}
