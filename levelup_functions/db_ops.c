#include "db_ops.h"

int nameToNuclIndex(const char * name, gdata *gd)
{
	int i;
	for(i=0;i<gd->numNucl;i++)
		if(strcmp(gd->nuclData[i].nuclName,name)==0)
			return i;
			
	return -1;//negative value indicates failure
}

void showCascadeData(int nucl, gdata *gd)
{
	//dump cascade data
	int m,n;
	
	printf("Printing cascade data for nucleus %s.  %i cascades found.\n",gd->nuclData[nucl].nuclName,gd->nuclData[nucl].numCascades);
	if(gd->nuclData[nucl].numCascades>=MAXCASCDESPERNUCL)
		printf("NOTE: cascade listing is truncated due to hitting the maximum size limit.  To increase the size limit, change the value of MAXCASCDESPERNUCL in 'levelup.h'.\n");
	
	if(nucl<MAXNUMNUCL)
		if(gd->nuclData[nucl].numCascades>0)
			{
				for(m=0;m<gd->nuclData[nucl].numCascades;m++)
					{
						printf("CASCADE %i:\nStep   Level Energy   Gamma energy\n   #          (keV)          (keV)\n",m+1);
						for(n=0;n<gd->nuclData[nucl].cascades[m].numLevels;n++)
							printf("%4i     %10.3f     %10.3f\n",n+1,gd->nuclData[nucl].cascades[m].energies[n],gd->nuclData[nucl].cascades[m].gammaEnergies[n]);
					}
				//if(strcmp(gd->nuclData[i].nuclName,"68SE")==0)
					//getc(stdin);
			}
}

//function to find cascade(s) which match the input cascade
void findCascade(gamma_cascade * c, int numToMatch, gdata *gd)
{

	int i,j,k,l,numMatching,numMatched;
	
	numMatched=0;
	for(i=0;i<gd->numNucl;i++)
		for(j=0;j<gd->nuclData[i].numCascades;j++)
			{
				numMatching=0;
				for(l=0;l<c->numLevels;l++)
					for(k=0;k<gd->nuclData[i].cascades[j].numLevels;k++)
						{
							if(fudgeNumbers(gd->nuclData[i].cascades[j].gammaEnergies[k],c->gammaEnergies[l],2.0))
								{
									//if(strcmp(gd->nuclData[i].nuclName,"68SE")==0)
									//	printf("Energy %f matches cascade energy of %f\n",c->gammaEnergies[l],gd->nuclData[i].cascades[j].gammaEnergies[k]);
									numMatching++;
									break;//don't have the same gamma match multiple gammas in a single cascade
								}
						}
				if(numMatching>=numToMatch)
					{
						numMatched++;
						printf("Cascade matches nucleus: %s\n",gd->nuclData[i].nuclName);
						//printf("Cascade matches nucleus: %s (N = %i, Z = %i)\n",gd->nuclData[i].nuclName,gd->nuclData[i].N,gd->nuclData[i].Z);
						break; //don't print out the same nucleus more than once
					}
			}
	
	if(numMatched==0)
		printf("Cascade didn't match any nuclei in database.\n");
}

void showLevelData(int nucl, gdata *gd)
{
	int i,j;
	double finalEnergy;
	
	printf("Printing level data for nucleus %s.  %i levels found.\n",gd->nuclData[nucl].nuclName,gd->nuclData[nucl].numLevels);
	if(gd->nuclData[nucl].numLevels>=MAXLEVELSPERNUCL)
		printf("NOTE: level listing is truncated due to hitting the maximum size limit.  To increase the size limit, change the value of MAXLEVELSPERNUCL in 'levelup.h'.\n");

	printf("Level Energy    Gamma Energy    Final energy\n       (keV)           (keV)           (keV)\n\n");
	if(nucl<MAXNUMNUCL)
		for(i=0;i<gd->nuclData[nucl].numLevels;i++)
			{
				if(gd->nuclData[nucl].levels[i].numGammas>0)
					{
						finalEnergy=gd->nuclData[nucl].levels[i].energy-gd->nuclData[nucl].levels[i].gamma_energies[0];
						if(finalEnergy<0)  finalEnergy=0.;
						printf("  %10.3f      %10.3f      %10.3f\n",gd->nuclData[nucl].levels[i].energy,gd->nuclData[nucl].levels[i].gamma_energies[0],finalEnergy);
					}
				else
					{
						printf("  %10.3f\n",gd->nuclData[nucl].levels[i].energy);
					}
				for(j=1;j<gd->nuclData[nucl].levels[i].numGammas;j++)
					{
						finalEnergy=gd->nuclData[nucl].levels[i].energy-gd->nuclData[nucl].levels[i].gamma_energies[j];
						if(finalEnergy<0)  finalEnergy=0.;
						printf("                  %10.3f      %10.3f\n",gd->nuclData[nucl].levels[i].gamma_energies[j],finalEnergy);
					}
			}
}

void showNuclNames(gdata *gd)
{
	int i;
	printf("Known nuclei: %s",gd->nuclData[0].nuclName);
	for(i=1;i<gd->numNucl;i++)
		printf(", %s",gd->nuclData[i].nuclName);
		printf("\n");
}

void rebuildDatabase(gdata *gd, FILE *db)
{
	int i;
	char fileName[256],str[8];
	
	initialize_database(gd);
	
	for(i=1;i<200;i++)
		{
			strcpy(fileName,"");
			strcat(fileName,getenv("ENSDF"));
			if(i<10)
				strcat(fileName,"ensdf.00");
			else if(i<100)
				strcat(fileName,"ensdf.0");
			else
				strcat(fileName,"ensdf.");
			sprintf(str,"%i",i);
			strcat(fileName,str);
			readENSDFFile(fileName,gd); //grab data from the ENSDF file (see parse_ENSDF.c)
		}
	printf("Data imported for %i nuclei.\n",gd->numNucl);
	printf("Generating cascade data.\n");
	generateCascadeData(gd);
	
	//write the database to disk
	if(gd->numNucl<=0)
		{
			printf("ERROR: no valid ENSDF data was found.\nPlease check that ENSDF files exist in the directory under the ENDSF environment variable.\n");
			exit(-1);
		}
	printf("Database generated.  Writing to disk...\n");
	strcpy(fileName,"");
	strcat(fileName,getenv("ENSDF"));
	strcat(fileName,"ensdf_db");
	if((db=fopen(fileName,"w"))==NULL)
		{
			printf("ERROR: cannot open the output file '%s'.  Exiting...\n",fileName);
		}
	fwrite(gd,sizeof(gdata),1,db);
	fclose(db);
	printf("Database rebuild finished.\n");
}
