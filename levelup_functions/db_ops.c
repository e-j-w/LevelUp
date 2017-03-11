#include "db_ops.h"

int nameToNuclIndex(const char * name, gdata *gd)
{
	int i;
	for(i=0;i<gd->numNucl;i++)
		if(strcmp(gd->nuclData[i].nuclName,name)==0)
			return i;
			
	return -1;//negative value indicates failure
}

int NZToNuclIndex(int N, int Z, gdata *gd)
{
	int i;
		for(i=0;i<gd->numNucl;i++)
			if(gd->nuclData[i].N==N)
				if(gd->nuclData[i].Z==Z)
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
						printf("Gamma energies match nucleus: %s (N = %i, Z = %i)\n",gd->nuclData[i].nuclName,gd->nuclData[i].N,gd->nuclData[i].Z);
						//printf("Cascade matches nucleus: %s (N = %i, Z = %i)\n",gd->nuclData[i].nuclName,gd->nuclData[i].N,gd->nuclData[i].Z);
						break; //don't print out the same nucleus more than once
					}
			}
	
	if(numMatched==0)
		printf("Gamma energies didn't match any nuclei in database.\n");
}

//function to find cascade(s) which match the input peak parameters
void findCascadesFromFit(peak_fit_par * par, casc_match_par * cpar, gdata *gd)
{
	int i,j,k,l;
	
	for(i=0;i<gd->numNucl;i++)
		for(j=0;j<gd->nuclData[i].numCascades;j++)
			{
				cpar->numMatchingGammas[i][j]=0;
				for(l=0;l<par->numPeaksFound;l++)
					for(k=0;k<gd->nuclData[i].cascades[j].numLevels;k++)
						{
							if(fudgeNumbers(gd->nuclData[i].cascades[j].gammaEnergies[k],par->centroid[l],2.0))
								{
									//if(strcmp(gd->nuclData[i].nuclName,"68SE")==0)
									//	printf("Energy %f matches cascade energy of %f\n",par->centroid[l],gd->nuclData[i].cascades[j].gammaEnergies[k]);
									cpar->matchCentroid[i][j][cpar->numMatchingGammas[i][j]]=par->centroid[l];
									cpar->numMatchingGammas[i][j]++;
									break;//don't have the same gamma match multiple gammas in a single cascade
								}
						}
			}
}

void findCascadeInSpec(peak_fit_par * par, gdata *gd)
{
	int i,j,k,l;
	int *cascFoundForNucl=(int*)calloc(MAXNUMNUCL,sizeof(int));//initialize to 0 using calloc
	int maxCascSize=par->numPeaksFound;
	
	casc_match_par *cpar=(casc_match_par*)malloc(sizeof(casc_match_par));
	findCascadesFromFit(par,cpar,gd);
	
	printf("NUCLEI INDENTIFIED FROM CASCADES:\n");
	
	int maxNumMatching=0;
	for(i=maxCascSize;i>=3;i--)
		for(j=0;j<gd->numNucl;j++)
			for(k=0;k<gd->nuclData[j].numCascades;k++)
				if(cpar->numMatchingGammas[j][k]==i)
					if(i>maxNumMatching)
						maxNumMatching=i;
	
	if(maxNumMatching>1)
		{
			for(i=maxNumMatching;i>=(maxNumMatching-1);i--)
				for(j=0;j<gd->numNucl;j++)
					if(cascFoundForNucl[j]==0)
						for(k=0;k<gd->nuclData[j].numCascades;k++)
							if(cascFoundForNucl[j]==0)
								if(cpar->numMatchingGammas[j][k]==i)
									{
										printf("%s:",gd->nuclData[j].nuclName);
										for(l=0;l<i;l++)
											if(l==0)
												printf(" %5.1f",cpar->matchCentroid[j][k][l]);
											else
												printf(", %5.1f",cpar->matchCentroid[j][k][l]);
										printf(" keV (matched %i gamma rays)\n",i);
										cascFoundForNucl[j]=1;
									}
		}
	else
		{
			printf("No matching nuclei found.\n");
		}
	free(cascFoundForNucl);
	free(cpar);
}

void findCascadeFromGammaEInSpec(peak_fit_par * par, gdata *gd, double ge)
{
	int i,j,k,l,m;
	//int *cascFoundForNucl=(int*)calloc(MAXNUMNUCL,sizeof(int));//initialize to 0 using calloc
	int maxCascSize=par->numPeaksFound;
	
	casc_match_par *cpar=(casc_match_par*)malloc(sizeof(casc_match_par));
	findCascadesFromFit(par,cpar,gd);
	
	printf("CASCADE CANDIDATES IN SPECTRUM MATCHING GAMMA ENERGY: %5.1f keV\n",ge);
	
	int maxNumMatching=0;
	for(i=maxCascSize;i>=3;i--)
		for(j=0;j<gd->numNucl;j++)
			for(k=0;k<gd->nuclData[j].numCascades;k++)
				if(cpar->numMatchingGammas[j][k]==i)
					if(i>maxNumMatching)
						maxNumMatching=i;
	
	int matched=0;
	if(maxNumMatching>1)
		{
			for(i=maxNumMatching;i>=(maxNumMatching-1);i--)
				for(j=0;j<gd->numNucl;j++)
					//if(cascFoundForNucl[j]==0)
						for(k=0;k<gd->nuclData[j].numCascades;k++)
							//if(cascFoundForNucl[j]==0)
								if(cpar->numMatchingGammas[j][k]==i)
									for(l=0;l<i;l++)
										if(fudgeNumbers(ge,cpar->matchCentroid[j][k][l],2.0))
											{
												printf("%s:",gd->nuclData[j].nuclName);
												for(m=0;m<i;m++)
													if(m==0)
														printf(" %5.1f",cpar->matchCentroid[j][k][m]);
													else
														printf(", %5.1f",cpar->matchCentroid[j][k][m]);
												printf(" keV (matched %i gamma rays)\n",i);
												matched=1;
												//cascFoundForNucl[j]=1;
												break;
											}
		}
	if(matched==0)
		{
			printf("No matching cascades found.\n");
		}
	//free(cascFoundForNucl);
	free(cpar);
}



void findOverlappingLevels(int nucl1, int nucl2, gdata *gd)
{
	int i,j,k,l;
	int overlapFound=0;
	int nearbyFound=0;
	
	printf("\n");
	for(i=0;i<gd->nuclData[nucl1].numLevels;i++)
		for(j=0;j<gd->nuclData[nucl2].numLevels;j++)
			for(k=0;k<gd->nuclData[nucl1].levels[i].numGammas;k++)
				for(l=0;l<gd->nuclData[nucl1].levels[j].numGammas;l++)
					{
						if(fudgeNumbers(gd->nuclData[nucl1].levels[i].gamma_energies[k],gd->nuclData[nucl2].levels[j].gamma_energies[l],2.0))
							{
								printf("Overlap at energies: %11.1f keV (%s lv#%2i), %6.1f keV (%s lv#%2i)\n",gd->nuclData[nucl1].levels[i].gamma_energies[k],gd->nuclData[nucl1].nuclName,i+1,gd->nuclData[nucl2].levels[j].gamma_energies[l],gd->nuclData[nucl2].nuclName,j+1);
								overlapFound++;
							}
					}
	printf("\n");
	for(i=0;i<gd->nuclData[nucl1].numLevels;i++)
		for(j=0;j<gd->nuclData[nucl2].numLevels;j++)
			for(k=0;k<gd->nuclData[nucl1].levels[i].numGammas;k++)
				for(l=0;l<gd->nuclData[nucl1].levels[j].numGammas;l++)
					{
						if(fudgeNumbers(gd->nuclData[nucl1].levels[i].gamma_energies[k],gd->nuclData[nucl2].levels[j].gamma_energies[l],10.0))
							if(!(fudgeNumbers(gd->nuclData[nucl1].levels[i].gamma_energies[k],gd->nuclData[nucl2].levels[j].gamma_energies[l],2.0)))
								{
									printf("Nearby peaks at energies: %6.1f keV (%s lv#%2i), %6.1f keV (%s lv#%2i)\n",gd->nuclData[nucl1].levels[i].gamma_energies[k],gd->nuclData[nucl1].nuclName,i+1,gd->nuclData[nucl2].levels[j].gamma_energies[l],gd->nuclData[nucl2].nuclName,j+1);
									nearbyFound++;
								}
					}
					
	if((overlapFound+nearbyFound)==0)
		printf("\nNo overlap in gamma ray energies found in the nuclei %s and %s.\n",gd->nuclData[nucl1].nuclName,gd->nuclData[nucl2].nuclName);
	else
		printf("\n%i overlapping and %i nearby gamma ray energies found in the nuclei %s and %s.\n",overlapFound,nearbyFound,gd->nuclData[nucl1].nuclName,gd->nuclData[nucl2].nuclName);

}


//finds levels of nucleus 1 which overlaps with levels of nuclei in the
//region of nucleus 2 (N,Z +/- searchDim)
void findOverlappingLevelsInRegion(int nucl1, int nucl2, int searchDim, gdata *gd)
{
	int i,j,n2ind;
	
	for(i=gd->nuclData[nucl2].N-searchDim;i<=gd->nuclData[nucl2].N+searchDim;i++)
		for(j=gd->nuclData[nucl2].Z-searchDim;j<=gd->nuclData[nucl2].Z+searchDim;j++)
			{
				n2ind=NZToNuclIndex(i,j,gd);
				if(n2ind>=0)
					if(nucl1!=n2ind)
						findOverlappingLevels(nucl1,n2ind,gd);
			}
}

//shows a table of levels
//numLevels is the number of levels to show (0 -> show all levels)
void showLevelData(int nucl, gdata *gd, int numLevels)
{
	int i,j;
	double finalEnergy;
	
	if(numLevels<0)
		printf("Printing level data for nucleus %s.  %i levels found.\n",gd->nuclData[nucl].nuclName,gd->nuclData[nucl].numLevels);
	else
		printf("Printing level data for nucleus %s.  Showing the first %i levels.\n",gd->nuclData[nucl].nuclName,numLevels);
	
	if((gd->nuclData[nucl].numLevels>=MAXLEVELSPERNUCL)||(numLevels>=MAXLEVELSPERNUCL))
		printf("NOTE: level listing is truncated due to hitting the maximum size limit.  To increase the size limit, change the value of MAXLEVELSPERNUCL in 'levelup.h'.\n");

	printf("Level Energy    Gamma Energy    Final energy\n       (keV)           (keV)           (keV)\n\n");
	if(nucl<MAXNUMNUCL)
		for(i=0;i<gd->nuclData[nucl].numLevels;i++)
			if((numLevels<=0)||(i<numLevels))
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

void showNZ(int nucl, gdata *gd)
{
	printf("For nucleus %s: N = %i, Z = %i\n",gd->nuclData[nucl].nuclName,gd->nuclData[nucl].N,gd->nuclData[nucl].Z);
}

void showNuclNames(gdata *gd)
{
	int i;
	printf("Known nuclei: %s",gd->nuclData[0].nuclName);
	for(i=1;i<gd->numNucl;i++)
		{
			printf(", %s",gd->nuclData[i].nuclName);
		}
	printf("\n");
}

void rebuildDatabase(gdata *gd, FILE *db)
{
	int i;
	char fileName[256],str[8];
	
	initialize_database(gd);
	
	for(i=1;i<350;i++)
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
