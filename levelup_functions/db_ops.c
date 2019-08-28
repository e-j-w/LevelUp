#include "db_ops.h"

int nameToNuclIndex(const char * name, ndata *nd)
{
	int i;
	for(i=0;i<nd->numNucl;i++)
		if(strcmp(nd->nuclData[i].nuclName,name)==0)
			return i;
			
	return -1;//negative value indicates failure
}

int NZToNuclIndex(int N, int Z, ndata *nd)
{
	int i;
	for(i=0;i<nd->numNucl;i++)
		if(nd->nuclData[i].N==N)
			if(nd->nuclData[i].Z==Z)
				return i;
	
	return -1;//negative value indicates failure
}

int isNuclRadioactive(int nucl, ndata *nd)
{
	int i;
	for(i=0;i<nd->nuclData[nucl].numLevels;i++){
		if(nd->nuclData[nucl].levels[i].energy <= 0.0){
			if(nd->nuclData[nucl].levels[i].lifetimeUnit == -1){
				return 0;
			}else{
				return 1;
			}
		}
	}
	printf("Could not find ground state data for nucleus %s, unknown if it is radioactive.\n",nd->nuclData[nucl].nuclName);
	return 0;
}

//returns:
//0=not on any shell closure
//1=on neutron shell closure
//2=on proton shell closure
//3=on neutron AND proton shell closure
//4=on neuton subshell closure
//5=on proton subshell closure
//6=on neutron AND proton subshell closure
//7=on neutron shell AND proton subshell closure
//8=on neutron subshell AND proton shell closure
int isNuclOnShellClosure(int nucl, ndata *nd)
{
	int protonShell=0;
	int neutronShell=0;

	if((nd->nuclData[nucl].Z == 2)||(nd->nuclData[nucl].Z == 8)||(nd->nuclData[nucl].Z == 20)||(nd->nuclData[nucl].Z == 28)||(nd->nuclData[nucl].Z == 50)||(nd->nuclData[nucl].Z == 82)){
		protonShell = 1;
	}else if((nd->nuclData[nucl].Z == 6)||(nd->nuclData[nucl].Z == 14)||(nd->nuclData[nucl].Z == 16)||(nd->nuclData[nucl].Z == 32)||(nd->nuclData[nucl].Z == 38)||(nd->nuclData[nucl].Z == 40)||(nd->nuclData[nucl].Z == 58)||(nd->nuclData[nucl].Z == 64)||(nd->nuclData[nucl].Z == 76)||(nd->nuclData[nucl].Z == 80)){
		protonShell = 2;
	}
		

	if((nd->nuclData[nucl].N == 2)||(nd->nuclData[nucl].N == 8)||(nd->nuclData[nucl].N == 20)||(nd->nuclData[nucl].N == 28)||(nd->nuclData[nucl].N == 50)||(nd->nuclData[nucl].N == 82)||(nd->nuclData[nucl].N == 126)){
		neutronShell = 1;
	}else if((nd->nuclData[nucl].N == 6)||(nd->nuclData[nucl].N == 14)||(nd->nuclData[nucl].N == 16)||(nd->nuclData[nucl].N == 32)||(nd->nuclData[nucl].N == 38)||(nd->nuclData[nucl].N == 40)||(nd->nuclData[nucl].N == 58)||(nd->nuclData[nucl].N == 64)||(nd->nuclData[nucl].N == 76)||(nd->nuclData[nucl].N == 80)){
		neutronShell = 2;
	}
	
	if((protonShell==0)&&(neutronShell==0))
		return 0;
	if((protonShell==0)&&(neutronShell==1))
		return 1;
	if((protonShell==1)&&(neutronShell==0))
		return 2;
	if((protonShell==1)&&(neutronShell==1))
		return 3;
	if((protonShell==0)&&(neutronShell==2))
		return 4;
	if((protonShell==2)&&(neutronShell==0))
		return 5;
	if((protonShell==2)&&(neutronShell==2))
		return 6;
	if((protonShell==2)&&(neutronShell==1))
		return 7;
	if((protonShell==1)&&(neutronShell==2))
		return 8;
	
	printf("WARNING: unknown shell closure for nucleus %s (%i %i)\n",nd->nuclData[nucl].nuclName,protonShell,neutronShell);
	return 0;
}
void fillShellClosure(int nucl, ndata *nd, char *scstr){

	int scval = isNuclOnShellClosure(nucl,nd);

	strcpy(scstr,""); //clear the string
	if(scval == 0)
		strcpy(scstr,"not on a shell closure");
	else if(scval == 1)
		strcpy(scstr,"on a neutron shell closure");
	else if(scval == 2)
		strcpy(scstr,"on a proton shell closure");
	else if(scval == 3)
		strcpy(scstr,"on a double shell closure");
	else if(scval == 4)
		strcpy(scstr,"on a neutron subshell closure");
	else if(scval == 5)
		strcpy(scstr,"on a proton subshell closure");
	else if(scval == 6)
		strcpy(scstr,"on a double subshell closure");
	else if(scval == 7)
		strcpy(scstr,"on a proton subshell and a neutron shell closure");
	else if(scval == 8)
		strcpy(scstr,"on a proton shell and a neutron subshell closure");

}

void showCascadeData(int nucl, ndata *nd)
{
	//dump cascade data
	int m,n;

	printf("Printing cascade data for the ");
	if(isNuclRadioactive(nucl,nd)){
		printf("radioactive");
	}else{
		printf("stable");
	}
	printf(" nucleus %s.  %i cascades found.\n",nd->nuclData[nucl].nuclName,nd->nuclData[nucl].numCascades);
	if(nd->nuclData[nucl].numCascades>=MAXCASCDESPERNUCL)
		printf("NOTE: cascade listing is truncated due to hitting the maximum size limit.  To increase the size limit, change the value of MAXCASCDESPERNUCL in 'levelup.h'.\n");
	
	if(nucl<MAXNUMNUCL)
		if(nd->nuclData[nucl].numCascades>0)
			{
				for(m=0;m<nd->nuclData[nucl].numCascades;m++)
					{
						printf("CASCADE %i:\nStep   Level Energy   Gamma energy\n   #          (keV)          (keV)\n",m+1);
						for(n=0;n<nd->nuclData[nucl].cascades[m].numLevels;n++)
							printf("%4i     %10.3f     %10.3f\n",n+1,nd->nuclData[nucl].cascades[m].energies[n],nd->nuclData[nucl].cascades[m].gammaEnergies[n]);
					}
				//if(strcmp(nd->nuclData[i].nuclName,"68SE")==0)
					//getc(stdin);
			}
}

//function to find cascade(s) which match the input cascade
void findCascade(gamma_cascade * c, int numToMatch, ndata *nd)
{

	int i,j,k,l,numMatching,numMatched;
	
	numMatched=0;
	for(i=0;i<nd->numNucl;i++)
		for(j=0;j<nd->nuclData[i].numCascades;j++)
			{
				numMatching=0;
				for(l=0;l<c->numLevels;l++)
					for(k=0;k<nd->nuclData[i].cascades[j].numLevels;k++)
						{
							if(fudgeNumbers(nd->nuclData[i].cascades[j].gammaEnergies[k],c->gammaEnergies[l],2.0))
								{
									//if(strcmp(nd->nuclData[i].nuclName,"68SE")==0)
									//	printf("Energy %f matches cascade energy of %f\n",c->gammaEnergies[l],nd->nuclData[i].cascades[j].gammaEnergies[k]);
									numMatching++;
									break;//don't have the same gamma match multiple gammas in a single cascade
								}
						}
				if(numMatching>=numToMatch)
					{
						numMatched++;
						printf("Gamma energies match nucleus: %s (N = %i, Z = %i)\n",nd->nuclData[i].nuclName,nd->nuclData[i].N,nd->nuclData[i].Z);
						//printf("Cascade matches nucleus: %s (N = %i, Z = %i)\n",nd->nuclData[i].nuclName,nd->nuclData[i].N,nd->nuclData[i].Z);
						break; //don't print out the same nucleus more than once
					}
			}
	
	if(numMatched==0)
		printf("Gamma energies didn't match any nuclei in database.\n");
}

//function to find cascade(s) which match the input peak parameters
void findCascadesFromFit(peak_fit_par * par, casc_match_par * cpar, ndata *nd)
{
	int i,j,k,l;
	
	for(i=0;i<nd->numNucl;i++)
		for(j=0;j<nd->nuclData[i].numCascades;j++)
			{
				cpar->numMatchingGammas[i][j]=0;
				for(l=0;l<par->numPeaksFound;l++)
					for(k=0;k<nd->nuclData[i].cascades[j].numLevels;k++)
						{
							if(fudgeNumbers(nd->nuclData[i].cascades[j].gammaEnergies[k],par->centroid[l],2.0))
								{
									//if(strcmp(nd->nuclData[i].nuclName,"68SE")==0)
									//	printf("Energy %f matches cascade energy of %f\n",par->centroid[l],nd->nuclData[i].cascades[j].gammaEnergies[k]);
									cpar->matchCentroid[i][j][cpar->numMatchingGammas[i][j]]=par->centroid[l];
									cpar->numMatchingGammas[i][j]++;
									break;//don't have the same gamma match multiple gammas in a single cascade
								}
						}
			}
}

void findCascadeInSpec(peak_fit_par * par, ndata *nd)
{
	int i,j,k,l;
	int *cascFoundForNucl=(int*)calloc(MAXNUMNUCL,sizeof(int));//initialize to 0 using calloc
	int maxCascSize=par->numPeaksFound;
	
	casc_match_par *cpar=(casc_match_par*)malloc(sizeof(casc_match_par));
	findCascadesFromFit(par,cpar,nd);
	
	printf("NUCLEI INDENTIFIED FROM CASCADES:\n");
	
	int maxNumMatching=0;
	for(i=maxCascSize;i>=3;i--)
		for(j=0;j<nd->numNucl;j++)
			for(k=0;k<nd->nuclData[j].numCascades;k++)
				if(cpar->numMatchingGammas[j][k]==i)
					if(i>maxNumMatching)
						maxNumMatching=i;
	
	if(maxNumMatching>1)
		{
			for(i=maxNumMatching;i>=(maxNumMatching-1);i--)
				for(j=0;j<nd->numNucl;j++)
					if(cascFoundForNucl[j]==0)
						for(k=0;k<nd->nuclData[j].numCascades;k++)
							if(cascFoundForNucl[j]==0)
								if(cpar->numMatchingGammas[j][k]==i)
									{
										printf("%s:",nd->nuclData[j].nuclName);
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

void findCascadeFromGammaEInSpec(peak_fit_par * par, ndata *nd, double ge)
{
	int i,j,k,l,m;
	//int *cascFoundForNucl=(int*)calloc(MAXNUMNUCL,sizeof(int));//initialize to 0 using calloc
	int maxCascSize=par->numPeaksFound;
	
	casc_match_par *cpar=(casc_match_par*)malloc(sizeof(casc_match_par));
	findCascadesFromFit(par,cpar,nd);
	
	printf("CASCADE CANDIDATES IN SPECTRUM MATCHING GAMMA ENERGY: %5.1f keV\n",ge);
	
	int maxNumMatching=0;
	for(i=maxCascSize;i>=3;i--)
		for(j=0;j<nd->numNucl;j++)
			for(k=0;k<nd->nuclData[j].numCascades;k++)
				if(cpar->numMatchingGammas[j][k]==i)
					if(i>maxNumMatching)
						maxNumMatching=i;
	
	int matched=0;
	if(maxNumMatching>1)
		{
			for(i=maxNumMatching;i>=(maxNumMatching-1);i--)
				for(j=0;j<nd->numNucl;j++)
					//if(cascFoundForNucl[j]==0)
						for(k=0;k<nd->nuclData[j].numCascades;k++)
							//if(cascFoundForNucl[j]==0)
								if(cpar->numMatchingGammas[j][k]==i)
									for(l=0;l<i;l++)
										if(fudgeNumbers(ge,cpar->matchCentroid[j][k][l],2.0))
											{
												printf("%s:",nd->nuclData[j].nuclName);
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



void findOverlappingLevels(int nucl1, int nucl2, ndata *nd)
{
	int i,j,k,l;
	int overlapFound=0;
	int nearbyFound=0;
	
	printf("\n");
	for(i=0;i<nd->nuclData[nucl1].numLevels;i++)
		for(j=0;j<nd->nuclData[nucl2].numLevels;j++)
			for(k=0;k<nd->nuclData[nucl1].levels[i].numGammas;k++)
				for(l=0;l<nd->nuclData[nucl1].levels[j].numGammas;l++)
					{
						if(fudgeNumbers(nd->nuclData[nucl1].levels[i].gamma_energies[k],nd->nuclData[nucl2].levels[j].gamma_energies[l],2.0))
							{
								printf("Overlap at energies: %11.1f keV (%s lv#%2i, ",nd->nuclData[nucl1].levels[i].gamma_energies[k],nd->nuclData[nucl1].nuclName,i+1);
								if(nd->nuclData[nucl1].levels[i].gamma_intensities[k]>0)
									printf("I=%4.0f",nd->nuclData[nucl1].levels[i].gamma_intensities[k]);
								else
									printf("I=unwn");
								printf("), %6.1f keV (%s lv#%2i, ",nd->nuclData[nucl2].levels[j].gamma_energies[l],nd->nuclData[nucl2].nuclName,j+1);
								if(nd->nuclData[nucl2].levels[j].gamma_intensities[l]>0)
									printf("I=%4.0f)\n",nd->nuclData[nucl2].levels[j].gamma_intensities[l]);
								else
									printf("I=unwn)\n");
								overlapFound++;
							}
					}
	printf("\n");
	for(i=0;i<nd->nuclData[nucl1].numLevels;i++)
		for(j=0;j<nd->nuclData[nucl2].numLevels;j++)
			for(k=0;k<nd->nuclData[nucl1].levels[i].numGammas;k++)
				for(l=0;l<nd->nuclData[nucl1].levels[j].numGammas;l++)
					{
						if(fudgeNumbers(nd->nuclData[nucl1].levels[i].gamma_energies[k],nd->nuclData[nucl2].levels[j].gamma_energies[l],10.0))
							if(!(fudgeNumbers(nd->nuclData[nucl1].levels[i].gamma_energies[k],nd->nuclData[nucl2].levels[j].gamma_energies[l],2.0)))
								{					
									printf("Nearby peaks at energies: %6.1f keV (%s lv#%2i, ",nd->nuclData[nucl1].levels[i].gamma_energies[k],nd->nuclData[nucl1].nuclName,i+1);
									if(nd->nuclData[nucl1].levels[i].gamma_intensities[k]>0)
										printf("I=%4.0f",nd->nuclData[nucl1].levels[i].gamma_intensities[k]);
									else
										printf("I=unwn");
									printf("), %6.1f keV (%s lv#%2i, ",nd->nuclData[nucl2].levels[j].gamma_energies[l],nd->nuclData[nucl2].nuclName,j+1);
									if(nd->nuclData[nucl2].levels[j].gamma_intensities[l]>0)
										printf("I=%4.0f)\n",nd->nuclData[nucl2].levels[j].gamma_intensities[l]);
									else
										printf("I=unwn)\n");
									nearbyFound++;
								}
					}
					
	if((overlapFound+nearbyFound)==0)
		printf("\nNo overlap in gamma ray energies found in the nuclei %s and %s.\n",nd->nuclData[nucl1].nuclName,nd->nuclData[nucl2].nuclName);
	else
		printf("\n%i overlapping and %i nearby gamma ray energies found in the nuclei %s and %s.\n",overlapFound,nearbyFound,nd->nuclData[nucl1].nuclName,nd->nuclData[nucl2].nuclName);

}


//finds levels of nucleus 1 which overlaps with levels of nuclei in the
//region of nucleus 2 (N,Z +/- searchDim)
void findOverlappingLevelsInRegion(int nucl1, int nucl2, int searchDim, ndata *nd)
{
	int i,j,n2ind;
	
	for(i=nd->nuclData[nucl2].N-searchDim;i<=nd->nuclData[nucl2].N+searchDim;i++)
		for(j=nd->nuclData[nucl2].Z-searchDim;j<=nd->nuclData[nucl2].Z+searchDim;j++)
			{
				n2ind=NZToNuclIndex(i,j,nd);
				if(n2ind>=0)
					if(nucl1!=n2ind)
						findOverlappingLevels(nucl1,n2ind,nd);
			}
}


//finds levels matching the specified energy in the region of the specified nucleus
//(region size is N,Z +/- searchDim)
void findLevelInRegion(double energy, int nucl, int searchDim, ndata *nd)
{
	int i,j,k,l,nind;
	
	for(i=nd->nuclData[nucl].N-searchDim;i<=nd->nuclData[nucl].N+searchDim;i++)
		for(j=nd->nuclData[nucl].Z-searchDim;j<=nd->nuclData[nucl].Z+searchDim;j++)
			{
				nind=NZToNuclIndex(i,j,nd);
				if(nind>=0)
					for(k=0;k<nd->nuclData[nind].numLevels;k++)
						for(l=0;l<nd->nuclData[nind].levels[k].numGammas;l++)
							if(fudgeNumbers(nd->nuclData[nind].levels[k].gamma_energies[l],energy,2.0))
								printf("%s: contains matching gamma with energy %6.1f keV (lv#%2i).\n",nd->nuclData[nind].nuclName,nd->nuclData[nind].levels[k].gamma_energies[l],k+1);
			}
	printf("\n");
	for(i=nd->nuclData[nucl].N-searchDim;i<=nd->nuclData[nucl].N+searchDim;i++)
		for(j=nd->nuclData[nucl].Z-searchDim;j<=nd->nuclData[nucl].Z+searchDim;j++)
			{
				nind=NZToNuclIndex(i,j,nd);
				if(nind>=0)
					for(k=0;k<nd->nuclData[nind].numLevels;k++)
						for(l=0;l<nd->nuclData[nind].levels[k].numGammas;l++)
							if(fudgeNumbers(nd->nuclData[nind].levels[k].gamma_energies[l],energy,10.0))
								if(!(fudgeNumbers(nd->nuclData[nind].levels[k].gamma_energies[l],energy,2.0)))
								printf("%s: contains nearby gamma with energy %6.1f keV (lv#%2i).\n",nd->nuclData[nind].nuclName,nd->nuclData[nind].levels[k].gamma_energies[l],k+1);
			}
}

void fillSpinPar(level *gl, char *spstr){

	int i;
	char val[256];

	strcpy(spstr,""); //clear the string

	for(i=0;i<gl->numspinparvals;i++){
		if((gl->spval[i].tentative == 1)||(gl->spval[i].tentative == 2)){
			if((i==0)||((i>0)&&((gl->spval[i-1].tentative != 1)&&(gl->spval[i-1].tentative != 2)))){
				strcat(spstr,"(");
			}
		}
		if(gl->spval[i].halfInt == 1){
			sprintf(val,"%i/2",gl->spval[i].spinVal);
		}else{
			sprintf(val,"%i",gl->spval[i].spinVal);
		}
		strcat(spstr,val);
		if(gl->spval[i].parVal == -1){
			strcat(spstr,"-");
		}else if(gl->spval[i].parVal == 1){
			strcat(spstr,"+");
		}
		if((gl->spval[i].tentative == 1)||(gl->spval[i].tentative == 2)){
			if((i==gl->numspinparvals-1)||((i<gl->numspinparvals-1)&&((gl->spval[i+1].tentative != 1)&&(gl->spval[i+1].tentative != 2)))){
				strcat(spstr,")");
			}
		}
		if(i!=gl->numspinparvals-1){
			strcat(spstr,",");
		}
	}
}

void fillLifetime(level *gl, char *ltstr){

	char val[256];

	strcpy(ltstr,""); //clear the string

	if(gl->lifetime >=0){
		if(gl->lifetimeUnit < 0){
			strcat(ltstr,"STABLE");
		}else{
			sprintf(val,"%0.3f ",gl->lifetime);
			strcat(ltstr,val);
			if(gl->lifetimeUnit == 0){
				strcat(ltstr,"years");
			}else if(gl->lifetimeUnit == 1){
				strcat(ltstr,"days");
			}else if(gl->lifetimeUnit == 2){
				strcat(ltstr,"hours");
			}else if(gl->lifetimeUnit == 3){
				strcat(ltstr,"minutes");
			}else if(gl->lifetimeUnit == 4){
				strcat(ltstr,"s");
			}else if(gl->lifetimeUnit == 5){
				strcat(ltstr,"ms");
			}else if(gl->lifetimeUnit == 6){
				strcat(ltstr,"us");
			}else if(gl->lifetimeUnit == 7){
				strcat(ltstr,"ns");
			}else if(gl->lifetimeUnit == 8){
				strcat(ltstr,"ps");
			}else if(gl->lifetimeUnit == 9){
				strcat(ltstr,"fs");
			}else if(gl->lifetimeUnit == 10){
				strcat(ltstr,"as");
			}else if(gl->lifetimeUnit == 11){
				strcat(ltstr,"eV");
			}else if(gl->lifetimeUnit == 12){
				strcat(ltstr,"keV");
			}else if(gl->lifetimeUnit == 13){
				strcat(ltstr,"MeV");
			}else{
				strcat(ltstr,"??");
			}
		}
	}

}

//shows a table of levels
//numLevels is the number of levels to show (0 -> show all levels)
void showLevelData(int nucl, ndata *nd, int numLevels)
{
	int i,j;
	double finalEnergy;
	char spstr[256],ltstr[256];
	
	printf("Printing level data for the ");
	if(isNuclRadioactive(nucl,nd)){
		printf("radioactive");
	}else{
		printf("stable");
	}
	printf(" nucleus %s.  ",nd->nuclData[nucl].nuclName);
	if(numLevels<=0)
		printf("%i levels found.\n",nd->nuclData[nucl].numLevels);
	else
		printf("Showing the first %i levels.\n",numLevels);
	
	if((nd->nuclData[nucl].numLevels>=MAXLEVELSPERNUCL)||(numLevels>=MAXLEVELSPERNUCL))
		printf("NOTE: level listing is truncated due to hitting the maximum size limit.  To increase the size limit, change the value of MAXLEVELSPERNUCL in 'levelup.h'.\n");

	printf("\n Level Energy  Gamma Energy  Intensity    Final E              Jpi       Lifetime\n        (keV)         (keV)                 (keV)\n\n");
	if(nucl<MAXNUMNUCL)
		for(i=0;i<nd->nuclData[nucl].numLevels;i++)
			if((numLevels<=0)||(i<numLevels))
				{
					if(nd->nuclData[nucl].levels[i].numGammas>0)
						{
							finalEnergy=nd->nuclData[nucl].levels[i].energy-nd->nuclData[nucl].levels[i].gamma_energies[0];
							if(finalEnergy<0)  
								finalEnergy=0.;
							fillSpinPar(&nd->nuclData[nucl].levels[i],spstr);
							fillLifetime(&nd->nuclData[nucl].levels[i],ltstr);
							if(nd->nuclData[nucl].levels[i].gamma_intensities[0]>0)
								printf("%13.3f %13.3f %10.1f %10.3f %16s %14s\n",nd->nuclData[nucl].levels[i].energy,nd->nuclData[nucl].levels[i].gamma_energies[0],nd->nuclData[nucl].levels[i].gamma_intensities[0],finalEnergy,spstr,ltstr);
							else
								printf("%13.3f %13.3f    unknown %10.3f %16s %14s\n",nd->nuclData[nucl].levels[i].energy,nd->nuclData[nucl].levels[i].gamma_energies[0],finalEnergy,spstr,ltstr);
						}
					else
						{
							fillSpinPar(&nd->nuclData[nucl].levels[i],spstr);
							fillLifetime(&nd->nuclData[nucl].levels[i],ltstr);
							printf("%13.3f                                     %16s %14s\n",nd->nuclData[nucl].levels[i].energy,spstr,ltstr);
						}
					for(j=1;j<nd->nuclData[nucl].levels[i].numGammas;j++)
						{
							finalEnergy=nd->nuclData[nucl].levels[i].energy-nd->nuclData[nucl].levels[i].gamma_energies[j];
							if(finalEnergy<0)  finalEnergy=0.;
							if(nd->nuclData[nucl].levels[i].gamma_intensities[j]>0)
								printf("              %13.3f %10.1f %10.3f\n",nd->nuclData[nucl].levels[i].gamma_energies[j],nd->nuclData[nucl].levels[i].gamma_intensities[j],finalEnergy);
							else
								printf("              %13.3f    unknown %10.3f\n",nd->nuclData[nucl].levels[i].gamma_energies[j],finalEnergy);
						}
				}
}

void showNZ(int nucl, ndata *nd)
{
	char scstr[256];
	fillShellClosure(nucl,nd,scstr);
	printf("For nucleus %s: N = %i, Z = %i\n",nd->nuclData[nucl].nuclName,nd->nuclData[nucl].N,nd->nuclData[nucl].Z);
	printf("This nucleus is %s.\n",scstr);
}

void showNuclNames(ndata *nd)
{
	int i;
	printf("Known nuclei: %s",nd->nuclData[0].nuclName);
	for(i=1;i<nd->numNucl;i++)
		{
			printf(", %s",nd->nuclData[i].nuclName);
		}
	printf("\n");
}

void rebuildDatabase(ndata *nd, FILE *db)
{
	int i;
	char fileName[256],str[8];
	
	initialize_database(nd);
	
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
			readENSDFFile(fileName,nd); //grab data from the ENSDF file (see parse_ENSDF.c)
		}
	printf("Data imported for %i nuclei.\n",nd->numNucl);
	printf("Generating cascade data.\n");
	generateCascadeData(nd);
	
	//write the database to disk
	if(nd->numNucl<=0)
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
	fwrite(nd,sizeof(ndata),1,db);
	fclose(db);
	printf("Database rebuild finished.\n");
}
