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
	if(nucl<MAXNUMNUCL)
		if(gd->nuclData[nucl].numCascades>0)
			{
				for(m=0;m<gd->nuclData[nucl].numCascades;m++)
					{
						printf("CASCADE %i:\nStep   Level Energy (keV)   Gamma energy (keV)\n",m+1);
						for(n=0;n<gd->nuclData[nucl].cascades[m].numLevels;n++)
							printf("%i      %f\n",n+1,gd->nuclData[nucl].cascades[m].energies[n]);
					}
				//if(strcmp(gd->nuclData[i].nuclName,"68SE")==0)
					//getc(stdin);
			}
}

void showLevelData(int nucl, gdata *gd)
{
	int i,j;
	double finalEnergy;

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
