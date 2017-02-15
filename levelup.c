//definitions
#include "levelup.h"
//functions and logic
#include "parse_ENSDF.c"

int main(int argc, char *argv[])
{
	if(!getenv("ENSDF")) 
		{
			printf("\nThe environment variable ENSDF is not set.\nWithout it, this program cannot find the ENSDF data files.  Set the envirnoment variable to point to the directory containing ENSDF data files using the following command (replace the path with the appropriate one):\n\n");
			printf ("export ENSDF=/path/to/ensdf/\n\n");
			exit(-1);
		}

  if(argc!=2)
    {
      printf("\nlevelup gamma_list\n----------------------\n\n");
      printf("Tells you what species are present in the\nlist of gamma ray energies you provide.\n\n");
      exit(-1);
    }
  
  //allocate structures
  gdata *gd=(gdata*)malloc(sizeof(gdata));
  initialize_database(gd);
	printf("Structure allocated\n");
	
	char fileName[256],str[8];
	int i;
	
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

  return 0; //great success
}
