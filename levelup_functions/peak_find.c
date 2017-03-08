#include "peak_find.h"

//function to find a peak, given an approximate value for its centroid
//uses a trapezoidal filter, large spikes in the filter correspond to peaks
peak_fit_par findPeak(const double * data, double contraction)
{ 

	peak_fit_par par;

	//trapezoidal filter parameters
	int windowSize=5; //increase to soften the effect of noise and peaks
	int windowSpacing=2; //spacing between windows in the filter

	double win1val,win2val;//holds values for trapezoidal filter window(s)
	double maxval,minval;//holds maximum and minimum values from the filter
	double maxThreshold,minThreshold;//maximum and minimum threshold values for peak detection
	int maxch;//holds index of the end of the filter output
	int i,j;


	//set bounds for trapezoidal filter
	int startCh=0;
	int endCh=S32K;

	maxval=0;
	minval=BIG_NUMBER;

	//construct the trapezoidal filter output
	for(i=startCh;i<endCh;i++)
		if((i+windowSize+windowSpacing)<S32K)//check that we won't go out of bounds
			{
				win1val=0;
				win2val=0;
				for(j=0;j<windowSize;j++)
					{
						win1val+=data[i+j];
						win2val+=data[i+j+windowSpacing];
					}
				par.filterValue[i]=win2val-win1val;// trapezoidal filter value
				if(par.filterValue[i]>maxval)
					maxval=par.filterValue[i];
				if(par.filterValue[i]<minval)
					minval=par.filterValue[i];
				maxch=i;
			}
		else
			{
				par.filterValue[i]=0.;
			}
			
	//set thresholds for peak detection
	maxThreshold=maxval*0.25;
	minThreshold=minval*0.25;
	
	//sweep through the filter output and find peaks
	int risingFlag=0;
	int fallingFlag=0;
	par.numPeaksFound=0;
	for(i=startCh;i<maxch;i++)
		{
			if(par.filterValue[i]>maxThreshold)
				risingFlag=1;
			if((risingFlag==1)&&(par.filterValue[i]<=0.))
				if(par.numPeaksFound<MAXPEAKSTOFIND)
					{
						par.centroid[par.numPeaksFound]=i;
						risingFlag=0;//reset the flag
						fallingFlag=1;
					}
			if((fallingFlag==1)&&(par.filterValue[i]<minThreshold))
				{
					fallingFlag=0;//reset the flag
					par.numPeaksFound++;//register that a peak was found
				}
		}
	
	//contract peak positions to the proper energy 
	for(i=0;i<par.numPeaksFound;i++)
		par.centroid[i]=par.centroid[i]*contraction;
	
	
	return par;
}

void reportPeakPositions(peak_fit_par * par)
{
	int i;
	for(i=0;i<par->numPeaksFound;i++)
		{
			if(i==0)
				printf("Peak(s) found at energy: %i", par->centroid[i]);
			else
				printf(", %i", par->centroid[i]);
		}
	printf(" keV.\n");
}
