#include "peak_find.h"

//function to find a peak, given an approximate value for its centroid
//uses a trapezoidal filter, large spikes in the filter correspond to peaks
//contraction - energy contraction (keV/ch)
//threshold - threshold (in sigma) for detection of a peak
//maxPeaks - maximum number of peaks to find (will be ordered by intensity)
peak_fit_par findPeak(const double * data, double contraction, double threshold, int maxPeaks)
{ 

	peak_fit_par par;

	//trapezoidal filter parameters
	int windowSize=1; //increase to soften the effect of noise and peaks
	int windowSpacing=1; //spacing between windows in the filter

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
	
	double avg=0.;
	double stdev=0.;
	for(i=startCh;i<maxch;i++)
		avg+=par.filterValue[i];
	if(maxch!=0)
		avg=avg/maxch;
	for(i=startCh;i<maxch;i++)
		stdev+=(par.filterValue[i]-avg)*(par.filterValue[i]-avg);
	if(maxch!=0)
		stdev=stdev/maxch;
	stdev=sqrt(stdev);
			
	//set thresholds for peak detection
	maxThreshold=threshold*stdev;
	minThreshold=-1*threshold*stdev;
	
	printf("Filter output mean: %10.3f, stdev: %10.3f, min: %10.3f, max: %10.3f\n",avg,stdev,minval,maxval);
	
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
					par.intensity[par.numPeaksFound]=data[par.centroid[par.numPeaksFound]];//get intensity of peak
					par.numPeaksFound++;//register that a peak was found
				}
		}
	
	//sort peaks by intensity in descending order
	double a;
	int ai;
	for(i=0;i<par.numPeaksFound;i++)
		for(j=i+1;j<par.numPeaksFound;j++)
			if(par.intensity[i]<par.intensity[j])
				{
					//swap values
					a=par.intensity[i];
					par.intensity[i]=par.intensity[j];
					par.intensity[j]=a;
					ai=par.centroid[i];
					par.centroid[i]=par.centroid[j];
					par.centroid[j]=ai;
				}
	
	//for(i=0;i<par.numPeaksFound;i++)
	//	printf("Intesity[%i]: %f\n",i,par.intensity[i]);
	
	//only retain the number of peaks requested
	par.numPeaksFound=maxPeaks;
	
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
	if(par->numPeaksFound>0)
		printf(" keV.\n");
	else
		printf("No peaks found.\n");
}
