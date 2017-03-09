#include "peak_find.h"

double getZero(int maxX, int minX, double maxY, double minY)
{
	double slope=(minY-maxY)/((double)(minX-maxX));
	double zeroPos=((-1.*maxY)/slope) + maxX;
	return zeroPos;
}

//function to find a peak, given an approximate value for its centroid
//takes the derivative, large spikes in the filter correspond to peaks
//the derivative is used to smooth out noise, and make the peak position correspond to a zero crossing in the filter output
//contraction - energy contraction (keV/ch)
//threshold - threshold (in sigma) for detection of a peak
//maxPeaks - maximum number of peaks to find (will be ordered by intensity)
peak_fit_par findPeak(const double * data, double contraction, double threshold, int maxPeaks)
{ 

	peak_fit_par par;

	//trapezoidal filter parameters
	int windowSize=1; //increase to soften the effect of noise and peaks
	int windowSpacing=1; //spacing between windows in the filter (mininum 1)
	int windowShift=(int)windowSize/2.;//number of channels to shift the windows by

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
	for(i=0;i<windowShift;i++)
		par.filterValue[i]=0.;
	for(i=windowShift;i<endCh;i++)
		if((i+windowSize+windowSpacing)<S32K)//check that we won't go out of bounds
			{
				win1val=0;
				win2val=0;
				for(j=0;j<windowSize;j++)
					{
						win1val+=data[i+j-windowShift];
						win2val+=data[i+j+windowSpacing-windowShift];
					}
				par.filterValue[i]=win2val-win1val;// trapezoidal filter value
			}
		else
			{
				par.filterValue[i]=0.;
			}
	
	
	//use smoother filter parameters for 2nd derivative 
	windowSize=3; //increase to soften the effect of noise and peaks
	windowSpacing=3; //spacing between windows in the filter (mininum 1)
	windowShift=(int)windowSize/2.;//number of channels to shift the windows by
	
	//take 2nd derivative
	for(i=0;i<windowShift;i++)
		par.filter2Value[i]=0.;
	for(i=windowShift;i<endCh;i++)
		if((i+windowSize+windowSpacing)<S32K)//check that we won't go out of bounds
			{
				win1val=0;
				win2val=0;
				for(j=0;j<windowSize;j++)
					{
						win1val+=par.filterValue[i+j-windowShift];
						win2val+=par.filterValue[i+j+windowSpacing-windowShift];
					}
				par.filter2Value[i]=win2val-win1val;// trapezoidal filter value
				if(par.filter2Value[i]>maxval)
					maxval=par.filter2Value[i];
				if(par.filter2Value[i]<minval)
					minval=par.filter2Value[i];
				maxch=i;
			}
		else
			{
				par.filter2Value[i]=0.;
			}
	
	double avg=0.;
	double stdev=0.;
	for(i=startCh;i<maxch;i++)
		avg+=par.filter2Value[i];
	if(maxch!=0)
		avg=avg/maxch;
	for(i=startCh;i<maxch;i++)
		stdev+=(par.filter2Value[i]-avg)*(par.filter2Value[i]-avg);
	if(maxch!=0)
		stdev=stdev/maxch;
	stdev=sqrt(stdev);
			
	//set thresholds for peak detection
	minThreshold=-1*(threshold+1)*stdev;
	maxThreshold=-1*threshold*stdev;
	
	printf("Filter output mean: %10.3f, stdev: %10.3f, min: %10.3f, max: %10.3f\n",avg,stdev,minval,maxval);
	//printf("Threshold 1: %10.3f, Threshold 2: %10.3f\n",maxThreshold,minThreshold);
		
	//sweep through the filter output and find peaks
	int fallingFlag=0;
	int risingFlag=0;
	int x1,x2;
	par.numPeaksFound=0;
	for(i=startCh;i<maxch;i++)
		{
			if(par.numPeaksFound<MAXPEAKSTOFIND)
				{
					if(par.filter2Value[i]<maxThreshold)
						if(i>0)
							{
								fallingFlag=1;
								x1=i-1;
							}
					if((fallingFlag==1)&&(par.filter2Value[i]<minThreshold))
						{
							fallingFlag=0;
							risingFlag=1;
						}
					if((risingFlag==1)&&(par.filter2Value[i]>maxThreshold))
						{
							risingFlag=0;
							x2=i;
							par.centroid[par.numPeaksFound]=(x1+x2)/2.;
							//printf("centroid: %f\n",par.centroid[par.numPeaksFound]);
							//getc(stdin);
							par.intensity[par.numPeaksFound]=data[(int)par.centroid[par.numPeaksFound]];//get intensity of peak
							par.numPeaksFound++;//register that a peak was found
						}
				}
		}
	
	//sort peaks by intensity in descending order
	double a;
	double ai;
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
	if(par.numPeaksFound>maxPeaks)
		par.numPeaksFound=maxPeaks;
	
	//contract peak positions to the proper energy 
	for(i=0;i<par.numPeaksFound;i++)
		par.centroid[i]=par.centroid[i]*contraction;
	
	return par;
}

void reportPeakPositions(peak_fit_par * par)
{
	int i;
	
	if(par->numPeaksFound>1)
		printf("%i peaks found at energies:",par->numPeaksFound);
	else if(par->numPeaksFound==1)
		printf("Peak found at energy:");
	else if(par->numPeaksFound==1)
		printf("No peaks found.\n");
		
	for(i=0;i<par->numPeaksFound;i++)
		{
			if(i==0)
				printf(" %5.1f", par->centroid[i]);
			else
				printf(", %5.1f", par->centroid[i]);
		}
	if(par->numPeaksFound>0)
		printf(" keV.\n");
}
