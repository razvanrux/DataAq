/*************************************************************************************
	PROIECT ACHIZITIA SI PRELUCRAREA DATELOR 
				2023-2024
			Ruxandari Razvan
				1307B
**************************************************************************************/
#include <utility.h>
#include <ansi_c.h>
#include <cvirte.h>		
#include <userint.h>
#include <stdio.h>
#include <stdlib.h>
#include <formatio.h>	//for formatting the file to an array -> FileToArray
#include <analysis.h>	//for the min, amx, dispersion etc.
#include "project.h"
#include <analysis.h>
#include <userint.h>
#include <string.h>
#include "toolbox.h"

/*************************************************************************************

	constante
**************************************************************************************/

#define SAMPLE_RATE			0
#define NPOINTS				1
#define NUM_POINTS_FOR_SKEWNESS 256 // numarul de esantioane pe care lucram cu skewness si kurtosis

// uhh da
#define CLAMP(value, min, max) ((value) < (min) ? (min) : ((value) > (max) ? (max) : (value)))


/*************************************************************************************
	Variabile globale
**************************************************************************************/

static int mainPanel;

	//wav32 -> used file
int waveInfo[2]; //waveInfo[0] = sampleRate
				 //waveInfo[1] = number of elements

	//frecventa de esantionare
double sampleRate = 0.0;
	//numarul de esantioane
int npoints = 0;
	//vector cu datelee semnalului audio -> signal in doc
double* waveData = 0;
double* anvelopa = 0;
double* deriv = 0;
	//vector cu datele semnalului filtrat in timp
double *filterTimeData = 0 ;
	//val. parametrului alpha
double* alpha;
	//intervalul de valori in care putem vedea
int startValue = 0;
int endValue = 6;
//int id=0;
double *filtered;
double windowSize;

int N=0;

WindowConst winConst;


		//valoarea maxima
	double maxVal = 0.0;
		//valoarea minima
	double minVal = 0.0;
		//indexul valorii maxime
	int maxIndex = 0;
		//indexul valorii minime
	int minIndex = 0;
		//valoarea mediei
	double mean = 0.0;
		//valoarea medianei
	double median = 0.0;
		//valoarea dispersiei
	double dispersion = 0.0;
		//numarul de 0
	int zeroCrossing = 0;
		int bitmapIDraw=0;
		int bitmapIDfiltered=0;
		int anvelopeCheck = 0;
		int derivCheck = 0;
		int interv=10;
double axis[100];
int hist[100];
int selectedFilter=1;
double *helps;

static int panel2;
/*************************************************************************************
	                      Functia ce arata graficul
**************************************************************************************/
void refreshWave (int panel, int minVal, int maxVal)
{
	waveData = (double*) calloc (npoints, sizeof(double));
	filtered= (double*) calloc (npoints, sizeof(double));
	anvelopa = (double*) calloc (npoints, sizeof(double));
	deriv = (double *) calloc(npoints, sizeof(double));
			//incarcare din fisierul .txt in memorie (vector)
			FileToArray("waveData.txt", waveData, VAL_DOUBLE, npoints, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			FileToArray("anvelopa.txt", anvelopa, VAL_DOUBLE, npoints, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			
			//incarcare/desenare pe graful RAW
			DeleteGraphPlot(mainPanel, MAIN_PANEL_IDC_GRAPH_RAW_DATA,-1, VAL_IMMEDIATE_DRAW);
			PlotY(mainPanel, MAIN_PANEL_IDC_GRAPH_RAW_DATA, waveData, npoints, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);	
}

/*************************************************************************************
	         Functia de filtrare (pe 16 sau 32, se mentioneaza windowsSize)
**************************************************************************************/
void DisplayFilteredWaveformInGraph(int panel, int start, int end, double *data, int numPoints) {
    // Calculate the indices corresponding to the specified start and end timestamps
    int startIndex = start * sampleRate; // Assuming sampleRate is in Hz
    int endIndex = end * sampleRate;     // Assuming sampleRate is in Hz

	 // Ensure startIndex and endIndex do not exceed the bounds of the waveData array
    startIndex = CLAMP(startIndex, 0, numPoints - 1);
    endIndex = CLAMP(endIndex, 0, numPoints - 1);
	
    // Clear the graph before adding new data
    DeleteGraphPlot(panel, MAIN_PANEL_IDC_GRAPH_FILTER_DATA, -1, VAL_IMMEDIATE_DRAW);

    // Plot the specified portion of the filtered waveform
    PlotY(panel, MAIN_PANEL_IDC_GRAPH_FILTER_DATA, &filtered[startIndex], endIndex - startIndex + 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
}

void DisplayDerivative(int panel, int start, int end, double *data, int numPoints) {
    // Calculate the indices corresponding to the specified start and end timestamps
    int startIndex = start * sampleRate; // Assuming sampleRate is in Hz
    int endIndex = end * sampleRate;     // Assuming sampleRate is in Hz

    // Ensure startIndex and endIndex do not exceed the bounds of the data array
    startIndex = CLAMP(startIndex, 1, numPoints - 1);
    endIndex = CLAMP(endIndex, 1, numPoints - 1);

    DifferenceEx(&data[startIndex], endIndex - startIndex + 1, 1.0, &data[startIndex], 1, &deriv[startIndex], 1, DIFF_FORWARD, deriv);
    // Plot the calculated derivative
	PlotY(panel, MAIN_PANEL_IDC_GRAPH_FILTER_DATA, &data[startIndex], endIndex - startIndex + 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_GREEN);
    //PlotY(panel, MAIN_PANEL_IDC_GRAPH_FILTER_DATA, data, endIndex - startIndex + 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_BLACK);
}


void DisplayWaveformInGraph(int panel, int start, int end, double *data, int numPoints) {
    // Calculate the indices corresponding to the specified start and end timestamps
    int startIndex = start * sampleRate; // Assuming sampleRate is in Hz
    int endIndex = end * sampleRate;     // Assuming sampleRate is in Hz

	 // Ensure startIndex and endIndex do not exceed the bounds of the waveData array
    startIndex = CLAMP(startIndex, 0, numPoints - 1);
    endIndex = CLAMP(endIndex, 0, numPoints - 1);
	
    // Clear the graph before adding new data
		 DeleteGraphPlot(panel, MAIN_PANEL_IDC_GRAPH_RAW_DATA, -1, VAL_IMMEDIATE_DRAW);
		 PlotY(panel, MAIN_PANEL_IDC_GRAPH_RAW_DATA, &data[startIndex], endIndex - startIndex + 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
    // Plot the specified portion of the filtered waveform
}

void DisplayAnvelope(int panel, int start, int end, double *data, int numPoints) {
    // Calculate the indices corresponding to the specified start and end timestamps
    int startIndex = start * sampleRate; // Assuming sampleRate is in Hz
    int endIndex = end * sampleRate;     // Assuming sampleRate is in Hz

	 // Ensure startIndex and endIndex do not exceed the bounds of the waveData array
    startIndex = CLAMP(startIndex, 0, numPoints - 1);
    endIndex = CLAMP(endIndex, 0, numPoints - 1);
	
    PlotY(panel, MAIN_PANEL_IDC_GRAPH_RAW_DATA, &data[startIndex], endIndex - startIndex + 1, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_BLUE); 
}

void DisplayHistogram(int panel, int start, int end, double *data, int numPoints) {
    // Calculate the indices corresponding to the specified start and end timestamps
    int startIndex = start * sampleRate; // Assuming sampleRate is in Hz
    int endIndex = end * sampleRate;     // Assuming sampleRate is in Hz

	 // Ensure startIndex and endIndex do not exceed the bounds of the waveData array
    startIndex = CLAMP(startIndex, 0, numPoints - 1);
    endIndex = CLAMP(endIndex, 0, numPoints - 1);
	Histogram(waveData+startIndex,endIndex-startIndex+1,minVal-1,maxVal+1, hist,axis,interv);
	DeleteGraphPlot (panel,MAIN_PANEL_GRAPH_HISTOGRAMA, -1, VAL_DELAYED_DRAW);
	PlotXY (panel,MAIN_PANEL_GRAPH_HISTOGRAMA, axis,  hist, interv, VAL_DOUBLE, VAL_SSIZE_T, VAL_VERTICAL_BAR, VAL_EMPTY_SQUARE, VAL_SOLID, 1, VAL_RED);
}

void AverageFilter(double *signal, double *filtered, int numPoints, int windowSize) {
    double sum = 0.0;

    // Calculate the initial sum for the first window
    for (int i = 0; i < windowSize; i++) {
        sum += signal[i];
    }

    // Apply the average filter for the first window
    for (int j = 0; j < windowSize; j++) {
        filtered[j] = sum / windowSize;
    }

    // Update the sum and apply the average filter for the subsequent windows
    for (int k = windowSize; k < numPoints; k++) {
        sum = sum - signal[k - windowSize] + signal[k];
        filtered[k] = sum / windowSize;
    }
}

void CustomFilter(double *signal, double *filtered, int numPoints) {
    filtered[0] = signal[0]; // Initialize the first element of the filtered signal

    for (int i = 1; i < numPoints; i++) {
        filtered[i] = (1 - *alpha) * filtered[i - 1] + *alpha * signal[i];
    }
}

void ApplyAveragingFilter(int panel, int option, double *filtered) {
	
    // Assuming waveData is your input audio signal and npoints is the number of data points
    if (option == 2) {
        CustomFilter(waveData, filtered, npoints);
    } else if (option == 1 || option == 0) {
        // Apply the averaging filter based on the specified window size
		if (option == 0)
		{
			AverageFilter(waveData, filtered, npoints, 32);
		}
        if (option == 1)
		{
			AverageFilter(waveData, filtered, npoints, 16);
		}
    }

    // Update the graph with the filtered signal
    DisplayFilteredWaveformInGraph(panel, startValue, endValue, filtered, npoints);
	if (derivCheck==1)
	{
		DisplayDerivative(panel, startValue, endValue, deriv, npoints);
	}
}

int ZeroCross(double *data, int points) {
    zeroCrossing = 0;

   for (int i = 1; i < points-1; i++) {
        // Check if the product of adjacent samples is negative (sign change) 
        if ((data[i - 1] >= 0 && data[i] <= 0) || (data[i - 1] <= 0 && data[i] >= 0) || data[i]==0) 
		{
            zeroCrossing++;
        }
    }
	
	return zeroCrossing;
}


void SkewnessAndKurtosis (double* data)
{
	 //skewness
    	double finSkew, skew1, skew2, maj; 
		skew1=0;
		skew2=0;
    	double tempSum = 0;
    	for(int i = 0; i <= 255; i++) {
        	tempSum += data[i];
    	}
    	maj = tempSum / 256;
    
    	for(int i = 0; i <= 255; i++) {
        	skew1 += (data[i] - maj) * (data[i] - maj) * (data[i] - maj);
        	skew2 += (data[i] - maj) * (data[i] - maj);
    	}
		
    	skew1 = skew1/256;
    	skew2 = skew2/256;
    
    	skew2 = skew2*skew2*skew2;
		skew2=sqrt(skew2);
    	finSkew = skew1 / skew2;
    
    	SetCtrlVal(mainPanel, MAIN_PANEL_COLORNUM1, finSkew); 

	//kurtosis
    double finKurt, kurt1, kurt2; 
    kurt1=0;
	kurt2=0;
    for(int i = 0; i <= 255; i++) {
        kurt1 += (data[i] - maj) * (data[i] - maj) * (data[i] - maj) * (data[i] - maj);
        kurt2 += (data[i] - maj) * (data[i] - maj);
    }
    kurt1 = kurt1/256;
    kurt2 = kurt2/256;
    
    kurt2 = kurt2*kurt2;
    finKurt = kurt1 / kurt2;
    
    SetCtrlVal(mainPanel, MAIN_PANEL_COLORNUM2, finKurt);
} 


/*************************************************************************************
	incarcarea fisierelor din wav pentru prelucrare ulterioara
**************************************************************************************/


double *loadWAV(const char* filePath, int *sampleRate, int *numSamples) {
	FILE *file = fopen(filePath, "rb");
	if(file == NULL)
		return NULL;
	// Citeste informatiile din antetul fisierului WAV
    fread(sampleRate, sizeof(int), 1, file);
    fread(numSamples, sizeof(int), 1, file);

    // Alocare memorie pentru semnal
    double *signal = (double*)malloc(sizeof(double) * (*numSamples));
    if (signal == NULL) {
        fclose(file);
        return NULL;
    }

    // Citire semnal
    fread(signal, sizeof(double), *numSamples, file);
    fclose(file);
    return signal;
}



/*************************************************************************************
								ZA main function
**************************************************************************************/

int main (int argc, char *argv[])
{
	int error = 0;
	
	/* initialize and load resources */
	nullChk (InitCVIRTE (0, argv, 0));
	errChk (mainPanel = LoadPanel (0, "project.uir", MAIN_PANEL));
	errChk (panel2 = LoadPanel (0, "project.uir", MAIN));
	
	/* display the panel and run the user interface */
	errChk (DisplayPanel (mainPanel));
	errChk (RunUserInterface ());

Error:
	/* clean up */
	if (mainPanel > 0)
		DiscardPanel (mainPanel);
	return 0;
}

/*************************************************************************************
								MAIN PANEL
**************************************************************************************/

int CVICALLBACK OnMainPanelCB (int panel, int event, void *callbackData,
							   int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_GOT_FOCUS:
			
			break;
		case EVENT_LOST_FOCUS:

			break;
		case EVENT_CLOSE:
			QuitUserInterface(0);
			break;
	}
	return 0;
}


/*************************************************************************************
								"Load Wave" le button
**************************************************************************************/

int CVICALLBACK OnLoadButtonCB (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{

	
	switch (event)
	{
		case EVENT_COMMIT:
			
			//incarcam informatiile din fisierele generate din wav.
			FileToArray("wafeInfo.txt", waveInfo, VAL_INTEGER, 2, 1, VAL_GROUPS_TOGETHER, VAL_GROUPS_AS_COLUMNS, VAL_ASCII);
			sampleRate = waveInfo[SAMPLE_RATE];
			npoints = waveInfo[NPOINTS];
			alpha = (double*)malloc(sizeof(double));
			helps = (double*)malloc(sizeof(double));
			*alpha = 0.01;
			
			refreshWave(MAIN_PANEL, maxVal, minVal);
			SetCtrlAttribute (mainPanel,MAIN_PANEL_IDC_GRAPH_RAW_DATA , ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_IDC_GRAPH_FILTER_DATA, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_PREV_BUTTON, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_IDC_TIME_START, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_IDC_TIME_STOP, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_NEXT_BUTTON, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_DIM_FEREASTRA, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_APPLY_BUTTON, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_MEAN_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_MIN_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_MAX_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_MIN_INDEX_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_MAX_INDEX_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_MEDIAN_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_DISPERSION_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_ZERO_CROSS_VAL, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_FILTER, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_COLORNUM1, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_COLORNUM2, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_COMMANDBUTTON, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_BINARYSWITCH, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_BINARYSWITCH2, ATTR_DIMMED, 0);
			SetCtrlAttribute (mainPanel, MAIN_PANEL_GRAPH_HISTOGRAMA, ATTR_DIMMED, 0);

			
			
			break;
	}
	return 0;
}

/*************************************************************************************
							DA main function
**************************************************************************************/

int CVICALLBACK OnPrevButtonCB (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	switch (event)
    {
        case EVENT_COMMIT:
		{
		GetCtrlVal(panel, MAIN_PANEL_IDC_TIME_STOP, &endValue);
        GetCtrlVal(panel, MAIN_PANEL_IDC_TIME_START, &startValue);
		GetCtrlVal(panel, MAIN_PANEL_BINARYSWITCH, &anvelopeCheck);
		GetCtrlVal(panel, MAIN_PANEL_BINARYSWITCH2, &derivCheck);
		GetCtrlVal(panel, MAIN_PANEL_DIM_FEREASTRA, &windowSize);
		GetCtrlVal(panel, MAIN_PANEL_FILTER, &selectedFilter);
            if (startValue > 0) {
                startValue--;
                endValue--;
                SetCtrlVal(panel, MAIN_PANEL_IDC_TIME_STOP, endValue);
                SetCtrlVal(panel, MAIN_PANEL_IDC_TIME_START, startValue);
				DisplayWaveformInGraph(panel, startValue, endValue, waveData, npoints);
				if (anvelopeCheck==1)
//				{
					DisplayAnvelope(panel, startValue, endValue, anvelopa, npoints);
				}
				ApplyAveragingFilter(panel, selectedFilter, filtered);
				//DisplayFilteredWaveformInGraph(panel, startValue, endValue, waveData, npoints);
				if (derivCheck==1)
				{
					DisplayDerivative(panel, startValue, endValue, deriv, npoints);
				}
            }
            // Calculate indices corresponding to the specified time range
            int startIndex = startValue * sampleRate;
            int endIndex = endValue * sampleRate;
            
            // Ensure the indices are within the bounds of the data
            startIndex = CLAMP(startIndex, 0, npoints - 1);
            endIndex = CLAMP(endIndex, 0, npoints - 1);
            
            // Calculate the number of points in the specified range
            int nrp = endIndex - startIndex + 1;
			//printf("%d, ", startIndex);
			//printf("%d\n", endIndex);
            
            // Perform analysis on the subset of the waveData array
            MaxMin1D(waveData + startIndex, nrp, &maxVal, &endIndex, &minVal, &startIndex);
            Mean(waveData + startIndex, nrp, &mean);
            Median(waveData + startIndex, nrp, &median);
            StdDev(waveData + startIndex, nrp, &mean, &dispersion);
			SkewnessAndKurtosis(waveData + startIndex);
			
			DisplayHistogram(panel, startValue, endValue, waveData, npoints);
            
            // Update the corresponding UI elements
            SetCtrlVal(panel, MAIN_PANEL_MEAN_VAL, mean);
            SetCtrlVal(panel, MAIN_PANEL_MIN_VAL, minVal);
            SetCtrlVal(panel, MAIN_PANEL_MAX_VAL, maxVal);
            SetCtrlVal(panel, MAIN_PANEL_MIN_INDEX_VAL, startIndex);
            SetCtrlVal(panel, MAIN_PANEL_MAX_INDEX_VAL, endIndex);
            SetCtrlVal(panel, MAIN_PANEL_MEDIAN_VAL, median);
            SetCtrlVal(panel, MAIN_PANEL_DISPERSION_VAL, dispersion);
            SetCtrlVal(panel, MAIN_PANEL_ZERO_CROSS_VAL,  ZeroCross(waveData + startIndex, nrp));
		
            
            break;
    }
    return 0;
}

int CVICALLBACK OnNextButtonCB (int panel, int control, int event,
								void *callbackData, int eventData1, int eventData2)
{
	switch (event)
    {
         case EVENT_COMMIT:
			  
            
            GetCtrlVal(panel, MAIN_PANEL_IDC_TIME_STOP, &endValue);
            GetCtrlVal(panel, MAIN_PANEL_IDC_TIME_START, &startValue);
			GetCtrlVal(panel, MAIN_PANEL_BINARYSWITCH, &anvelopeCheck);
			GetCtrlVal(panel, MAIN_PANEL_BINARYSWITCH2, &derivCheck);
			GetCtrlVal(panel, MAIN_PANEL_DIM_FEREASTRA, &windowSize);
			GetCtrlVal(panel, MAIN_PANEL_FILTER, &selectedFilter);
            if (endValue < 6) {
				
                startValue++;
                endValue++;
                SetCtrlVal(panel, MAIN_PANEL_IDC_TIME_STOP, endValue);
                SetCtrlVal(panel, MAIN_PANEL_IDC_TIME_START, startValue);
				DisplayWaveformInGraph(panel, startValue, endValue, waveData, npoints);
				if (anvelopeCheck==1)
				{
					DisplayAnvelope(panel, startValue, endValue, anvelopa, npoints);
				}
				ApplyAveragingFilter(panel, selectedFilter, filtered);
				//DisplayFilteredWaveformInGraph(panel, startValue, endValue, waveData, npoints);
				if (derivCheck==1)
				{
					DisplayDerivative(panel, startValue, endValue, deriv, npoints);
				}
            }
            // Calculate indices corresponding to the specified time range
            int startIndex = startValue * sampleRate;
            int endIndex = endValue * sampleRate;
			//printf("%d, ", startIndex);
			//printf("%d\n", endIndex);
            
            // Ensure the indices are within the bounds of the data
            startIndex = CLAMP(startIndex, 0, npoints - 1);
            endIndex = CLAMP(endIndex, 0, npoints - 1);
            
            // Calculate the number of points in the specified range
            int nrp = endIndex - startIndex + 1;
			
            int zeros=ZeroCross(waveData + startIndex, nrp);
			
            // Perform analysis on the subset of the waveData array
            MaxMin1D(waveData + startIndex, nrp, &maxVal, &endIndex, &minVal, &startIndex);
            Mean(waveData + startIndex, nrp, &mean);
            Median(waveData + startIndex, nrp, &median);
            StdDev(waveData + startIndex, nrp, &mean, &dispersion);
			SkewnessAndKurtosis(waveData + startIndex);
			DisplayHistogram(panel, startValue, endValue, waveData, npoints);
			
            
            // Update the corresponding UI elements
            SetCtrlVal(panel, MAIN_PANEL_MEAN_VAL, mean);
            SetCtrlVal(panel, MAIN_PANEL_MIN_VAL, minVal);
            SetCtrlVal(panel, MAIN_PANEL_MAX_VAL, maxVal);
            SetCtrlVal(panel, MAIN_PANEL_MIN_INDEX_VAL, startIndex);
            SetCtrlVal(panel, MAIN_PANEL_MAX_INDEX_VAL, endIndex);
            SetCtrlVal(panel, MAIN_PANEL_MEDIAN_VAL, median);
            SetCtrlVal(panel, MAIN_PANEL_DISPERSION_VAL, dispersion);
            SetCtrlVal(panel, MAIN_PANEL_ZERO_CROSS_VAL, zeros);
			
			

           break;
    }
    return 0;
}

int CVICALLBACK OnApplyButtonCB(int panel, int control, int event,
                                 void *callbackData, int eventData1, int eventData2)
{
    switch (event)
    {
        case EVENT_COMMIT:			
            // Calculate the alpha coefficient based on the window size
            *alpha = 1.0 / (double)(*helps + 1); // Adjust the calculation as needed
            break;
    }
    return 0;
}




int CVICALLBACK OnFilterClickCB(int panel, int control, int event,
                                 void *callbackData, int eventData1, int eventData2)
{
    switch (event)
    {
        case EVENT_COMMIT:
            GetCtrlVal(panel, control, &selectedFilter);
            switch (selectedFilter) {
                case 0: // Mediere 32
                    ApplyAveragingFilter(panel, 0, filtered);
                    break;
                case 1: // Mediere 16
                    ApplyAveragingFilter(panel, 1, filtered);
                    break;
                case 2: // Filtru
                    ApplyAveragingFilter(panel, 2, filtered);
                    break;
                default:
                    break;
            }
            break;
    }
    return 0;
}

int CVICALLBACK SaveGraphs(int panel, int control, int event,
                            void *callbackData, int eventData1, int eventData2)
{
    switch (event)
    {
        case EVENT_COMMIT:
			{
            // Calculate indices corresponding to the current second
            // Ensure the indices are within the bounds of the data

            // Capture and save the bitmap for the current interval
            GetCtrlDisplayBitmap(panel, MAIN_PANEL_IDC_GRAPH_RAW_DATA, 1, &bitmapIDraw);
			GetCtrlDisplayBitmap(panel, MAIN_PANEL_IDC_GRAPH_FILTER_DATA, 1, &bitmapIDfiltered);
			char fileName[256];
			char fileName1[256];
            sprintf(fileName, "RAW_DATA(%d-%d).jpg", startValue, startValue+1);
			sprintf(fileName1, "FILTERED_DATA(%d-%d).jpg", startValue, startValue+1);
            SaveBitmapToJPEGFileAnsi(bitmapIDraw, fileName, JPEG_PROGRESSIVE, 100);
			SaveBitmapToJPEGFileAnsi(bitmapIDfiltered, fileName1, JPEG_PROGRESSIVE, 100);
            break;
    		}
	}
    return 0;
}

int CVICALLBACK GenAnv (int panel, int control, int event,
						void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			break;
	}
	return 0;
}

int CVICALLBACK GenDeriv (int panel, int control, int event,
						  void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:

			break;
	}
	return 0;
}

int CVICALLBACK GetVal (int panel, int control, int event,
						void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
				GetCtrlVal(panel, MAIN_PANEL_DIM_FEREASTRA, helps);
			break;
	}
	return 0;
}



int CVICALLBACK ChangePanel (int panel, int control, int event,
							 void *callbackData, int eventData1, int eventData2)
{
	switch (event)
	{
		case EVENT_COMMIT:
			// Hide the current panel

            // Determine which panel to display next (mainpanel or panel2)
            if (panel == mainPanel)
            {
                DisplayPanel(panel2);
				HidePanel(mainPanel);
            }
            else
            {
                DisplayPanel(mainPanel);
				HidePanel(panel2);
            }
			break;
	}
	return 0;
}

int CVICALLBACK OnFrequencyButtonCB (int panel, int control, int event,
								  void *callbackData, int eventData1, int eventData2)
{
	GetCtrlVal(panel, MAIN_POINTS, &N);
	WindowConst winConst;
	
	double *convertedSpectrum; //vectorul ce contine spectrul semnalului convertit in volti
	double *autoSpectrum; //spectrul de putere
	double df=0.0; //pasul in domeniul frecventei
	double freqPeak=0.0; //valoarea maxima din spectrul de putere
	double powerPeak=0.0; //frecventa estimata pentru spectrum de putere
	
	char unit[32]="V";  //voltage semnal
	convertedSpectrum=(double*)calloc(npoints/6,sizeof(double));
	autoSpectrum=(double*)calloc(npoints/6,sizeof(double));
	
	switch (event)
	{
		case EVENT_COMMIT:
			
			ScaledWindowEx (waveData,N, RECTANGLE_, 0, &winConst);
			//se calculeaza partea pozitiva a spectrului scalat de putere pentru un semnal esantionat
			AutoPowerSpectrum(waveData,npoints/6, 1.0/sampleRate, autoSpectrum, &df);
			//calculeaza puterea si frecventa corespunzatoare varfului din spectrul semnalului
			PowerFrequencyEstimate(autoSpectrum,npoints/6,-1.0,winConst,df,7,&freqPeak,&powerPeak);
			
			SetCtrlVal(panel,MAIN_FREQ_PEAK,freqPeak);
			SetCtrlVal(panel,MAIN_POWER_PEAK,powerPeak);
			
			//se converteste spectrul de intrare în formate ce permit o reprezentare grafica mai convenabila
			SpectrumUnitConversion(autoSpectrum, N,0, SCALING_MODE_LINEAR, DISPLAY_UNIT_VRMS, df, winConst, convertedSpectrum, unit);
			
			//afisam spectrul pe grafic
			DeleteGraphPlot(panel,MAIN_GRAPH_SPECTRU,-1,VAL_IMMEDIATE_DRAW);
			PlotWaveform(panel,MAIN_GRAPH_SPECTRU, convertedSpectrum, npoints/12 ,VAL_DOUBLE, 1.0, 0.0, 0.0, df,VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID,  VAL_CONNECTED_POINTS, VAL_RED);
			
			break;
	}
	return 0;
}

int CVICALLBACK OnFilterButtonCB (int panel, int control, int event,
								  void *callbackData, int eventData1, int eventData2)
{
int secunda;
	double *raw;
	int winType;
	double final[npoints/6];
	double window[npoints/6];
	
	int fcut;
	int order;
	int fpass;
	int signalType;
	static WindowConst winConst;
	
	double *powerSpectrum;
	double *frequencyArray;
	
	char unit[32]="V";
	double df=0.0; //pasul in domeniul frecventei
		
	powerSpectrum=(double*)calloc(npoints/12,sizeof(double));
	frequencyArray=(double*)calloc(npoints/12,sizeof(double));
	
		switch (event)
	{
			
		case EVENT_COMMIT:
			
			//ferestruirea
			DeleteGraphPlot (panel,MAIN_GRAPH_WINDOW, -1, VAL_IMMEDIATE_DRAW);
			GetCtrlVal(panel, MAIN_POINTS, &N);
			
			GetCtrlVal(panel,MAIN_IDC_SECUNDA,&secunda);
			
			GetCtrlVal(panel,MAIN_IDC_CUT_FREQ,&fcut);
			GetCtrlVal(panel,MAIN_IDC_ORDER,&order);
			GetCtrlVal(panel,MAIN_IDC_FPASS,&fpass);
			
			raw=(double*)calloc(npoints/6,sizeof(double));
			for(int i=0;i<npoints/6;i++)
			{
				raw[i]=waveData[secunda*npoints/6+i];
			}	 
			
			GetCtrlVal(panel,MAIN_WINDOW_TYPE,&winType);
			
			
			//afisam semnalul pe secunde
			DeleteGraphPlot (panel, MAIN_GRAPH_RAW_DATA, -1, VAL_IMMEDIATE_DRAW);
			PlotY (panel, MAIN_GRAPH_RAW_DATA, raw, npoints/6, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			
			
			LinEv1D(raw,npoints/6,0,1,window);	   
			
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			switch (winType)
			{
				case 0: //CosWin
					GenCosWin (window, N, window, N);
					break;
					
				case 1: //Kaiser
					KsrWin (window, N, 4.86);
					break;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			DeleteGraphPlot (panel,MAIN_GRAPH_WINDOW, -1, VAL_IMMEDIATE_DRAW);
			PlotY (panel,MAIN_GRAPH_WINDOW, window, npoints/6, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			
			
			Mul1D(raw,window,npoints/6,final);
			DeleteGraphPlot (panel,MAIN_GRAPH_RAW_WINDOW, -1, VAL_IMMEDIATE_DRAW);
			PlotY (panel,MAIN_GRAPH_RAW_WINDOW, final, npoints/6, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			
			
			
			GetCtrlVal(panel,MAIN_FILTER_TYPE,&signalType);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			double *filteredSignal;
			filteredSignal = (double *) calloc(npoints, sizeof(double));
			
			switch(signalType)
			{
				case 0: //butterworth 5
				
					Bw_LPF(window, npoints/18, sampleRate, fcut, 5, filteredSignal);
					for (int i = 0; i < npoints/18; i++) {
   						 final[i] = filteredSignal[i];
					}
			
					break;
					
				case 1: // chebysev 5
 				
					Ch_HPF (final, npoints/24, npoints/24, fcut, 0.1, order, filteredSignal);
				   	break;

						   
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			DeleteGraphPlot (panel, MAIN_GRAPH_FILTER, -1, VAL_IMMEDIATE_DRAW);
			PlotY (panel, MAIN_GRAPH_FILTER, final, npoints/6, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);
			
				
			//afisam spectrul semnalului dupa ce am aplicat filtrul
			DeleteGraphPlot (panel,MAIN_GRAPH_SPECTRUM_FILTER, -1, VAL_IMMEDIATE_DRAW);
        	AutoPowerSpectrum (filteredSignal, npoints/6, 1.0/sampleRate, powerSpectrum, &df);
			SpectrumUnitConversion(powerSpectrum, npoints/12, 0, SCALING_MODE_LINEAR, DISPLAY_UNIT_VPK, df, winConst,frequencyArray, unit);
			PlotY (panel,MAIN_GRAPH_SPECTRUM_FILTER, frequencyArray, npoints/12, VAL_DOUBLE, VAL_THIN_LINE, VAL_EMPTY_SQUARE, VAL_SOLID, VAL_CONNECTED_POINTS, VAL_RED);		
	}
	return 0;
}
