
#include <fftw3.h>
#include <iostream>
#include <utils/microphone.hpp>
#include <math/math.hpp>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <thread>


#include <map>
#include <vector>
#include <cmath>

#include "gnuplot-iostream.h"

#define OUTPUT_ENABLED 1
#define NUM_POINTS 1028
#define SAMPLE_DURATION 0.02
#define NUMBER_OF_FRAMES 200
#define MINIMAL_SIGNAL_THRESHOLD 0.1 //this is the threshold for the noise reduction
#define NUMBER_OF_COMBFILTERS 100
#define COMFILTER_INCREMENT 0.5
#define COMBFILTER_START_BPM 70
#define GENERAL_GAIN 1

using namespace std;

#define REAL 0
#define IMAG 1

float GAIN_PER_FREQUENCY[] = {2, 2, 1, 1, 1, 1};
float SMOOTHING_COEFFICIENT_PER_FREQ[] = {0.84, 0.84, 0.6, 0.6, 0.6, 0.6};
float SMOOTHING_COEFFICIENT_PER_FREQ_DIFF[] = {0.5, 0.5, 0.7, 0.7, 0.8, 0.8};
//float SMOOTHING_COEFFICIENT_PER_FREQ[] = {0, 0, 0, 0, 0, 0};
// float SMOOTHING_COEFFICIENT_PER_FREQ_DIFF[] = {0, 0, 0, 0, 0, 0};

float sampleFreq = NUM_POINTS / SAMPLE_DURATION;
float scaledFreqBand[6];
float unprocessedMagnitudes[NUM_POINTS/2]; 
float magnitues[6][NUMBER_OF_FRAMES];
float magnitudeDifferentials[6][NUMBER_OF_FRAMES];

float alignedSumDifferentialBuffer[NUMBER_OF_FRAMES]; //this is the plotting buffer fo the ringbufferdifferentialovertime
float ringBufferDifferentialOverTime[NUMBER_OF_FRAMES]; // this stores the summed up differential over time

float rollingAverage[6];
float alignedRollingDifferentialSumAverage[NUMBER_OF_FRAMES]; // this buffer is only for plotting. Old Averages have to be kept over time. This is the buffer that gets aligned
float rollingDifferentialSumAverage[NUMBER_OF_FRAMES]; //this buffer is only for plotting. Old Averages have to be kept over time

float combFilterbank[6][NUMBER_OF_COMBFILTERS][NUMBER_OF_FRAMES]; //this buffers stores the summed up energy of each comb filter
float combFilterEnergySum[NUMBER_OF_COMBFILTERS];
int combFilterPeakBin[6];
float alignedCombFilterONEBANDONLY[NUMBER_OF_FRAMES];//ONLY FOR PLOTTING
int voteHistory[NUMBER_OF_FRAMES];
int voteCounter[NUMBER_OF_COMBFILTERS];

float loopFrequencyDelta = 1; //this value gets updated each loop. It is the frequency Delta
float updatesPerSecond = 24; //number of updates of gnuplot per second

float frequencydelta  = sampleFreq / NUM_POINTS;

fftwf_complex magDifferentialResults[NUMBER_OF_FRAMES];

//the buffer index for the ringbuffers
int magBufferIdx = 0;
//just and counter that counted how often there were not enough audiosamples when accessed. It is no longer required.
int bufferMisses = 0;


int peakBin = 0;
int votedPeakBin = 0;
 Gnuplot gp;
/**
 * @brief This method can be used to generate debug signals. Otherwise it has no purpose
 * 
 * @param signal the signal buffer
 */
void acquire_from_somewhere(float* signal) {

    int i;
    for (i = 0; i < NUM_POINTS; i++) {
        float theta = ((float)i / (float)NUM_POINTS) * SAMPLE_DURATION;
        signal[i] = 1.0 * cos(510.0 * 2*M_PI * theta) +  2.0 * cos(300.0 * 2*M_PI * theta );
    }
}


void plotDifferentialBuffer()
{
   
    std::vector<std::pair<double, double> > xy_pts_A;
    for(int x=0; x< NUMBER_OF_FRAMES; x++) {
        double y = alignedSumDifferentialBuffer[x] -0.5;
        xy_pts_A.push_back(std::make_pair(x, y));
    }
    std::vector<std::pair<double, double> > xy_pts_B;
    for(int x=0; x< NUMBER_OF_FRAMES; x++) {
        double y = alignedRollingDifferentialSumAverage[x] - 0.5;
        xy_pts_B.push_back(std::make_pair(x, y));
    }

    std::vector<std::pair<float, float> > xy_pts_C;
    float barwidth = NUMBER_OF_FRAMES/8;
    for(int i = 0; i < 6; i ++)
    {
        float y = magnitues[i][(magBufferIdx - 1 + NUMBER_OF_FRAMES)%NUMBER_OF_FRAMES];
        xy_pts_C.push_back(std::make_pair((i * barwidth + barwidth/2.0), y  * 0.04 + 0.5));
    }

    std::vector<std::pair<float, float> > xy_pts_D;
    for(int x=0; x< NUMBER_OF_FRAMES; x++)
    {
        float y = alignedCombFilterONEBANDONLY[x];
        xy_pts_D.push_back(std::make_pair(x,y));
    }

    std::vector<std::pair<float, float> > xy_pts_E;
    barwidth = NUMBER_OF_FRAMES/NUMBER_OF_COMBFILTERS;
    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i ++)
    {
        float y = combFilterEnergySum[i] / 300.0;
        xy_pts_E.push_back(std::make_pair((i * barwidth + barwidth/2.0), y ));
    }



    gp << "set xrange [0:" + to_string(NUMBER_OF_FRAMES) + "]\nset yrange [-1:7]\n";
   
    gp << "plot" 
    << gp.file1d(xy_pts_A) << "with lines title 'Sum differential',"  
    << gp.file1d(xy_pts_B) << "with lines title 'avg'," 
    //<< gp.file1d(xy_pts_D) << "with lines title 'comb'," 
    << gp.file1d(xy_pts_C) << "with boxes title 'Amplitudes',"
    << gp.file1d(xy_pts_E) << "with boxes title 'comb',"
        << std::endl;

}

void suppressNoiseAmp(float* signal)
{
    //this might not be ideal. Since the (percieved) sounds at high frequency are spread out over a larger frequency band
    //this method of noise cancelation reduces high frequencies more. A possible solution would be to scale the signal threshold
    //down for higher frequencies.
    for (int i = 0; i < NUM_POINTS / 2; i ++)
    {
        if(signal[i] < MINIMAL_SIGNAL_THRESHOLD)
        {
            signal[i] = 0;
        }
    }
}


void alignMagnitudeRingBuffers()
{
    for(int i = 0; i < NUMBER_OF_FRAMES; i ++)
    {
        alignedSumDifferentialBuffer[i] = ringBufferDifferentialOverTime[(i + magBufferIdx) % NUMBER_OF_FRAMES];
        alignedRollingDifferentialSumAverage[i] = rollingDifferentialSumAverage[(i + magBufferIdx) % NUMBER_OF_FRAMES];
        alignedCombFilterONEBANDONLY[i] = combFilterbank[3][0][(i + magBufferIdx) % NUMBER_OF_FRAMES];
         
    }
}

void clearCombFilterEnergySum()
{
    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i++)
    {
        combFilterEnergySum[i] = 0;
    }
}

int findPeakBin()
{   
    int bin = 0;
    float peakMax = 0;
    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i++)
    {
        if(peakMax < combFilterEnergySum[i])
        {
            peakMax = combFilterEnergySum[i];
            bin = i;
        }
    }   
    return bin;
}

int voteLogic()
{
    //remove old vote
    voteCounter[voteHistory[magBufferIdx]] --;
    //add new vote
    voteCounter[peakBin] ++;
    voteHistory[magBufferIdx] = peakBin;

    //find overall Winner
    int maxCount = 0;
    int maxBin = 0;
    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i ++)
    {
        if(voteCounter[i] > maxCount)
        {
            maxCount = voteCounter[i];
            maxBin = i;
        }
    }

    votedPeakBin = maxBin;
}

float getCombFilterRecursive(float currentOffset, float generalOffset, int band)
{
    if(currentOffset > NUMBER_OF_FRAMES)
        return 0;
    
    int magIndexWithOffsetBottom = (magBufferIdx - (int)(currentOffset) + NUMBER_OF_FRAMES)%NUMBER_OF_FRAMES;
    int magIndexWithOffsetTop = (magBufferIdx - (int)(currentOffset) + 1 + NUMBER_OF_FRAMES)%NUMBER_OF_FRAMES;
    float residue = currentOffset - (int)currentOffset;

    float interpolMagDiffValue = magnitudeDifferentials[band][magIndexWithOffsetBottom] * (1 - residue) + magnitudeDifferentials[band][magIndexWithOffsetTop] * (residue);

    float alpha = pow(0.5, (float)currentOffset/(float)generalOffset);
    return (1-alpha)* interpolMagDiffValue + alpha*getCombFilterRecursive(currentOffset + generalOffset, generalOffset, band);
}

void applyCombFilterToFrequencyBand(int bandId)
{
   

    //apply comb filters
    //the comb filter with delay T: y_t = alpha * y_(t-T) + ( 1 - alpha ) * x_t 

    float delay_bpm = COMBFILTER_START_BPM; //in bpm
    float delay_steps_bpm = COMFILTER_INCREMENT;
    float alpha = 0;
    float numberOfTimeSteps = 0;    

    float combMax = 0;
    int combBin = 0;

    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i++)
    {
        delay_bpm = COMBFILTER_START_BPM + ( i * delay_steps_bpm );
        float delay_in_seconds = 60.0/delay_bpm;

        numberOfTimeSteps = delay_in_seconds/SAMPLE_DURATION;
        combFilterbank[bandId][i][magBufferIdx] = getCombFilterRecursive(0, numberOfTimeSteps, bandId);
        float localMax = 0;
        float combSum = 0;
        for(int a = 0; a < NUMBER_OF_FRAMES; a++)
        {
            if(combFilterbank[bandId][i][a] > combMax)
            {
                combMax = combFilterbank[bandId][i][a];
                combBin = i;
            }
            if(localMax < combFilterbank[bandId][i][a])
            {
                localMax = combFilterbank[bandId][i][a];
            }
            combSum += combFilterbank[bandId][i][a];
        }
        combFilterEnergySum[i] += localMax/combSum;
        combFilterPeakBin[bandId] = combBin; 
            


    }
    cout << "Peak bin for band" <<  bandId << " : " <<combFilterPeakBin[bandId] << "with " << (combFilterPeakBin[bandId] * delay_steps_bpm + COMBFILTER_START_BPM) <<endl;


}


void findDominantFrequency()
{
    for(int i = 0; i < NUMBER_OF_FRAMES; i ++)
    {
        alignedSumDifferentialBuffer[i] = ringBufferDifferentialOverTime[(i + magBufferIdx) % NUMBER_OF_FRAMES];
        alignedRollingDifferentialSumAverage[i] = rollingDifferentialSumAverage[(i + magBufferIdx) % NUMBER_OF_FRAMES];

       //alignedSumDifferentialBuffer[i] *= math::applyHanningFunction(i, NUMBER_OF_FRAMES);

    }
    fftwf_plan plan = fftwf_plan_dft_r2c_1d( NUMBER_OF_FRAMES,
                                            alignedSumDifferentialBuffer,
                                            magDifferentialResults,
                                            0);

    //acquire_from_somewhere(signal);
    fftwf_execute(plan);


    //now find the dominant frequency!!
    float diffFreqDelta = (1/SAMPLE_DURATION) / (NUMBER_OF_FRAMES);

    int dominantBin = 0;
    float dominantAmp = 0;
    std::system("clear");
    for(int i = 60; i <  90; i ++)
    {
        float curAmp = sqrt(magDifferentialResults[i][REAL] * magDifferentialResults[i][REAL] + magDifferentialResults[i][IMAG] * magDifferentialResults[i][IMAG]);
        if(curAmp > dominantAmp){
            dominantAmp = curAmp;
            dominantBin = i;
        }

        string out = to_string((float)i * diffFreqDelta  * 60);
        for(int a = 0; a < curAmp/3 ; a ++)
        {
            out = out + "#";
        }
        cout << out << endl;
    }
    cout << "dominant frequency " << ((float)dominantBin * diffFreqDelta * 60.0) << endl;
    cout << "time window length:" << NUMBER_OF_FRAMES * SAMPLE_DURATION << endl;

    fftwf_destroy_plan(plan);

}



void processResults(fftwf_complex* result)
{
    //std::system("clear");
    //first compute the magnitudes
    for(int i = 0; i < NUM_POINTS/2; i++)
    {
        unprocessedMagnitudes[i] = (sqrt(result[i][REAL] * result[i][REAL] + result[i][IMAG] * result[i][IMAG]));
        
    }
    //Now apply noise reduction to the unprocessed magnitudes
    suppressNoiseAmp(unprocessedMagnitudes);
    //reduce the frequency bands to 6 channels
    math::getScaledFrequencyBands(unprocessedMagnitudes, scaledFreqBand, frequencydelta);
    int previousMagBuffer = (magBufferIdx - 1 + NUMBER_OF_FRAMES)%NUMBER_OF_FRAMES;
    //now apply some rolling average smoothing and store them in the buffer
    for(int i = 0; i < 6; i++)
    {
        
        magnitues[i][magBufferIdx] = ( ( GAIN_PER_FREQUENCY[i] * scaledFreqBand[i] * (1 - SMOOTHING_COEFFICIENT_PER_FREQ[i])) 
        + (SMOOTHING_COEFFICIENT_PER_FREQ[i] * magnitues[i][previousMagBuffer]));
    }

    //now calculate the discrete derivative of the magnitudes(and apply a small rolling average to the signal)
    for(int i = 0; i < 6; i++)
    {
        float currentDifferential =  max(0.0f, magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - (rollingAverage[i] * 2) );
        //float currentDifferential =  magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - rollingAverage[i];
        rollingAverage[i] = rollingAverage[i] * 0.999 + currentDifferential * 0.001;
        // magnitudeDifferentials[i][magBufferIdx] = magnitudeDifferentials[i][previousMagBuffer] * 0.005 + currentDifferential * 0.995 ;
        magnitudeDifferentials[i][magBufferIdx] =  currentDifferential* (1-SMOOTHING_COEFFICIENT_PER_FREQ_DIFF[i])
             + magnitudeDifferentials[i][previousMagBuffer] * SMOOTHING_COEFFICIENT_PER_FREQ_DIFF[i] ;
    }

    //just for fun add up all current differentials
    float differentialSum = 0;
    float averageSum = 0;

    for(int a = 0; a < 6; a ++)
    {
        averageSum += rollingAverage[a];
        differentialSum+= magnitudeDifferentials[a][magBufferIdx];
        
    }
    differentialSum /= 6;
    string out = "Diff:";
    for(int i = 0; i < (differentialSum  * 18) + 5; i ++)
    {
        out = out + "#";
    }
    
    //cout << out << endl;

    rollingDifferentialSumAverage[magBufferIdx] = averageSum/6;
    ringBufferDifferentialOverTime[magBufferIdx] = differentialSum;
    
    alignMagnitudeRingBuffers();
    // findDominantFrequency();
    clearCombFilterEnergySum();
    system("clear");
    for(int i = 0; i < 5; i++)
        applyCombFilterToFrequencyBand(i);


    int peakbin = findPeakBin();
    peakBin = peakbin;
    voteLogic();
    cout << "peak bin: " << peakbin << " at freq: " << (peakbin * COMFILTER_INCREMENT + COMBFILTER_START_BPM) << endl;
    cout << "vote peak bin: " << votedPeakBin << " at freq: " << (votedPeakBin * COMFILTER_INCREMENT + COMBFILTER_START_BPM) << endl;
    magBufferIdx = (magBufferIdx +1 ) % NUMBER_OF_FRAMES;

}

void preprocessInputSignal(float* signal)
{
    for(int i = 0; i < NUM_POINTS ; i ++)
    {
        signal[i] = signal[i] * math::applyHanningFunction(i, NUM_POINTS) * GENERAL_GAIN;
    }

}

int startMicRoutine()
{
    float* buffer = (float*) malloc(sizeof(float) * 512);
    initMicrophone(sampleFreq, NUM_POINTS);

    float signal[NUM_POINTS];
    fftwf_complex result[NUM_POINTS];

    auto plotUpdateTimer = chrono::high_resolution_clock::now();

    while(1){

        Pa_Sleep(1);
        if (!readCurrentBuffer(signal))
        {
            auto start = chrono::high_resolution_clock::now();
            preprocessInputSignal(signal);

            fftwf_plan plan = fftwf_plan_dft_r2c_1d(NUM_POINTS,
                                            signal,
                                            result,
                                            0);

            //acquire_from_somewhere(signal);
            fftwf_execute(plan);
            processResults(result);
            fftwf_destroy_plan(plan);
            bufferMisses = 0;

            auto stop = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
            loopFrequencyDelta = loopFrequencyDelta * 0.995 +  (duration.count())* 0.005;
  
            // To get the value of duration use the count()
            // member function on the duration object
            //if(OUTPUT_ENABLED)
                // cout << "millisec per cycle"  << loopFrequencyDelta << endl;
                // cout << "millisec per cycle"  << duration.count() << endl;


        }
        else {
            auto stopUpdate = std::chrono::high_resolution_clock::now();
            auto lastUpdate = std::chrono::duration_cast<chrono::milliseconds>(stopUpdate - plotUpdateTimer);
            
            bufferMisses ++;
            if(lastUpdate.count() > 1000/updatesPerSecond)
            {
                plotUpdateTimer = chrono::high_resolution_clock::now();
                plotDifferentialBuffer();
            }
        }

    }
    releaseMicrophone();
}

int main()
{
    startMicRoutine();   

    return 0;
}