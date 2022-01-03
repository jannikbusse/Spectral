#pragma once


#include <fftw3.h>
#include <iostream>
#include "../math/math.hpp"
#include "../utils/microphone.hpp"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <thread>


#include <map>
#include <vector>
#include <cmath>

#include "gnuplot-iostream.h"

#define OUTPUT_ENABLED 1
#define NUM_STORED_SIGNAL_FRAMES 3
#define NUM_AUDIO_BUFFER_POINTS 512
#define NUM_AUDIO_SAMPLE_DURATION 0.01
#define NUM_POINTS (NUM_AUDIO_BUFFER_POINTS * NUM_STORED_SIGNAL_FRAMES)
#define SAMPLE_DURATION (NUM_AUDIO_SAMPLE_DURATION * NUM_STORED_SIGNAL_FRAMES)
#define NUMBER_OF_FRAMES 800
#define MINIMAL_SIGNAL_THRESHOLD 0.1 //this is the threshold for the noise reduction
#define NUMBER_OF_COMBFILTERS 150
#define COMFILTER_INCREMENT 0.3
#define COMBFILTER_START_BPM 60
#define GENERAL_GAIN 0.6

using namespace std;

#define REAL 0
#define IMAG 1

float GAIN_PER_FREQUENCY[] = {1.5, 1.5, 1, 1, 0.7, 0.7};
float SMOOTHING_COEFFICIENT_PER_FREQ[] = {0.74, 0.74, 0.74, 0.74, 0.74, 0.74};
float SMOOTHING_COEFFICIENT_PER_FREQ_DIFF[] = {0.65, 0.65, 0.65, 0.65, 0.7, 0.7};


float scaledFreqBand[6]; //this are the 6 frequency bands that are used after the first processing steps
float unprocessedMagnitudes[NUM_POINTS/2]; 
float magnitues[6][NUMBER_OF_FRAMES]; //this are the magnitudes over time
float magnitudeDifferentials[6][NUMBER_OF_FRAMES]; //this are the differentials over time

float alignedSumDifferentialBuffer[NUMBER_OF_FRAMES]; //this is the plotting buffer fo the ringbufferdifferentialovertime
float ringBufferDifferentialOverTime[NUMBER_OF_FRAMES]; // this stores the summed up differential over time only for plotting

float rollingAverage[6];
float alignedRollingDifferentialSumAverage[NUMBER_OF_FRAMES]; // this buffer is only for plotting. Old Averages have to be kept over time. This is the buffer that gets aligned
float rollingDifferentialSumAverage[NUMBER_OF_FRAMES]; //this buffer is only for plotting. Old Averages have to be kept over time

float combFilterbank[6][NUMBER_OF_COMBFILTERS][NUMBER_OF_FRAMES]; //this buffers stores the summed up energy of each comb filter for each frame
float combFilterEnergySum[NUMBER_OF_COMBFILTERS];
int combFilterPeakBin[6];
float alignedCombFilterONEBANDONLY[NUMBER_OF_FRAMES];//ONLY FOR PLOTTING

int voteHistory[NUMBER_OF_FRAMES];
int voteCounter[NUMBER_OF_COMBFILTERS];


bool onsetsDetected[3]; //bass high and one for later


//fft signal and result buffers
float signal[NUM_STORED_SIGNAL_FRAMES][NUM_AUDIO_BUFFER_POINTS];
float alignedSignalBuffer[NUM_POINTS];
int currentSignalBuffer = 0;
fftwf_complex result[NUM_POINTS];
//------------------------------


float loopFrequencyDelta = 1; //this value gets updated each loop. It is the frequency Delta
float updatesPerSecond = 24; //number of updates of gnuplot per second
float sampleFreq = NUM_POINTS / SAMPLE_DURATION;
float frequencydelta  = sampleFreq / NUM_POINTS;

fftwf_complex magDifferentialResults[NUMBER_OF_FRAMES];

//the buffer index for the ringbuffers
int magBufferIdx = 0;
//just and counter that counted how often there were not enough audiosamples when accessed. It is no longer required.

float currentVolume = 0;

int peakBin = 0;
int votedPeakBin = 0;

void (*bassOnsetCallback)(void);


static const auto& cout_alias = system("pkill -f gnu");


 Gnuplot gp, gp2;
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
    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i ++)
    {
        float y = combFilterEnergySum[i] / 200;
        xy_pts_E.push_back(std::make_pair((i * NUMBER_OF_FRAMES/NUMBER_OF_COMBFILTERS + NUMBER_OF_FRAMES/NUMBER_OF_COMBFILTERS/2.0), y ));
    }

    std::vector<std::pair<float, float> > volume;
    barwidth = NUMBER_OF_COMBFILTERS/8;
    float y = currentVolume * 4;
    volume.push_back(std::make_pair((3 * barwidth + barwidth/2.0), y  * 0.04 + 0.5));
    volume.push_back(std::make_pair((2 * barwidth + barwidth/2.0), y  * 0.04 + 0.5));



    gp << "set xrange [0:" + to_string(NUMBER_OF_FRAMES) + "]\nset yrange [-1:7]\n";
    gp2 << "set xrange [0:" + to_string(NUMBER_OF_COMBFILTERS) + "]\nset yrange [-1:7]\n";
   
    gp << "plot" 
    << gp.file1d(xy_pts_A) << "with lines title 'Sum differential',"  
    << gp.file1d(xy_pts_B) << "with lines title 'avg'," 
    //<< gp.file1d(xy_pts_D) << "with lines title 'comb'," 
    //<< gp.file1d(xy_pts_C) << "with boxes title 'Amplitudes',"
    //<< gp.file1d(xy_pts_E) << "with boxes title 'comb',"
        << std::endl;
    gp2 << "plot" 
    //<< gp.file1d(xy_pts_D) << "with lines title 'comb'," 
    << gp2.file1d(volume) << "with boxes title 'Amplitudes',"
    //<< gp.file1d(xy_pts_E) << "with boxes title 'comb',"
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
        alignedCombFilterONEBANDONLY[i] = combFilterbank[4][0][(i + magBufferIdx) % NUMBER_OF_FRAMES];
         
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
    return maxBin;
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

    float combMaxEnergy = 0;
    int combBin = 0;

    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i++)
    {
        delay_bpm = COMBFILTER_START_BPM + ( i * delay_steps_bpm );
        float delay_in_seconds = 60.0/delay_bpm;

        numberOfTimeSteps = delay_in_seconds/NUM_AUDIO_SAMPLE_DURATION;
        combFilterbank[bandId][i][magBufferIdx] = getCombFilterRecursive(0, numberOfTimeSteps, bandId);
        float energysum = 0;
        
        for(int a = 0; a < NUMBER_OF_FRAMES; a++)
        {
            energysum += combFilterbank[bandId][i][a]*combFilterbank[bandId][i][a];
        }
        combFilterEnergySum[i] += energysum;
        combFilterPeakBin[bandId] = combBin; 
        if(energysum > combMaxEnergy)
        {
            combMaxEnergy = energysum;
            combBin = i;
        }
            
    }
    // cout << "Peak bin for band" <<  bandId << " : " <<combFilterPeakBin[bandId] << "with " << (combFilterPeakBin[bandId] * delay_steps_bpm + COMBFILTER_START_BPM) 
    //  << " and energy peak at " << combMaxEnergy<<endl;


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
    float diffFreqDelta = (1/NUM_AUDIO_SAMPLE_DURATION) / (NUMBER_OF_FRAMES);

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
        //cout << out << endl;
    }
    //cout << "dominant frequency " << ((float)dominantBin * diffFreqDelta * 60.0) << endl;
    //cout << "time window length:" << NUMBER_OF_FRAMES * NUM_AUDIO_SAMPLE_DURATION << endl;

    fftwf_destroy_plan(plan);

}



void processResults(fftwf_complex* result)
{
    //std::system("clear");
    //first compute the magnitudes
    for(int i = 0; i < NUM_POINTS/2; i++)
    {
        unprocessedMagnitudes[i] = (sqrt(result[i][REAL] * result[i][REAL] + result[i][IMAG] * result[i][IMAG]))/NUM_POINTS * 1000;
        
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
        float currentDifferential =  max(0.0f, magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - (rollingAverage[i] * 4) -2 );
        //float currentDifferential =  magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - rollingAverage[i];
        rollingAverage[i] = rollingAverage[i] * 0.999 + currentDifferential * 0.001;
        // magnitudeDifferentials[i][magBufferIdx] = magnitudeDifferentials[i][previousMagBuffer] * 0.005 + currentDifferential * 0.995 ;
        magnitudeDifferentials[i][magBufferIdx] =  currentDifferential* (1-SMOOTHING_COEFFICIENT_PER_FREQ_DIFF[i])
             + magnitudeDifferentials[i][previousMagBuffer] * SMOOTHING_COEFFICIENT_PER_FREQ_DIFF[i] ;
    }

    //just for fun add up all current differentials
    float differentialSum = 0;
    float averageSum = 0;

    for(int a = 0; a < 2; a ++)
    {
        averageSum += rollingAverage[a];
        differentialSum+= magnitudeDifferentials[a][magBufferIdx];
        
    }
    differentialSum /= 2;
    string out = "Diff:";
    for(int i = 0; i < (differentialSum  * 18) + 5; i ++)
    {
        out = out + "#";
    }
    
    //cout << out << endl;

    rollingDifferentialSumAverage[magBufferIdx] = averageSum/2;
    ringBufferDifferentialOverTime[magBufferIdx] = differentialSum;
    
    alignMagnitudeRingBuffers();
    // findDominantFrequency();
    clearCombFilterEnergySum();
    // system("clear");
    for(int i = 0; i < 5; i++)
        applyCombFilterToFrequencyBand(i);


    int peakbin = findPeakBin();
    peakBin = peakbin;
    voteLogic();
    // cout << "frequencydelta "  << frequencydelta << endl;
    // cout << "time duration " << SAMPLE_DURATION << endl;

    // cout << "peak bin: " << peakbin << " at freq: " << (peakbin * COMFILTER_INCREMENT + COMBFILTER_START_BPM) << endl;
    // cout << "vote peak bin: " << votedPeakBin << " at freq: " << (votedPeakBin * COMFILTER_INCREMENT + COMBFILTER_START_BPM) << endl;

}

void preprocessInputSignal(float* signal)
{
    for(int i = 0; i < NUM_POINTS ; i ++)
    {
        signal[i] = signal[i] * math::applyHanningFunction(i, NUM_POINTS) * GENERAL_GAIN;
    }

}


float calculateVolume(float* signal)
{
    float sum = 0;
    for(int i = 0; i < NUM_POINTS;i++ )
    {
        sum += signal[i] * signal[i];
    }
    return sqrt(sum);
}

void handleCallbacks()
{
    float bassSum = magnitudeDifferentials[1][magBufferIdx] + magnitudeDifferentials[1][magBufferIdx];
    float bassAvg = rollingAverage[1] + rollingAverage[1] ;
    if(bassSum > 6*bassAvg)
    {
        if(!onsetsDetected[0]){
            bassOnsetCallback();
            onsetsDetected[0] = true;
        }
    }else
    {
        onsetsDetected[0] = false;
    }


}

int init(void bassCallback(void))
{
    bassOnsetCallback = bassCallback;
    return initMicrophone(sampleFreq, NUM_AUDIO_BUFFER_POINTS);
}

int releaseMusic()
{
    //delete all the arrays...
    return releaseMicrophone();
}

int update()
{
        Pa_Sleep(1);
        if (!readCurrentBuffer(signal[currentSignalBuffer]))
        {
            //align the buffers
            for(int i = 0; i < NUM_STORED_SIGNAL_FRAMES; i ++)
            {
                memcpy(alignedSignalBuffer + (i*NUM_AUDIO_BUFFER_POINTS), signal[(currentSignalBuffer + i)%NUM_STORED_SIGNAL_FRAMES], 
                                                    sizeof(signal[(currentSignalBuffer + i)%NUM_STORED_SIGNAL_FRAMES]));
            }

            currentVolume = calculateVolume(alignedSignalBuffer);
            preprocessInputSignal(alignedSignalBuffer);

            fftwf_plan plan = fftwf_plan_dft_r2c_1d(NUM_POINTS,
                                            alignedSignalBuffer,
                                            result,
                                            0);

            //acquire_from_somewhere(signal);
            fftwf_execute(plan);
            processResults(result);
            handleCallbacks();
            fftwf_destroy_plan(plan);

  
            // To get the value of duration use the count()
            // member function on the duration object
            //if(OUTPUT_ENABLED)
                // cout << "millisec per cycle"  << loopFrequencyDelta << endl;
                // cout << "millisec per cycle"  << duration.count() << endl;
            magBufferIdx = (magBufferIdx +1 ) % NUMBER_OF_FRAMES;
            currentSignalBuffer = (currentSignalBuffer + 1 ) %NUM_STORED_SIGNAL_FRAMES;

        }
        else {
            
                plotDifferentialBuffer();
        }

    return 0;
}

