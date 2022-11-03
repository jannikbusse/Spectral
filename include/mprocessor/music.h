#pragma once


#include <fftw3.h>
#include <iostream>
#include "../math/math.hpp"
#include "../utils/microphone.hpp"
#include "../utils/serial.hpp"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <thread>
#include<unistd.h>
#include <map>
#include <vector>
#include <cmath>
#include "gnuplot-iostream.h"


#define GNUPLOT_FRAMES 900

#define OUTPUT_ENABLED 1
#define NUM_STORED_SIGNAL_FRAMES 5
#define NUM_AUDIO_BUFFER_POINTS 512
#define NUM_AUDIO_SAMPLE_DURATION 0.02
#define NUM_POINTS (NUM_AUDIO_BUFFER_POINTS * NUM_STORED_SIGNAL_FRAMES)
#define SAMPLE_DURATION (NUM_AUDIO_SAMPLE_DURATION * NUM_STORED_SIGNAL_FRAMES)
#define MINIMAL_SIGNAL_THRESHOLD 0.4 //this is the threshold for the noise reduction
#define GENERAL_GAIN 4.6

using namespace std;

#define REAL 0
#define IMAG 1


float scaledFreqBand[6]; //this are the 6 frequency bands that are used after the first processing steps

float unprocessedMagnitudes[NUM_POINTS/2]; 

//fft signal and result buffers
float signal[NUM_STORED_SIGNAL_FRAMES][NUM_AUDIO_BUFFER_POINTS];
float alignedSignalBuffer[NUM_POINTS];
float avrgFreqMagn[NUM_POINTS/2];
int currentSignalBuffer = 0;
fftwf_complex result[NUM_POINTS];
//------------------------------


float sampleFreq = NUM_POINTS / SAMPLE_DURATION;
float frequencydelta  = sampleFreq / NUM_POINTS;

float longBassAvrg = 0;
float shortBassAvrg = 0;
float bassEnergy = 0;
float bassDiff = 0;
float highEnergy = 0;
float shortHighAvg = 0;
float longHighAvg = 0;
float newsum = 0;
float totalsum = 0;
float avgSmoothedSum = 0;
float unsmoothed = 0;
float currentVolume = 0;

sendStruct sendBuffer;


const int NUM_PLOTS = 7;
std::vector<std::pair<double, double> > xy_plotpoints[NUM_PLOTS];
float gnuPoints[NUM_PLOTS][GNUPLOT_FRAMES];
int skipFrames = 3;
int curFrame = 0;

static const auto& cout_alias = system("pkill -f gnu");
Gnuplot gp, gp2;

void plotGnu()
{   
    for(int i = 0; i < NUM_PLOTS; i++)
        xy_plotpoints[i].clear();
    //GNULOGIC
    //align buffer 
    for(int i = 0; i < NUM_PLOTS; i++)
        memmove(gnuPoints[i],gnuPoints[i]+1,sizeof(float)*(GNUPLOT_FRAMES-1));
    gnuPoints[0][GNUPLOT_FRAMES -1]= totalsum;
    gnuPoints[1][GNUPLOT_FRAMES -1]= scaledFreqBand[0];
    gnuPoints[2][GNUPLOT_FRAMES -1]= scaledFreqBand[1];
    gnuPoints[3][GNUPLOT_FRAMES -1]= scaledFreqBand[2];
    gnuPoints[4][GNUPLOT_FRAMES -1]= scaledFreqBand[3];
    gnuPoints[5][GNUPLOT_FRAMES -1]= scaledFreqBand[4];
    gnuPoints[6][GNUPLOT_FRAMES -1]= scaledFreqBand[5];
    cout <<currentVolume << endl;
    for(int x=0; x< GNUPLOT_FRAMES; x++) {
        for(int i=0; i < NUM_PLOTS; i++){
            double y = gnuPoints[i][x];
            xy_plotpoints[i].push_back(std::make_pair(x, y));
        }
    }
    gp << "plot"  
   // << gp.file1d(xy_plotpoints[0]) << "with lines title 'Sum differential',"  
    << gp.file1d(xy_plotpoints[1]) << "with lines title '0',"  
    << gp.file1d(xy_plotpoints[2]) << "with lines title '1',"  
    << gp.file1d(xy_plotpoints[3]) << "with lines title '2',"  
    << gp.file1d(xy_plotpoints[4]) << "with lines title '3',"  
    << gp.file1d(xy_plotpoints[5]) << "with lines title '4',"  
    << gp.file1d(xy_plotpoints[6]) << "with lines title '5',"  
    << std::endl;
}

int init()
{
    return initMicrophone(sampleFreq, NUM_AUDIO_BUFFER_POINTS);
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

float calculateVolume(float* signal)
{
    float sum = 0;
    for(int i = 0; i < NUM_POINTS;i++ )
    {
        sum += signal[i] * signal[i];
    }
    return sqrt(sum);
}

float RSME_smooth(float newVal, float oldVal, float newWeight)
{
    return sqrt((newVal * newVal * newWeight + oldVal *  oldVal *(1-newWeight)));
}
float ROLAVG_smooth(float newVal, float oldVal, float newWeight)
{
    return newVal * newWeight +  oldVal *(1-newWeight);
}


void processResults(fftwf_complex* result)
{
    printf("\033c");
    newsum = 0;
    for(int i = 1; i < NUM_POINTS/2; i++)
    {
        unprocessedMagnitudes[i] = (sqrt(result[i][REAL] * result[i][REAL] + result[i][IMAG] * result[i][IMAG]))/NUM_POINTS * 1000;
        
    }
    suppressNoiseAmp(unprocessedMagnitudes);

    math::getScaledFrequencyBands(unprocessedMagnitudes, scaledFreqBand, frequencydelta);

    float con = 0.99f;
    for(int i = 1; i < NUM_POINTS/2; i++)
    {

        newsum += max(0.0f, unprocessedMagnitudes[i] - avrgFreqMagn[i]);
        avrgFreqMagn[i] = avrgFreqMagn[i] * con + (1-con)*unprocessedMagnitudes[i];
        
    }
    avgSmoothedSum = ROLAVG_smooth(newsum, totalsum, 0.3f);

    totalsum = RSME_smooth(newsum, totalsum, 0.8f);
    //totalsum = max(pow(totalsum,0.99f), max(0.0f, totalsum));

    sendBuffer.vol = currentVolume;
    sendBuffer.energychange = totalsum;
    cout << "VOL " << currentVolume;
    cout << "SUM " << totalsum;
    //writeToArduino(sendBuffer);


}

void preprocessInputSignal(float* signal)
{
    for(int i = 0; i < NUM_POINTS ; i ++)
    {
        signal[i] = signal[i] * math::applyHanningFunction(i, NUM_POINTS) * GENERAL_GAIN;
    }
}


int releaseMusic()
{
    //delete all the arrays...
    return releaseMicrophone();
}

int update()
{
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
            fftwf_destroy_plan(plan);

  
            // To get the value of duration use the count()
            // member function on the duration object
            //if(OUTPUT_ENABLED)
                // cout << "millisec per cycle"  << loopFrequencyDelta << endl;
                // cout << "millisec per cycle"  << duration.count() << endl;
            currentSignalBuffer = (currentSignalBuffer + 1 ) %NUM_STORED_SIGNAL_FRAMES;

        }else{
            if(curFrame <= 0){
                plotGnu();
                curFrame = skipFrames;
            }
            curFrame --;
        }
    return 0;
}

