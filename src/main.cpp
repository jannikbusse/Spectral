
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
#define NUM_POINTS 128
#define SAMPLE_DURATION 0.05
#define NUMBER_OF_FRAMES 500
#define MINIMAL_SIGNAL_THRESHOLD 0.05 //this is the threshold for the noise reduction
#define NUMBER_OF_COMBFILTERS 100

using namespace std;

#define REAL 0
#define IMAG 1

float GAIN_PER_FREQUENCY[] = {4, 3, 2, 1, 1, 1};
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

float combFilterbank[6][NUMBER_OF_COMBFILTERS]; //this buffers stores the summed up energy of each comb filter

float loopFrequencyDelta = 1; //this value gets updated each loop. It is the frequency Delta
float updatesPerSecond = 24; //number of updates of gnuplot per second

float frequencydelta  = sampleFreq / NUM_POINTS;

fftwf_complex magDifferentialResults[NUMBER_OF_FRAMES];

//the buffer index for the ringbuffers
int magBufferIdx = 0;
//just and counter that counted how often there were not enough audiosamples when accessed. It is no longer required.
int bufferMisses = 0;

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
        double y = alignedSumDifferentialBuffer[x];
        xy_pts_A.push_back(std::make_pair(x, y));
    }
    std::vector<std::pair<double, double> > xy_pts_B;
    for(int x=0; x< NUMBER_OF_FRAMES; x++) {
        double y = alignedRollingDifferentialSumAverage[x];
        xy_pts_B.push_back(std::make_pair(x, y));
    }
    gp << "set xrange [0:" + to_string(NUMBER_OF_FRAMES) + "]\nset yrange [-1:3]\n";
   
    gp << "plot" << gp.file1d(xy_pts_A) << "with lines title 'Sum differential',"  << gp.file1d(xy_pts_B) << "with lines title 'avg',"
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
         
    }
}

void applyCombFilterToFrequencyBand(int bandId)
{
   

    //apply comb filters
    //the comb filter with delay T: y_t = alpha * y_(t-T) + ( 1 - alpha ) * x_t 
    //for alpha = 0.5^(t/T)
    //note that the frequencyDelta is given by SAMPLE TIME
    //loopfrequencydelta is given in sec
    
    float delay_T = 1; //in seconds
    float alpha = 0;
    float numberOfTimeSteps = 0;
    for(int i = 0; i < NUMBER_OF_COMBFILTERS; i++)
    {
        delay_T = 1 * i;
        numberOfTimeSteps = delay_T/SAMPLE_DURATION;
        for(int t = 0; t < NUMBER_OF_FRAMES; t++)
        {
            alpha = pow(0.5, (t/delay_T));
            
            combFilterbank[bandId][magBufferIdx] = magnitudeDifferentials[bandId][magBufferIdx];
            //add the comb filter signal
            combFilterbank[bandId][magBufferIdx] += (1-alpha)*combFilterbank[bandId][magBufferIdx];
            //att the first iteration of the filter
            
        }
        cout << numberOfTimeSteps << endl;


    }
    cout << NUMBER_OF_FRAMES * SAMPLE_DURATION <<endl;


}


void findDominantFrequency()
{
    for(int i = 0; i < NUMBER_OF_FRAMES; i ++)
    {
        alignedSumDifferentialBuffer[i] = ringBufferDifferentialOverTime[(i + magBufferIdx) % NUMBER_OF_FRAMES];
        alignedRollingDifferentialSumAverage[i] = rollingDifferentialSumAverage[(i + magBufferIdx) % NUMBER_OF_FRAMES];

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
    for(int i = 4; i <  NUMBER_OF_FRAMES/10; i ++)
    {
        float curAmp = sqrt(magDifferentialResults[i][REAL] * magDifferentialResults[i][REAL] + magDifferentialResults[i][IMAG] * magDifferentialResults[i][IMAG]);
        if(curAmp > dominantAmp){
            dominantAmp = curAmp;
            dominantBin = i;
        }

        string out = to_string((float)i * diffFreqDelta  * 60);
        for(int a = 0; a < curAmp ; a ++)
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
        unprocessedMagnitudes[i] = ((sqrt(result[i][REAL] * result[i][REAL] + result[i][IMAG] * result[i][IMAG])));
        
    }
    //Now apply noise reduction to the unprocessed magnitudes
    suppressNoiseAmp(unprocessedMagnitudes);
    //reduce the frequency bands to 6 channels
    math::getScaledFrequencyBands(unprocessedMagnitudes, scaledFreqBand);
    int previousMagBuffer = (magBufferIdx - 1 + NUMBER_OF_FRAMES)%NUMBER_OF_FRAMES;
    //now apply some rolling average smoothing and store them in the buffer
    for(int i = 0; i < 6; i++)
    {
        
        magnitues[i][magBufferIdx] = ( ( GAIN_PER_FREQUENCY[i] * scaledFreqBand[i] * 0.1) + (0.9 * magnitues[i][previousMagBuffer]));
    }

    //now calculate the discrete derivative of the magnitudes(and apply a small rolling average to the signal)
    for(int i = 0; i < 6; i++)
    {
        //float currentDifferential =  max(0.0f, magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - rollingAverage[i] - 1);
        float currentDifferential =  magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - rollingAverage[i];
        rollingAverage[i] = rollingAverage[i] * 0.997 + currentDifferential * 0.003;

        // magnitudeDifferentials[i][magBufferIdx] = magnitudeDifferentials[i][previousMagBuffer] * 0.005 + currentDifferential * 0.995 ;
        magnitudeDifferentials[i][magBufferIdx] = magnitudeDifferentials[i][previousMagBuffer] * 0.87 + currentDifferential * 0.13 ;
    }

    //just for fun add up all current differentials
    float differentialSum = 0;
    float averageSum = 0;

    for(int a = 0; a < 6; a ++)
    {
        averageSum += rollingAverage[i];
        differentialSum+= magnitudeDifferentials[a][magBufferIdx];
        
    }
    string out = "Diff:";
    for(int i = 0; i < (differentialSum  * 18) + 5; i ++)
    {
        out = out + "#";
    }
    cout << out << endl;

    rollingDifferentialSumAverage[magBufferIdx] = averageSum;
    ringBufferDifferentialOverTime[magBufferIdx] = differentialSum;
    
    alignMagnitudeRingBuffers();
    findDominantFrequency();

    //applyCombFilterToFrequencyBand(0);

    magBufferIdx = (magBufferIdx +1 ) % NUMBER_OF_FRAMES;


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
    for(int i = 0; i < 6; i ++)
    {
        
    }

   


    return 0;
}