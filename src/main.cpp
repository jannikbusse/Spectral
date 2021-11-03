
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
#define SAMPLE_DURATION 0.01

#define NUMBER_OF_FRAMES 400

#define MINIMAL_SIGNAL_THRESHOLD 0.05 //this is the threshold for the noise reduction

using namespace std;

#define REAL 0
#define IMAG 1

float sampleFreq = NUM_POINTS / SAMPLE_DURATION;
float scaledFreqBand[6];
float unprocessedMagnitudes[NUM_POINTS/2];
float magnitues[6][NUMBER_OF_FRAMES];
float magnitudeDifferentials[6][NUMBER_OF_FRAMES];

float allignedBuffer[NUMBER_OF_FRAMES];
float ringBufferDifferentialOverTime[NUMBER_OF_FRAMES];

float rollingAverage[6];
float alignedRollingDifferentialSumAverage[NUMBER_OF_FRAMES];
float rollingDifferentialSumAverage[NUMBER_OF_FRAMES];


float loopFrequency = 4000;

fftwf_complex magDifferentialResults[NUMBER_OF_FRAMES];

float updatesPerSecond = 10;





float frequencydelta  = sampleFreq / NUM_POINTS;

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
        double y = allignedBuffer[x];
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

int detectSignalOnset(int frequencyBin)
{
   

}

void findDominantFrequency()
{
    for(int i = 0; i < NUMBER_OF_FRAMES; i ++)
    {
        allignedBuffer[i] = ringBufferDifferentialOverTime[(i + magBufferIdx) % NUMBER_OF_FRAMES];
        alignedRollingDifferentialSumAverage[i] = rollingDifferentialSumAverage[(i + magBufferIdx) % NUMBER_OF_FRAMES];
         
    }
    fftwf_plan plan = fftwf_plan_dft_r2c_1d( NUMBER_OF_FRAMES,
                                            allignedBuffer,
                                            magDifferentialResults,
                                            0);

    //acquire_from_somewhere(signal);
    fftwf_execute(plan);


    //now find the dominant frequency!!
    float diffSampleFreq = 1000000 / loopFrequency;
    float diffFreqDelta = loopFrequency / (NUMBER_OF_FRAMES);

    int dominantBin = 0;
    float dominantAmp = 0;
    //std::system("clear");
    for(int i = 0; i <  NUMBER_OF_FRAMES/130; i ++)
    {
        float curAmp = sqrt(magDifferentialResults[i][REAL] * magDifferentialResults[i][REAL] + magDifferentialResults[i][IMAG] * magDifferentialResults[i][IMAG]);
        if(curAmp > dominantAmp){
            dominantAmp = curAmp;
            dominantBin = i;
        }

        string out = to_string(i * diffFreqDelta  * 60);
        for(int a = 0; a < curAmp ; a ++)
        {
            out = out + "#";
        }
        //cout << out << endl;
    }
    //cout << "dominant frequency " << ((float)dominantBin * diffFreqDelta * 60.0) << endl;

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
    math::getScaledFrequencyBands(unprocessedMagnitudes, scaledFreqBand);
    int previousMagBuffer = (magBufferIdx - 1 + NUMBER_OF_FRAMES)%NUMBER_OF_FRAMES;
    //now apply some rolling average smoothing and store them in the buffer
    for(int i = 0; i < 6; i++)
    {
        
        magnitues[i][magBufferIdx] = ( (scaledFreqBand[i] * 0.2) + (0.8 * magnitues[i][previousMagBuffer]));
    }

    //now calculate the discrete derivative of the magnitudes(and apply a small rolling average to the signal)
    for(int i = 0; i < 6; i++)
    {
        float currentDifferential =  max(0.0f, magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - rollingAverage[i] - 1);
        //float currentDifferential =  magnitues[i][magBufferIdx] - magnitues[i][previousMagBuffer] - rollingAverage[i];
        rollingAverage[i] = rollingAverage[i] * 0.999 + currentDifferential * 0.001;

        // magnitudeDifferentials[i][magBufferIdx] = magnitudeDifferentials[i][previousMagBuffer] * 0.005 + currentDifferential * 0.995 ;
        magnitudeDifferentials[i][magBufferIdx] = magnitudeDifferentials[i][previousMagBuffer] * 0.8 + currentDifferential * 0.2 ;
    }

    //just for fun add up all current differentials
    float differentialSum = 0;
    for(int a = 0; a < 6; a ++)
    {
        
        differentialSum+= magnitudeDifferentials[a][magBufferIdx];
        
    }
    string out = "Diff:";
    for(int i = 0; i < (differentialSum  * 18) + 5; i ++)
    {
        out = out + "#";
    }
    cout << out << endl;

    rollingDifferentialSumAverage[magBufferIdx] = differentialSum * 0.01 + rollingDifferentialSumAverage[previousMagBuffer] *  0.99;
    ringBufferDifferentialOverTime[magBufferIdx] = differentialSum;
    

    findDominantFrequency();

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
            

            if(duration < std::chrono::microseconds(16)){
                //std::this_thread::sleep_for(std::chrono::milliseconds(16) - (duration));
            }

            loopFrequency = loopFrequency * 0.99 +  (duration.count())* 0.01;
  
            // To get the value of duration use the count()
            // member function on the duration object
            //if(OUTPUT_ENABLED)
                // cout << "millisec per cycle"  << loopFrequency << endl;
                // cout << "millisec per cycle"  << duration.count() << endl;


        }
        else {
            auto stopUpdate = std::chrono::high_resolution_clock::now();
            auto lastUpdate = std::chrono::duration_cast<chrono::milliseconds>(stopUpdate - plotUpdateTimer);
            
            bufferMisses ++;
        if(lastUpdate.count() > 1000/updatesPerSecond)
            plotDifferentialBuffer();
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