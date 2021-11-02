
#include <fftw3.h>
#include <iostream>
#include <utils/microphone.hpp>
#include <math/math.hpp>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <thread>

#define OUTPUT_ENABLED 1

#define NUM_POINTS 128
#define SAMPLE_DURATION 0.01

#define NUMBER_OF_FRAMES 100

#define MINIMAL_SIGNAL_THRESHOLD 1 //this is the threshold for the noise reduction
#define MINIMAL_DIFFERENTIAL_OFFSET 0.05

using namespace std;

#define REAL 0
#define IMAG 1

float sampleFreq = NUM_POINTS / SAMPLE_DURATION;
float magnitudes[6][NUMBER_OF_FRAMES];
float magDifferentials[6][NUMBER_OF_FRAMES];
float scaledFreqBand[6];

float unprocessedMagnitudes[NUM_POINTS/2];

float rollingMagnitudeDifferential[NUM_POINTS/2];
float rollingMagnitudeDifferentialAverage[NUM_POINTS/2];
float rollingMagnitude[NUM_POINTS/2];


float frequencydelta  = sampleFreq / NUM_POINTS;

//the buffer index for the ringbuffers
int magBufferIdx = 0;

//just and counter that counted how often there were not enough audiosamples when accessed. It is no longer required.
int bufferMisses = 0;


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

void storeMagnitudeValues(int frequencyBin, float mag)
{
    
    //take the rolling average of the magnitude differential
    rollingMagnitude[frequencyBin] = rollingMagnitude[frequencyBin] * (0.2f) + mag * (0.8f);
    magnitudes[frequencyBin][magBufferIdx] = rollingMagnitude[frequencyBin];

    float magDiff = magnitudes[frequencyBin][magBufferIdx] - magnitudes[frequencyBin][(magBufferIdx + NUMBER_OF_FRAMES-1)%NUMBER_OF_FRAMES];
    rollingMagnitudeDifferential[frequencyBin] = rollingMagnitudeDifferential[frequencyBin] * (0.4f) + magDiff * (0.6f);
    magDifferentials[frequencyBin][magBufferIdx] = rollingMagnitudeDifferential[frequencyBin];

    rollingMagnitudeDifferentialAverage[frequencyBin] = rollingMagnitudeDifferentialAverage[frequencyBin] * 0.999f + rollingMagnitude[frequencyBin] * 0.001f; 
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
    float onsetThreshold = rollingMagnitudeDifferentialAverage[frequencyBin] + MINIMAL_DIFFERENTIAL_OFFSET;
    if( rollingMagnitudeDifferential[frequencyBin]>onsetThreshold)
        return 1;
    return 0;

}

void processResults(fftwf_complex* result)
{

    

    magBufferIdx = (magBufferIdx +1 ) % NUMBER_OF_FRAMES;

}


int startMicRoutine()
{
    float* buffer = (float*) malloc(sizeof(float) * 512);
    initMicrophone(sampleFreq, NUM_POINTS);

    float signal[NUM_POINTS];
    fftwf_complex result[NUM_POINTS];


    while(1){
        auto start = chrono::high_resolution_clock::now();

        Pa_Sleep(1);
        if (!readCurrentBuffer(buffer))
        {

            fftwf_plan plan = fftwf_plan_dft_r2c_1d(NUM_POINTS,
                                            signal,
                                            result,
                                            0);

            for(int i = 0; i < NUM_POINTS; i ++)
            {
                signal[i] = buffer[i];
            }
            //acquire_from_somewhere(signal);
            fftwf_execute(plan);
            processResults(result);
            fftwf_destroy_plan(plan);
            bufferMisses = 0;

            auto stop = chrono::high_resolution_clock::now();
            

            auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
            if(duration < std::chrono::milliseconds(16)){
                //std::this_thread::sleep_for(std::chrono::milliseconds(16) - (duration));
            }
  
            // To get the value of duration use the count()
            // member function on the duration object
            if(OUTPUT_ENABLED)
                cout << "millisec per cycle"  << duration.count() << endl;


        }
        else 
            bufferMisses ++;

    }
    

    releaseMicrophone();
}


int main()
    {

    startMicRoutine();   
    for(int i = 0; i < 6; i ++)
    {
        rollingMagnitudeDifferential[i] = 0;
        rollingMagnitudeDifferentialAverage[i] = 0;
        rollingMagnitude[i] = 0;
    }

   


    return 0;
}