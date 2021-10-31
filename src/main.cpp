
#include <fftw3.h>
#include <iostream>
#include <utils/microphone.hpp>
#include <math/math.hpp>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define NUM_POINTS 128
#define SAMPLE_DURATION 0.01

using namespace std;

#define REAL 0
#define IMAG 1

float sampleFreq = NUM_POINTS / SAMPLE_DURATION;
float *magnitudes = (float*) malloc(sizeof(float) * NUM_POINTS/2);

float scaledFreqBand[6][2];

void acquire_from_somewhere(fftw_complex* signal) {
    /* Generate two sine waves of different frequencies and
     * amplitudes.
     */

    int i;
    for (i = 0; i < NUM_POINTS; i++) {
        float theta = ((float)i / (float)NUM_POINTS) * SAMPLE_DURATION;

        signal[i][REAL] = 1.0 * cos(510.0 * 2*M_PI * theta) +  2.0 * cos(300.0 * 2*M_PI * theta);


        signal[i][IMAG] = 0;
    }
}

void do_something_with(fftw_complex* result) {
    float frequencydelta  = sampleFreq / NUM_POINTS;

    math::getScaledFrequencyBands(result, scaledFreqBand);
    scaledFreqBand[6][REAL] = 5;


    for (int i = 0; i < 6; i++) {
        float mag = sqrt(scaledFreqBand[i][REAL] * scaledFreqBand[i][REAL] +
                          scaledFreqBand[i][IMAG] * scaledFreqBand[i][IMAG]);

        magnitudes[i] = mag;
    }
    //math::normBuffer(magnitudes, NUM_POINTS/2);
    std::system("clear");

    for (int i = 0; i < 6; i++) {
        string out = "";
       
        out = to_string((int)(frequencydelta *     pow(2,i+1)))  ;
        //out = out + " " + to_string(magnitudes[i]);
        for(int a = 0; a < log10(magnitudes[i]) * 5; a ++)
        {
            out =  out + "#";
        }
        cout << out << endl;

    }
    cout << "frequency max: " << (sampleFreq/2) << endl;

}


/* Resume reading here */


int startMicRoutine()
{
    float* buffer = (float*) malloc(sizeof(float) * 512);
    initMicrophone(sampleFreq, NUM_POINTS);

    fftw_complex signal[NUM_POINTS];
    fftw_complex result[NUM_POINTS];

    for(int i = 0; i < NUM_POINTS; i ++)
    {
        signal[i][IMAG] = 0;
    }

    while(1){
        Pa_Sleep(1);
        if (!readCurrentBuffer(buffer))
        {
            
            fftw_plan plan = fftw_plan_dft_1d(NUM_POINTS,
                                            signal,
                                            result,
                                            FFTW_FORWARD,
                                            FFTW_ESTIMATE);

            for(int i = 0; i < NUM_POINTS; i ++)
            {
                signal[i][REAL] = buffer[i];
            }
            fftw_execute(plan);
            do_something_with(result);
            fftw_destroy_plan(plan);

        }
        //cout << "buffer miss" << endl;

    }
    

    releaseMicrophone();
}


int main() {
    startMicRoutine();   

    fftw_complex signal[NUM_POINTS];
    fftw_complex result[NUM_POINTS];

    fftw_plan plan = fftw_plan_dft_1d(NUM_POINTS,
                                      signal,
                                      result,
                                      FFTW_FORWARD,
                                      FFTW_ESTIMATE);

    acquire_from_somewhere(signal);
    fftw_execute(plan);
    do_something_with(result);

    fftw_destroy_plan(plan);


    return 0;
}