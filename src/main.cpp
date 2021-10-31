
#include <fftw3.h>
#include <iostream>
#include <utils/microphone.hpp>
#include <math/math.hpp>
#include <string.h>

#define NUM_POINTS 512
#define SAMPLE_DURATION 3.14


float sampleFreq = NUM_POINTS / SAMPLE_DURATION;


/* Never mind this bit */

#include <stdio.h>
#include <math.h>

#define REAL 0
#define IMAG 1

using namespace std;

void acquire_from_somewhere(fftw_complex* signal) {
    /* Generate two sine waves of different frequencies and
     * amplitudes.
     */

    int i;
    for (i = 0; i < NUM_POINTS; ++i) {
        double theta = (double)i / (double)NUM_POINTS * SAMPLE_DURATION;

        signal[i][REAL] = 1.0 * cos(10.0 * theta) +
                          0.5 * cos(25.0 * theta);

        cout << signal[i][REAL] << endl;

        signal[i][IMAG] = 0;
    }
}

void do_something_with(fftw_complex* result) {
    int i;
    float frequencydelta  = sampleFreq / NUM_POINTS;

    for (i = 0; i < NUM_POINTS / 8; ++i) {
        double mag = sqrt(result[i][REAL] * result[i][REAL] +
                          result[i][IMAG] * result[i][IMAG]);

        printf("%g: %g\n", (frequencydelta * i), mag);
    }
}


/* Resume reading here */

using namespace std;

int main() {
    float* buffer = (float*) malloc(sizeof(float) * 512);
    initMicrophone();
    while(1){
        Pa_Sleep(2);
        if (!readCurrentBuffer(buffer))
        {
            float peak = math::getPeakAmplitude(buffer, 512) * 20;
            string out = "";
            for(int i = 0; i < peak; i ++)
            {
                out = out + "#";
            }
            std::system("clear");
            cout << out << endl;


        }
    }
    releaseMicrophone();
    exit(0);

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

    cout << "frequency max: " << (sampleFreq/2) << endl;

    return 0;
}