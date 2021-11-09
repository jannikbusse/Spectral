#pragma once

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

namespace math
{

    // float SAMPLE_POINTS = 1024;
    // float SAMPLE_TIME = 0.015;
    // float SAMPLE_FREQUENCY = SAMPLE_POINTS / SAMPLE_TIME;
    // float FREQUENCY_DELTA = SAMPLE_FREQUENCY / SAMPLE_POINTS;

    float getPeakAmplitude(const float* signal, int bufferSize)
    {
        float peak = 0;
        for(int i = 0; i < bufferSize; i ++)
        {
            if(abs(signal[i]) > peak)
            {
                peak = abs(signal[i]);
            }
        }
        return peak;
    }

    void normBuffer(float *buffer, int bufferSize)
    {
        double peak = 0;
        for(int i = 0; i < bufferSize; i ++)
        {
            if(abs(buffer[i]) > peak)
            {
                peak = abs(buffer[i]);
            }
        }
        for(int i = 0; i < bufferSize; i ++)
        {
            buffer[i] /= peak;
        }
    }

    void getScaledFrequencyBands(float rawSignal[], float (&sequencyBands)[6], float FREQUENCY_DELTA)
    {
        //this will be hardcoded for 512 size input and 6 size output.. this should be implemented better while having 0.04 sec time window
        //0 - 200
        //200 - 400
        //400 - 800
        //800 - 1600
        //1600 - 3200
        //3200 - Rest

        sequencyBands[0] = 0;
        for(int i = 1; i < 200.0/FREQUENCY_DELTA; i++){
            sequencyBands[0] += rawSignal[i];
        }
        sequencyBands[1] = 0;
        for(int i = 200.0/FREQUENCY_DELTA; i <400.0/FREQUENCY_DELTA; i++){
            sequencyBands[1] += rawSignal[i];
        }
        sequencyBands[2] = 0;
        for(int i = 400.0/FREQUENCY_DELTA; i <800.0/FREQUENCY_DELTA; i++){
            sequencyBands[2] += rawSignal[i];
        }
        sequencyBands[3] = 0;
        for(int i = 800.0/FREQUENCY_DELTA; i <1600.0/FREQUENCY_DELTA; i++){
            sequencyBands[3] += rawSignal[i];
        }
        sequencyBands[4] = 0;
        for(int i = 1600.0/FREQUENCY_DELTA; i <3200.0/FREQUENCY_DELTA; i++){
            sequencyBands[4] += rawSignal[i];
        }
        sequencyBands[5] = 0;
        for(int i = 3200.0/FREQUENCY_DELTA; i <5000.0/FREQUENCY_DELTA; i++){
            sequencyBands[5] += rawSignal[i];
        }
    }

    float applyHanningFunction(int curId, int maxId)
    {
        return 0.5f * ( 1 - cos((2* M_PI * (float) curId)/((float)maxId - 1.0f)));
    }
    
    int applyCombFilter()
    {
        return 0;
    }
    
}