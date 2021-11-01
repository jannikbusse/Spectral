#pragma once

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <iostream>

namespace math
{

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

    void getScaledFrequencyBands(float rawSignal[], float (&sequencyBands)[6])
    {
        //this will be hardcoded for 64 size input and 6 size output.. this should be implemented better
        //the deltafrequency is 100 with these values (which is nice)
        sequencyBands[0] = 0;
        for(int i = 0; i <2; i++){
            sequencyBands[0] += rawSignal[i];
        }
        sequencyBands[1] = 0;
        for(int i = 2; i <4; i++){
            sequencyBands[1] += rawSignal[i];
        }
        sequencyBands[2] = 0;
        for(int i = 4; i <8; i++){
            sequencyBands[2] += rawSignal[i];
        }
        sequencyBands[3] = 0;
        for(int i = 8; i <16; i++){
            sequencyBands[3] += rawSignal[i];
        }
        sequencyBands[4] = 0;
        for(int i = 16; i <32; i++){
            sequencyBands[4] += rawSignal[i];
        }
        sequencyBands[5] = 0;
        for(int i = 32; i <64; i++){
            sequencyBands[5] += rawSignal[i];
        }
    }
    
    int applyCombFilter()
    {
        return 0;
    }
    
}