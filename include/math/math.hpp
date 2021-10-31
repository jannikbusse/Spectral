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

    void getScaledFrequencyBands(fftw_complex rawSignal[], float (&sequencyBands)[6][2])
    {
        //this will be hardcoded for 32 size input and 6 size output.. this should be implemented better
        //the deltafrequency is 100 with these values (which is nice)
            for(int a = 0; a < 2; a++)
                {
                sequencyBands[0][a] = 0;
                for(int i = 0; i <2; i++){
                    sequencyBands[0][a] += rawSignal[i][a];
                }
                sequencyBands[1][a] = 0;
                for(int i = 2; i <4; i++){
                    sequencyBands[1][a] += rawSignal[i][a];
                }
                sequencyBands[2][a] = 0;
                for(int i = 4; i <8; i++){
                    sequencyBands[2][a] += rawSignal[i][a];
                }
                sequencyBands[3][a] = 0;
                for(int i = 8; i <16; i++){
                    sequencyBands[3][a] += rawSignal[i][a];
                }
                sequencyBands[4][a] = 0;
                for(int i = 16; i <32; i++){
                    sequencyBands[4][a] += rawSignal[i][a];
                }
                sequencyBands[5][a] = 0;
                for(int i = 32; i <50; i++){
                    sequencyBands[5][a] += rawSignal[i][a];
                }
            }            
    }
    

    
}