#pragma once

#include <string.h>
#include <stdio.h>
#include <iostream>

namespace math
{

    float getPeakAmplitude(float* signal, int bufferSize)
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
    

    
}