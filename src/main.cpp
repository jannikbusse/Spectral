
#include <mprocessor/music.h>

#include <fftw3.h>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <thread>
#include "ws2811/ws2811.h"

#include <map>
#include <vector>
#include <cmath>
#include <bitset>

using namespace std;

#define REAL 0
#define IMAG 1


//led---------------

// defaults for cmdline options
#define TARGET_FREQ             WS2811_TARGET_FREQ
#define GPIO_PIN                18
#define DMA                     10
//#define STRIP_TYPE            WS2811_STRIP_RGB		// WS2812/SK6812RGB integrated chip+leds
#define STRIP_TYPE              WS2811_STRIP_GRB		// WS2812/SK6812RGB integrated chip+leds
//#define STRIP_TYPE            SK6812_STRIP_RGBW		// SK6812RGBW (NOT SK6812RGB)

#define LED_COUNT                   150

int led_count = LED_COUNT;

int clear_on_exit = 0;

ws2811_t ledstring =
{
    .freq = TARGET_FREQ,
    .dmanum = DMA,
    .channel =
    {
        [0] =
        {
            .gpionum = GPIO_PIN,
            .invert = 0,
	    .count = LED_COUNT,
	    .strip_type = STRIP_TYPE,
            .brightness = 255,
        },
        [1] =
        {
            .gpionum = 0,
            .invert = 0,
            .count = 0,
            .brightness = 0,
        },
    },
};

//---------------
uint32_t rgb2hex(uint8_t r, uint8_t g, uint8_t b)
{
    uint32_t res = b |g << 8 | r << 16 ;
    return res;
}

void basscallback()
{
    cout << "bass detected" << endl;
}


int main()
{

    // ws2811_init(&ledstring);
    // for(int i = 0; i < 80;i ++)
    //     ledstring.channel[0].leds[i] = rgb2hex(0,0,255);
    
//    ws2811_render(&ledstring);

  //  return 0;
    
    init(&basscallback);  
    while(1)
        update();

    return 0;
}
