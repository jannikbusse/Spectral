/*
    * $Id$
    *
    * This program uses the PortAudio Portable Audio Library.
    * For more information see: http://www.portaudio.com
    * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
    *
    * Permission is hereby granted, free of charge, to any person obtaining
    * a copy of this software and associated documentation files
    * (the "Software"), to deal in the Software without restriction,
    * including without limitation the rights to use, copy, modify, merge,
    * publish, distribute, sublicense, and/or sell copies of the Software,
    * and to permit persons to whom the Software is furnished to do so,
    * subject to the following conditions:
    *
    * The above copyright notice and this permission notice shall be
    * included in all copies or substantial portions of the Software.
    *
    * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
    * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    */
   
   /*
    * The text above constitutes the entire PortAudio license; however, 
    * the PortAudio community also makes the following non-binding requests:
    *
    * Any person wishing to distribute modifications to the Software is
    * requested to send the modifications to the original developer so that
    * they can be incorporated into the canonical version. It is also 
    * requested that these non-binding requests be included along with the 
    * license above.
    */
   
   #include <stdio.h>
   #include <stdlib.h>
   #include "portaudio.h"
   #include <iostream>
   
   /* #define SAMPLE_RATE  (17932) // Test failure to open with this value. */
   #define SAMPLE_RATE  (44100)
   #define FRAMES_PER_BUFFER (512)
   #define NUM_SECONDS     (5)
   #define NUM_CHANNELS    (1)
   /* #define DITHER_FLAG     (paDitherOff) */
   #define DITHER_FLAG     (0) 
   
   
   /* Select sample format. */
   #define PA_SAMPLE_TYPE  paFloat32
   typedef float SAMPLE;
   #define SAMPLE_SILENCE  (0.0f)
   #define PRINTF_S_FORMAT "%.8f"
   
   
   typedef struct
   {
       int          frameIndex;  /* Index into sample array. */
       int          maxFrameIndex;
       SAMPLE      *recordedSamples;
   }
   paTestData;

    PaStreamParameters  inputParameters;
    PaStream*           stream;
    PaError             err = paNoError;
    paTestData          bufferdata;
    int                 i;
    int                 totalFrames;
    int                 numSamples;
    int                 numBytes;
    SAMPLE              maxframes, val;
    double              average;

    int sleepCycles = 0;

   
   /* This routine will be called by the PortAudio engine when audio is needed.
   ** It may be called at interrupt level on some machines so don't do anything
   ** that could mess up the system like calling malloc() or free().
   */
   static int recordCallback( const void *inputBuffer, void *outputBuffer,
                              unsigned long framesPerBuffer,
                              const PaStreamCallbackTimeInfo* timeInfo,
                              PaStreamCallbackFlags statusFlags,
                              void *userData )
   {
       paTestData *bufferdata = (paTestData*)userData;
       const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
       SAMPLE *wptr = &bufferdata->recordedSamples[bufferdata->frameIndex * NUM_CHANNELS];
       long framesToCalc;
       long i;
       int finished;
       unsigned long framesLeft = bufferdata->maxFrameIndex - bufferdata->frameIndex;
   
       (void) outputBuffer; /* Prevent unused variable warnings. */
       (void) timeInfo;
       (void) statusFlags;
       (void) userData;
   
       if( framesLeft < framesPerBuffer )
       {
           framesToCalc = framesLeft;
           finished = paComplete;
       }
       else
       {
           framesToCalc = framesPerBuffer;
           finished = paContinue;
       }
   
       if( inputBuffer == NULL )
       {
           for( i=0; i<framesToCalc; i++ )
           {
               *wptr++ = SAMPLE_SILENCE;  /* left */
               if( NUM_CHANNELS == 2 ) *wptr++ = SAMPLE_SILENCE;  /* right */
           }
       }
       else
       {
           for( i=0; i<framesToCalc; i++ )
           {
               *wptr++ = *rptr++;  /* left */
               if( NUM_CHANNELS == 2 ) *wptr++ = *rptr++;  /* right */
           }
       }
       bufferdata->frameIndex += framesToCalc;
       return finished;
   }

   int initMicrophone()
   {
       return 0;
   }

   int releaseMicrophone()
   {
       return 0;
   }
   
  
   paTestData psudomain(void)
   {
      
   
       bufferdata.maxFrameIndex = FRAMES_PER_BUFFER; /* Record for a few seconds. */
       bufferdata.frameIndex = 0;
       bufferdata.recordedSamples = (SAMPLE *) malloc( FRAMES_PER_BUFFER * sizeof(SAMPLE)); /* From now on, recordedSamples is initialised. */
      
       for( i=0; i<numSamples; i++ ) bufferdata.recordedSamples[i] = 0;
   
       err = Pa_Initialize();
   
       inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
       if (inputParameters.device == paNoDevice) {
           fprintf(stderr,"Error: No default input device.\n");
           goto done;
       }
       inputParameters.channelCount = 1;                    /* stereo input */
       inputParameters.sampleFormat = PA_SAMPLE_TYPE;
       inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency;
       inputParameters.hostApiSpecificStreamInfo = NULL;
   
       /* Record some audio. -------------------------------------------- */
       err = Pa_OpenStream(
                 &stream,
                 &inputParameters,
                 NULL,                  /* &outputParameters, */
                 SAMPLE_RATE,
                 FRAMES_PER_BUFFER,
                 paClipOff,      /* we won't output out of range samples so don't bother clipping them */
                 NULL,
                 &bufferdata );
   
       err = Pa_StartStream( stream );
       if( err != paNoError ) goto done;
       printf("\n=== Now recording!! Please speak into the microphone. ===\n"); fflush(stdout);
   
       while( ( err = Pa_IsStreamActive( stream ) ) == 1 )
       {
             Pa_Sleep(1);
        //    printf("index = %d\n", bufferdata.frameIndex ); fflush(stdout);
            unsigned long availableData =  Pa_GetStreamReadAvailable(stream);
            if (availableData <= FRAMES_PER_BUFFER) 
            {   
                sleepCycles ++;
                continue;
            }
           std::cout << sleepCycles << std::endl;
           sleepCycles = 0;

            err = Pa_ReadStream(stream, bufferdata.recordedSamples,FRAMES_PER_BUFFER);




          
          
        //    std::cout << bufferdata.recordedSamples[bufferdata.frameIndex - 1] << std::endl;
           //printf("index = %g\n", bufferdata.recordedSamples[bufferdata.frameIndex] ); fflush(stdout);
           
       }
       if( err < 0 ) goto done;
   
       err = Pa_CloseStream( stream );
       
     
   done:
       Pa_Terminate();
       if( bufferdata.recordedSamples )       /* Sure it is NULL or valid. */
           free( bufferdata.recordedSamples );
       if( err != paNoError )
       {
           fprintf( stderr, "An error occured while using the portaudio stream\n" );
           fprintf( stderr, "Error number: %d\n", err );
           fprintf( stderr, "Error message: %s\n", Pa_GetErrorText( err ) );
           err = 1;          /* Always return 0 or 1, but no other return codes. */
       }
       return bufferdata;
   }
   