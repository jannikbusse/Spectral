   #include <stdio.h>
   #include <stdlib.h>
   #include "portaudio.h"
   #include <iostream>
   
   /* #define SAMPLE_RATE  (17932) // Test failure to open with this value. */
   unsigned int SAMPLE_RATE = 44100;
   unsigned int FRAMES_PER_BUFFER = 512;

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

   int initMicrophone(unsigned int sample_rate = 44100, unsigned int buffer_size = 512)
   {
       SAMPLE_RATE = sample_rate;
       FRAMES_PER_BUFFER = buffer_size;

    bufferdata.maxFrameIndex = FRAMES_PER_BUFFER; /* Record for a few seconds. */
       bufferdata.frameIndex = 0;
       bufferdata.recordedSamples = (SAMPLE *) malloc( FRAMES_PER_BUFFER * sizeof(SAMPLE)); /* From now on, recordedSamples is initialised. */
      
       for( i=0; i<numSamples; i++ ) bufferdata.recordedSamples[i] = 0;
   
       err = Pa_Initialize();
   
       inputParameters.device = Pa_GetDefaultInputDevice(); /* default input device */
       if (inputParameters.device == paNoDevice) {
           fprintf(stderr,"Error: No default input device.\n");
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

       return 0;
   }

   int releaseMicrophone()
   {

    err = Pa_CloseStream( stream );
    
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

    return 0;
   }
   
  
   int readCurrentBuffer(float *buffer)
   {
   
        if( ( err = Pa_IsStreamActive( stream ) ) == 0 )
        {
            std::cout << "ERR: audio stream is not open" << std::endl;
            return 1;
        }
        unsigned long availableData =  Pa_GetStreamReadAvailable(stream);
        if(availableData >= FRAMES_PER_BUFFER)
        {
            std::cout <<"buffer: " << Pa_GetStreamReadAvailable(stream) << std::endl;

           
            while(Pa_GetStreamReadAvailable(stream)>= FRAMES_PER_BUFFER )
                 err = Pa_ReadStream(stream, buffer,FRAMES_PER_BUFFER);
            

        }
        else
            return 1;
           
       return 0;
   
   }
   