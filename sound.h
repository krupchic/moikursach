#include <iostream>
#include <windows.h>
#include <mmsystem.h>
#include <cmath>
#pragma comment(lib,"winmm.lib") 
  
using namespace std;

static int quiet = 0;
static int bits = 16;
static int endian = 0;
static int raw = 0;
static int sign = 1;
const int NUMPTS = 44100 * 2 *0.5;   // 0,5 seconds
int sampleRate = 44100;
#define M_PI  3.14159265358979323846
#define PI	M_PI	
#define TWOPI	(2.0*PI)

void four1(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;
    
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}


double **
paintarray_gen (short int waveIn[NUMPTS])
{
	int i;
	int Nx;
	int NFFT;
	double *x;
	double *X;


	Nx = NUMPTS;

	x = (double *) malloc(Nx * sizeof(double));
	for(i=0; i<Nx; i++)
	{
		x[i] = waveIn[i];
	}

	NFFT = (int)pow(2.0, ceil(log((double)Nx)/log(2.0)));

	X = (double *) malloc((2*NFFT+1) * sizeof(double));


	for(i=0; i<Nx; i++)
	{
		X[2*i+1] = x[i];
		X[2*i+2] = 0.0;
	}

	for(i=Nx; i<NFFT; i++)
	{
		X[2*i+1] = 0.0;
		X[2*i+2] = 0.0;
	}
	double **paintarray = new double*[2];
	paintarray[0]=new double[NFFT+1];
	paintarray[1]=new double[NFFT+1];

	four1(X, NFFT, 1);
	double delta=((float)sampleRate)/(float)NFFT;
	double cur_freq=0;
	for (i=0;i<NFFT-2;i+=2)
	{
		cur_freq+=delta;
		paintarray[0][i/2]=cur_freq;
		paintarray[1][i/2]=(sqrt(X[i]*X[i]+X[i+1]*X[i+1]));
	}

  return paintarray;
}

double ** rec()
		{
			short int waveIn[NUMPTS];
		WAVEFORMATEX pFormat;
		pFormat.wFormatTag = WAVE_FORMAT_PCM;     // simple, uncompressed format
		pFormat.nChannels = 2;                    //  1=mono, 2=stereo
		pFormat.wBitsPerSample = 16;              //  16 for high quality, 8 for telephone-grade
		pFormat.nSamplesPerSec = sampleRate;     
		pFormat.nAvgBytesPerSec = sampleRate * pFormat.nChannels * pFormat.wBitsPerSample / 8;
		pFormat.nBlockAlign = pFormat.nChannels * pFormat.wBitsPerSample / 8;                 
		pFormat.cbSize = 0;
 
		HWAVEIN hWaveIn;
		WAVEHDR waveInHdr;

		waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT); 
 
		waveInHdr.lpData = (LPSTR)waveIn;
		waveInHdr.dwBufferLength = NUMPTS * 2;
		waveInHdr.dwBytesRecorded = 0;
		waveInHdr.dwUser = 0L;
		waveInHdr.dwFlags = 0L;
		waveInHdr.dwLoops = 0L;
 
		waveInPrepareHeader(hWaveIn, &waveInHdr, sizeof(WAVEHDR));
 
		waveInAddBuffer(hWaveIn, &waveInHdr, sizeof(WAVEHDR));
 

		if (waveInStart(hWaveIn)==MMSYSERR_NOERROR) 
			printf("no error1\n");
 
		// Wait until finished recording
		do{} 
		while(waveInUnprepareHeader(hWaveIn, &waveInHdr, sizeof(WAVEHDR)) == WAVERR_STILLPLAYING);
 
		waveInClose(hWaveIn);

		return paintarray_gen(waveIn);
		}