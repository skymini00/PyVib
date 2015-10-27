#include "extcode.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef uint16_t  OCTSetup;
#define OCTSetup_SS200Khz 0
#define OCTSetup_SS50KHz 1
#define OCTSetup_SD26Khz 2
#define OCTSetup_SS200KHzNew 3
#define OCTSetup_SS50KHzNew 4
typedef struct {
	uint16_t Ch0Shift;
	uint16_t numKlinPts;
	uint16_t StartTriggerOfffset;
	int16_t PDScaling;
	LVBoolean BGSubtract;
	LVBoolean DispCorr;
	LVBoolean InterpData;
	LVBoolean FFT;
	uint16_t procRoiBegin;
	uint16_t procRoiEnd;
	LVBoolean InterpFilter;
	LVBoolean InterpDownsample;
	uint16_t klinRoiBegin;
	uint16_t klinRoiEnd;
	int16_t SampleOffset;
	int16_t MZIScaleFactor;
	LVBoolean MZIHPFilter;
	uint16_t DownsampleFactor;
	int16_t SamplesPerTrig;
	int16_t Ch0OffsetTweak;
	LVBoolean SyncTrig;
	LVBoolean ProcessMZI;
	LVBoolean FFT16b;
	LVBoolean Polar;
	LVBoolean MagOnly;
	int8_t FFTReImgScale;
	LVBoolean InterpDSAvg;
	uint16_t OnDutyCycleTrigLast;
	uint16_t OffDutyCycleTrigLast;
	LVBoolean UseDutyCycle;
} FPGASSOCTAcquistionOpts;

/*!
 * InitFPGA
 */
int32_t __cdecl InitFPGA(OCTSetup OCTSetup, 
	FPGASSOCTAcquistionOpts *FPGAProcessingOpts, int32_t klinOut[], int32_t *len);
/*!
 * ConfigureFPGAAcquisition
 */
int32_t __cdecl ConfigureFPGAAcquisition(OCTSetup OCTSetup, 
	uint32_t NumTriggers, FPGASSOCTAcquistionOpts *FPGAProcessingOpts, 
	int32_t *ROIOut, int32_t *NumTrigsOut, 
	FPGASSOCTAcquistionOpts *FPGAProcessingOptsOut);
/*!
 * AcquireFPGAData
 */
int32_t __cdecl AcquireFPGAData(OCTSetup OCTSetup, uint32_t NumTriggers, 
	uint32_t NumSamples, uint32_t TriggerOffset, LVBoolean UnpackData, 
	float data_re_out[], float data_im_out[], int32_t *data_re_im_len, 
	uint64_t PackedDataOut[], int32_t *packed_data_len, uint32_t *TimeElapsedMs, 
	uint32_t *TransferTimeMs, uint32_t *UnpackTimeMs);
/*!
 * CloseFPGA
 */
int32_t __cdecl CloseFPGA(OCTSetup OCTSetup);
/*!
 * DAQIOWaitUntilDone
 */
int32_t __cdecl DAQIOWaitUntilDone(uintptr_t *taskChannelsIn, 
	double timeoutSec);
/*!
 * DAQIOXYScan
 */
int32_t __cdecl DAQIOXYScan(char ScanXChannel[], char ScanYChannel[], 
	char TrigChannel[], double outputRate, double XCommand[], double YCommand[], 
	int32_t x_len, int32_t y_len, uintptr_t *taskOut);
/*!
 * LoadDispersionFPGA
 */
int32_t __cdecl LoadDispersionFPGA(OCTSetup OCTSetup, double MagWindow[], 
	double PhaseCorr[], int32_t len);
/*!
 * ResetFPGA
 */
int32_t __cdecl ResetFPGA(OCTSetup OCTSetup, 
	FPGASSOCTAcquistionOpts *FPGAProcessingOpts, int32_t klinOut[], int32_t *len);
/*!
 * DAQIOAudioTwoSpeakersMicInput
 */
int32_t __cdecl DAQIOAudioTwoSpeakersMicInput(double OutputRate, 
	double InputRate, char MicChan[], char TrigSrcIn[], char TrigSrcOut[], 
	char RightSpkChan[], char LeftSpkChan[], double RightOutput[], 
	double LeftOutput[], uint32_t *numMicSamples, uintptr_t *inputTaskOut, 
	uintptr_t *outputTaskOut, int32_t rightOutput_len, int32_t leftOutput_len);
/*!
 * DAQIOSetAttenuatorLevel
 */
int32_t __cdecl DAQIOSetAttenuatorLevel(uint8_t AttenLevelDB, 
	char AttenuatorLines[]);
/*!
 * DAQIOSetMirrorPostion
 */
int32_t __cdecl DAQIOSetMirrorPostion(double XCommand, double YCommand, 
	char XChan[], char YChan[]);
/*!
 * DAQIOAudioStartOutputAndGrabMicData
 */
int32_t __cdecl DAQIOAudioStartOutputAndGrabMicData(uintptr_t *micTaskIn, 
	char TrigOutLine[], uintptr_t *speakerTaskIn, int32_t numMicSamples, 
	double MicDataOut[]);
/*!
 * DAQIOCamTriggerOnly
 */
int32_t __cdecl DAQIOCamTriggerOnly(double camTrigCommand[], int32_t len, 
	double OutputRate, char CamTrigChan[], char DIOTrigIn[], 
	uintptr_t *daqTaskOut);
/*!
 * DAQIOSendTrigger
 */
int32_t __cdecl DAQIOSendTrigger(char TrigOutLine[]);
/*!
 * DAQIOXYScanWithCamTrig
 */
int32_t __cdecl DAQIOXYScanWithCamTrig(char ScanXChannel[], 
	char ScanYChannel[], char CamTrigCh[], double XCommand[], double YCommand[], 
	double CamTrigCommand[], char DAQTrigChannel[], double outputRate, 
	int32_t len, uintptr_t *taskOut);
/*!
 * DAQReset
 */
int32_t __cdecl DAQReset(char DeviceName[]);
/*!
 * DAQIOAudioGrabMicData
 */
int32_t __cdecl DAQIOAudioGrabMicData(uintptr_t *micTaskIn, double timeout, 
	int32_t numMicSamples, double MicDataOut[], int32_t *len);
/*!
 * RecalcKlin
 */
int32_t __cdecl RecalcKlin(OCTSetup OCTSetup, 
	FPGASSOCTAcquistionOpts *FPGAProcessingOpts, int32_t klinOut[], int32_t *len);

MgErr __cdecl LVDLLStatus(char *errStr, int errStrLen, void *module);

#ifdef __cplusplus
} // extern "C"
#endif

