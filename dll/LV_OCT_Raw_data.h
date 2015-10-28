#include "extcode.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef uint16_t  OCTSetupRaw;
#define OCTSetupRaw_SS50KHz 0
#define OCTSetupRaw_SS200KHz 1

/*!
 * InitFPGARaw
 */
int32_t __cdecl InitFPGARaw(OCTSetupRaw OCTSetup);
/*!
 * ConfigureFPGAAcquisitionRaw
 */
int32_t __cdecl ConfigureFPGAAcquisitionRaw(uint32_t NumTriggers, 
	uint16_t SamplesPerTrigger);
/*!
 * AcquireFPGADataRaw
 */
int32_t __cdecl AcquireFPGADataRaw(uint32_t NumTriggers, int16_t ch0Data[], 
	int16_t ch1Data[], int32_t data_len_in, int32_t *data_len_out, 
	uint32_t *SamplesRemaining, uint32_t *TimeElapsedMs);
/*!
 * CloseFPGARaw
 */
int32_t __cdecl CloseFPGARaw(void);

MgErr __cdecl LVDLLStatus(char *errStr, int errStrLen, void *module);

#ifdef __cplusplus
} // extern "C"
#endif

