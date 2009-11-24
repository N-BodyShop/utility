#ifndef HILBERT_HINCLUDED
#define HILBERT_HINCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

uint64_t hilbert2d(float x,float y);
uint64_t hilbert3d(float x,float y,float z);
void ihilbert3d(uint64_t s,float *px,float *py,float *pz);
void ihilbert2d(uint64_t s,float *px,float *py,float *pz);

#if defined(__cplusplus)
}
#endif

#endif
