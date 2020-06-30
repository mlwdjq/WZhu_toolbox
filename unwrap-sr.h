//-------------------------------------------------------------------------------
// �ļ���unwrap-sr.h
//
// ��������λ��� -- ���ڿɿ�������Ķ�ά��λ��� SR (Sorting Reliability)
//       �ο����� "Fast two-dimensional phase-unwrapping algorithm based on 
//       sorting by reliability following a noncontinuous path"
//
// ���ߣ����
//
// �޸ģ�2009��03��21�� ��һ��
//-------------------------------------------------------------------------------

#ifndef __UNWRAP_SR_H_
#define __UNWRAP_SR_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <malloc.h>
#include <string.h>

#if (defined(__cplusplus) || defined(c_plusplus))
extern "C" {
#endif

// ����ʱ������ʾ������Ϣ
#define __UNWRAP_SR_DEBUG_ENABLE__
//#define __UNWRAP_SR_DLL_ENABLE__

#ifdef __UNWRAP_SR_DLL_ENABLE__
#define UNWRAP_SR_API __declspec(dllexport)
#else
#define UNWRAP_SR_API
#endif

#define UNWRAP_SR_SUCCESS					0	// �������óɹ�
#define	UNWRAP_SR_ERROR_PARAM_PHASE			-1  // ���������λ������Ч
#define	UNWRAP_SR_ERROR_PARAM_LO_HI_VAL		-2  // ���������λ���ݵ�����ֵ������ֵ��Ч
#define	UNWRAP_SR_ERROR_PARAM_HEIGHT		-3  // ��λ���ݵ�������Ч 
#define	UNWRAP_SR_ERROR_PARAM_WIDTH			-4  // ��λ���ݵ�������Ч 
#define	UNWRAP_SR_ERROR_MEMORY				-5  // �����ڴ���Ч
#define	UNWRAP_SR_ERROR_PARAM_MASK			-6  // ��λ���ݵ���ģ������Ч 

#define UNWRAP_SR_PI			3.141592653589793
#define UNWRAP_SR_TWO_PI		6.28318530717959
#define UNWRAP_SR_SWAP(x,y)		{UNWRAP_SR_EDGE t; t=x; x=y; y=t;}
#define UNWRAP_SR_ORDER(x,y)	if (x.reliab > y.reliab) UNWRAP_SR_SWAP(x,y)
#define UNWRAP_SR_O2(x,y)		UNWRAP_SR_ORDER(x,y)
#define UNWRAP_SR_O3(x,y,z)		UNWRAP_SR_O2(x,y); UNWRAP_SR_O2(x,z); UNWRAP_SR_O2(y,z)
typedef enum {UNWRAP_SR_YES, UNWRAP_SR_NO} UNWRAP_SR_YES_NO;


// UNWRAP_SR_PIXELM information
typedef struct PixElmStruct
{
    int increment;				// No. of 2*pi to add to the pixel to unwrap it
    int PixelNumInGroup;		// No. of pixel in the pixel group
    double value;				// value of the pixel
	double reliability;	
	unsigned char InputMask;	// 0 pixel is masked. 255 pixel is not masked
	unsigned char ExtendedMask;	// 0 pixel is masked. 255 pixel is not masked
    int group;					// group No.
    int NewGroup;
    struct PixElmStruct *head;	// pointer to the first pixel in the group in the linked list
    struct PixElmStruct *last;	// pointer to the last pixel in the group
    struct PixElmStruct *next;	// pointer to the next pixel in the group
} UNWRAP_SR_PIXELM;

// the UNWRAP_SR_EDGE is the line that connects two pixels.
// if we have S pixels, then we have S horizontal edges and S vertical edges
typedef struct 
{    
	double reliab;				// reliabilty of the edge and it depends on the two pixels
	UNWRAP_SR_PIXELM *pointer1;	// pointer to the first pixel
    UNWRAP_SR_PIXELM *pointer2;	// pointer to the second pixel
    int increment;				// No. of 2*pi to add to one of the pixels to unwrap it with respect to the second 
} UNWRAP_SR_EDGE; 


//-------------------------------------------------------------------------------
// UnwrapPhaseBySortingReliabilityWithMask()
// ����������Ƿ�����ͼ��λ���ݣ���Ҫ��ģͼ��
// ������phase -- ���������λ���ݣ������������Ϊ�ѽ��������
//       LowVal -- ���������λ���ݵ�����ֵ
//       HighVal -- ���������λ���ݵ�����ֵ
//       mask -- ��λ���ݵ���ģ����
//       ValidPoint -- ��λ��ģ��Ч����ı��ֵ
//       height -- ��λ���ݵ�����
//       width -- ��λ���ݵ�����
// ����ֵ�����ΪUNWRAP_SR_SUCCESS����ʾ�ú�������ȷ���á����Ϊ��ֵ����ʾ�ú�
//         ������ʱ���ִ��󣬳��ִ�����������ͨ�������ֵ���жϣ��������
//         UNWRAP�������ĺ궨�塣
//-------------------------------------------------------------------------------
UNWRAP_SR_API int UnwrapPhaseBySortingReliabilityWithMask(double *phase, double LowVal, double HighVal, 
										                  unsigned char *mask, int ValidPoint, int height, 
										                  int width);


#if (defined(__cplusplus) || defined(c_plusplus))
}
#endif

#endif
