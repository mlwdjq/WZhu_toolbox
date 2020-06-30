//-------------------------------------------------------------------------------
// 文件：unwrap-sr.h
//
// 描述：相位解包 -- 基于可靠性排序的二维相位解包 SR (Sorting Reliability)
//       参考文献 "Fast two-dimensional phase-unwrapping algorithm based on 
//       sorting by reliability following a noncontinuous path"
//
// 作者：李建欣
//
// 修改：2009年03月21日 第一版
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

// 调试时用来显示错误信息
#define __UNWRAP_SR_DEBUG_ENABLE__
//#define __UNWRAP_SR_DLL_ENABLE__

#ifdef __UNWRAP_SR_DLL_ENABLE__
#define UNWRAP_SR_API __declspec(dllexport)
#else
#define UNWRAP_SR_API
#endif

#define UNWRAP_SR_SUCCESS					0	// 函数调用成功
#define	UNWRAP_SR_ERROR_PARAM_PHASE			-1  // 待解包的相位数据无效
#define	UNWRAP_SR_ERROR_PARAM_LO_HI_VAL		-2  // 待解包的相位数据的下限值或上限值无效
#define	UNWRAP_SR_ERROR_PARAM_HEIGHT		-3  // 相位数据的行数无效 
#define	UNWRAP_SR_ERROR_PARAM_WIDTH			-4  // 相位数据的列数无效 
#define	UNWRAP_SR_ERROR_MEMORY				-5  // 分配内存无效
#define	UNWRAP_SR_ERROR_PARAM_MASK			-6  // 相位数据的掩模数据无效 

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
// 描述：处理非方形整图相位数据，需要掩模图像
// 参数：phase -- 待解包的相位数据，函数处理完后为已解包的数据
//       LowVal -- 待解包的相位数据的下限值
//       HighVal -- 待解包的相位数据的上限值
//       mask -- 相位数据的掩模数据
//       ValidPoint -- 相位掩模有效区域的标记值
//       height -- 相位数据的行数
//       width -- 相位数据的列数
// 返回值：如果为UNWRAP_SR_SUCCESS，表示该函数已正确调用。如果为负值，表示该函
//         数调用时出现错误，出现错误的情况可以通过具体的值来判断，见上面的
//         UNWRAP错误代码的宏定义。
//-------------------------------------------------------------------------------
UNWRAP_SR_API int UnwrapPhaseBySortingReliabilityWithMask(double *phase, double LowVal, double HighVal, 
										                  unsigned char *mask, int ValidPoint, int height, 
										                  int width);


#if (defined(__cplusplus) || defined(c_plusplus))
}
#endif

#endif
