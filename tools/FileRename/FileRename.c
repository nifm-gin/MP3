// FileRename.c
// Rename file or folder
// This function renames the existing file or folder specified by the string
// Source to the name given by the string Dest. You can use FileRename to move
// a file from one folder to another folder or drive, but folders can be renamed
// only, not moved.
//
// Files and folders can be renamed by Matlab's MOVEFILE also, but this C-Mex is
// faster (timings vary with the size and number of the files due to the
// caching of write operations by the hard disk and the OS):
//    Matlab 2009a: 4 to 50 times faster,
//    Matlab 6.5:   1600 times faster (!).
// The fast C-Mex file must be compiled before using.
//
// [Status, Msg] = FileRename(Source, Dest, [Mode])
// INPUT:
//   Source: String, name of the source file or folder.
//           Unicode and UNC paths are considered.
//   Dest:   String, name of the destination file or folder.
//   Mode:   String, if 'forced' an existing Dest file is overwritten,
//           if it is not write protected. Folders are *not* overwritten.
//           Optional, default: 'DoNotOverwrite'.
//
// OUTPUT:
//   Status: Scalar DOUBLE. Optional.
//            0: Success
//           -1: Source is not existing
//           -2: Dest is existing already
//           -3: Dest is write protected, in forced [Mode] only
//           -4: Unknown problems:
//               Source or Dest is accessed from another program,
//               Source is a folder and Dest is on another drive.
//   Msg: String, empty on success, some information in case of problems.
//
// COMPILE:
//   mex -O FileRename.c
// Linux: consider C99 comments:
//   mex -O CFLAGS="\$CFLAGS -std=c99" FileRename.c
// This function cannot be compiled with LCC2.4 shipped with Matlab 6.5,
// but with LCC included in Matlab 2009a.
// Pre-compiled Mex: http://www.n-simon.de/mex
//
// NOTE: The behaviour of the underlying C function _wrename() depends on the
//   implementation, if the destination is existing. Run the included test
//   function uTest_FileRename after compiling to check if the function works
//   as exspected.
//
// Tested: Matlab 6.5, 7.7, 7.8, WinXP, 32bit
//         Compiler: LCC (Matlab 2009a), OWC1.8, BCC5.5, MSVC2008
// Assumed Compatibility: higher Matlab versions, Mac, Linux, 64bit
// Author: Jan Simon, Heidelberg, (C) 2006-2010 matlab.THISYEAR(a)nMINUSsimon.de

/*
% $JRev: R0d V:003 Sum:xgfOguzq3o0q Date:28-Nov-2010 03:31:30 $
% $License: BSD $
% $UnitTest: uTest_FileRename $
% $File: Tools\Mex\Source\FileRename.c $
% History:
% 001: 27-Nov-2010 00:12, First version.
*/

#include <wchar.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"

#if defined(__BORLANDC__)  //  needed for _wrename
#include <io.h>
#endif

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_ID   "JSimon:FileRename:"
#define ERR_HEAD "*** FileRename[mex]: "

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  wchar_t    *Source, *Dest, ModeChar;
  mwSize     SourceLen, DestLen;
  int        result, forced = 0;
  const char *ErrMsg;
  
  // Check number and type of arguments: ---------------------------------------
  if (nrhs != 2 && nrhs != 3) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNInput",
                       ERR_HEAD "2 or 3 inputs required.");
  }
  if (nlhs > 2) {
     mexErrMsgIdAndTxt(ERR_ID   "BadNOutput",
                       ERR_HEAD "2 outputs allowed.");
  }
  
  // Type of input arguments:
  if (!mxIsChar(prhs[0]) || !mxIsChar(prhs[1])) {
     mexErrMsgIdAndTxt(ERR_ID   "BadInputType",
                       ERR_HEAD "Inputs must be strings.");
  }
  
  // Parse 3rd input - check first character only:
  if (nrhs == 3) {
     if (!mxIsChar(prhs[2])) {
        mexErrMsgIdAndTxt(ERR_ID   "BadTypeInput3",
                          ERR_HEAD "3rd input [Mode] must be a string.");
     }
     
     if (!mxIsEmpty(prhs[2])) {
        ModeChar = *(mxChar *) mxGetData(prhs[2]);
        forced   = (ModeChar == L'f') || (ModeChar == L'F');
     }
  }
  
  // Obtain names as Unicode strings: ------------------------------------------
  SourceLen = mxGetNumberOfElements(prhs[0]);
  Source    = (wchar_t *) mxMalloc((SourceLen + 1) * sizeof(mxChar));
  if (Source == NULL) {
     mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                       ERR_HEAD "Cannot get memory for file name.");
  }
  memcpy(Source, mxGetData(prhs[0]), SourceLen * sizeof(mxChar));
  Source[SourceLen] = L'\0';
  
  DestLen = mxGetNumberOfElements(prhs[1]);
  Dest    = (wchar_t *) mxMalloc((DestLen + 1) * sizeof(mxChar));
  if (Dest == NULL) {
     mexErrMsgIdAndTxt(ERR_ID   "NoMemory",
                       ERR_HEAD "Cannot get memory for file name.");
  }
  memcpy(Dest, mxGetData(prhs[1]), DestLen * sizeof(mxChar));
  Dest[DestLen] = L'\0';
  
  // Renaming (some compilers reply [errno_t], some [int]): --------------------
  // Deleting Dest in forced [Mode] is implemented in the problem handling.
  result = (int) _wrename(Source, Dest);
  
  // Handle problems: ----------------------------------------------------------
  if (result == 0) {                       // Success:
     if (nlhs > 0) {
        plhs[0] = mxCreateDoubleScalar(1.0);
     
        if (nlhs > 1) {
           plhs[1] = mxCreateString("");
        }
     }
     
  } else {                                 // Problems:
     if (_waccess(Source, 0) != 0) {       // Source is not existing:
        result = -1;
     } else if (_waccess(Dest, 0) == 0) {  // Dest is existing already:
        // Try to delete Dest for forced [Mode]:
        if (forced) {                      // Forced [Mode]:
           result = (int) _wremove(Dest);  // Delete Dest file (not folder!)
           if (result == 0) {              // Dest was deleted:
              result = (int) _wrename(Source, Dest);
              if (result != 0) {
                 result = -4;              // Source is blocked?!
              }
           } else {                        // Dest is write protected:
              result = -3;
           }
        } else {
          result = -2;
        }
     } else {                              // Unknwon problem:
        result = -4;
     }
     
     if (result == 0) {                    // Success after removing Dest:
        if (nlhs > 0) {
           plhs[0] = mxCreateDoubleScalar(1.0);
     
           if (nlhs > 1) {
              plhs[1] = mxCreateString("");
           }
        }
        
     } else if (nlhs > 0) {                // Caller catchs outputs:
        plhs[0] = mxCreateDoubleScalar((double) result);
        
        // Reply an error message:
        if (nlhs > 1) {
           switch(result) {
              case -1:
                 plhs[1] = mxCreateString("Source is not existing.");
                 break;
              case -2:
                 plhs[1] = mxCreateString("Destination is existing already.");
                 break;
              case -3:
                 plhs[1] = mxCreateString(
                                "Destination is write protected or a folder.");
                 break;
              default:  // File or folder is accessed byother program?
                 plhs[1] = mxCreateString("Unknown problem.");
           }
        }
        
     } else {           // Output is not caught - stop with an error: ----------
        switch(result) {
           case -1:
              ErrMsg = "Source is not existing.";
              break;
           case -2:
              ErrMsg = "Destination is existing already.";
              break;
           case -3:
              ErrMsg = "Destination is write protected or a folder.";
              break;
           default:     // File or folder is accessed byother program?
              ErrMsg = "Unknown problem.";
        }
        
        // Stop with an error:
        // Rely on automatic cleanup of [Source] and [Dest] names!
        mexErrMsgIdAndTxt(ERR_ID "Failed",
                       ERR_HEAD "Cannot rename file: %s\n"
                       "    Source: %ws\n    Dest:   %ws",
                       ErrMsg, Source, Dest);
                      
     }  // end: if (nlhs > 0)
  }  // end: if (result == 0)
  
  // Release memory:
  mxFree(Source);
  mxFree(Dest);
    
  return;
}
