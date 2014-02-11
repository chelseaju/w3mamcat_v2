// Programmed by Hsiao-Te Su
// 12/26/93       
//#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>   
#include <iostream>
#include "iia_mem.h"
#include "iia_err.h" 

#define CPP  1

#if CPP
#endif

// allocate ByteNum to the return pointer
// includes simple error handling                            
// usage:  (Type *) alloc( ByteNum, __FILE__, __LINE__ );
void*  alloc( int ByteNum, const char* FileName, int LineNum )
{
  void*  ret;

#if CPP
  ret = (void *) new char[ ByteNum ];  
#else
  ret = malloc( ByteNum );
#endif

  if( ret == NULL )                   
  {      
      std::cerr << "dynamic memory allocation failure." << std::endl;
      std::cerr << "attempting to allocate " << ByteNum << "bytes" << std::endl;
      std::cerr << "in file " << FileName << " at line " << LineNum << std::endl;
    crash( "dynamic memory allocation failure" );
    exit( EXIT_FAILURE );
  }  
  return( ret );
}                                

// deallocate the ptr
void  afree( void*  ptr )
{
#if CPP
  delete [] ptr;
#else
  if( ptr != NULL )
    free( ptr );
#endif
}
