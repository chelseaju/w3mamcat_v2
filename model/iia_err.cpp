//#include "stdafx.h"
#include <iostream>
#include <stdlib.h>
#include "iia_err.h"


void  crash( const char*  Msg )
{   
#if 0
  MessageBox( NULL, Msg, "Crashing",
              MB_OK | MB_ICONHAND );
	throw DivideByZeroException();
#endif
    std::cerr << Msg;
	throw DivideByZeroException();
}

void  warn( const char* Msg )
{                                 
#if 0
  MessageBox( NULL, Msg, "Warning",
              MB_OK | MB_ICONEXCLAMATION );
#endif
    std::cout << Msg;
}
