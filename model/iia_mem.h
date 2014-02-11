// Programmed by Hsiao-Te Su 
// 12/26/93
#ifndef IIA_MEM_H
#define IIA_MEM_H

void*  alloc( int  ByteNum, const char* FileName, int LineNum ); 
void  afree( void* ptr );

#endif  // IIA_MEM_H
