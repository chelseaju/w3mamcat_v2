// Programmed by Hsiao-Te Su
// 12/26/93

#ifndef IIA_ERR_H
#define IIA_ERR_H
void  crash( const char*  Msg );
void  warn( const char* Msg );

class DivideByZeroException {
public:
	DivideByZeroException()
		: message ("attempted to divide by zero") {}
	const char *what() const { return message; }
private:
	const char * message;
};
#endif  // IIA_ERR_H
