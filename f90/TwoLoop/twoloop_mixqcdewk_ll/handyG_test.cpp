#include <iostream>
#include "handyG_cpp.h"

int main()
{
	std::cout << handyG::G(1.0,2.0,1.0) << std::endl;
	std::cout << handyG::G(1.0,0.0,0.5,0.3) << std::endl;
	std::cout << handyG::G(1.0,0.0,0.5,std::complex<double>(1.0,1.0),0.3) << std::endl;
	std::cout << handyG::G(std::complex<double>(1.0,0.0000000000001),0.0,5.0,1.0/0.3) << std::endl;
	std::cout << handyG::G(std::complex<double>(1.0,-0.0000000000001),0.0,5.0,1.0/0.3) << std::endl;
	return 0;
}
