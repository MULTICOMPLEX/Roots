
#include <MULTICOMPLEX.hpp>

int main(int argc, char** argv)
{
	MX0 x;

	x.random(-1.5, 1.5);

	std::cout.setf(std::ios::fixed, std::ios::floatfield);
	std::cout.precision(5);
	std::cout << "multicomplex roots..." << std::endl;
	std::cout << "ini" << x << std::endl << std::endl;

	auto fx1 = [](const auto& x) {
		return ((pow(8, x) - pow(2, x)) / (pow(6, x) - pow(3, x))) - 2; };

	//auto fx1 = [](const auto& x) { return pow(x, 3.365) - 10.2 / x; };

	//auto fx1 = [](const auto& x) { return exp(x)+2; };

	root(fx1, x, 30);
	
	std::cout << std::endl;

	return 0;
}
