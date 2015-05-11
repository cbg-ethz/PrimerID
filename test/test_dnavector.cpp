#include <src/pidalyse/statistics.hpp>

#include <iostream>

int main()
{
	const int L = 2;
	DNAvector<int> LOLvector(L);
	
	std::string DNA1;
	std::string DNA2;
	int number1;
	int number2;
	int TEMP;
	
	for (int i = 0; i < 4*4; ++i)
	{
		TEMP = 2*i;
		LOLvector[i] = TEMP;
		
		DNA1 = number_To_DNA(i, L);
		number1 = DNA_to_number(DNA1);
		DNA2 = number_To_DNA(number1, L);
		number2 = DNA_to_number(DNA2);
		
		std::cout << "i: " << i << "\tDNA1: " << DNA1 << "\tnumber1: " << number1 << "\tDNA2: " << DNA2 << "\tnumber2: " << number2 << "\tTEMP: " << TEMP << "  \tVECTOR: " << LOLvector[DNA2] << '\n';
		
		if ((i != number1) || (number1 != number2))
			std::cout << "ERROR!\n";
		if (DNA1 != DNA2)
			std::cout << "ERROR!\n";
	}
	
	return 0;
}