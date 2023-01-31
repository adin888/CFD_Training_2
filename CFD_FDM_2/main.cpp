#include"../Header/cfd_head.h"

int main(int argc, char* argv[])
{
	int choice;

	cout << "Please enter the serial number of the simulation you wish to run: " << endl;
	cin >> choice;
	switch (choice)
	{
	case(1):
		Roe_Riemann_Solver();
		break;
	case(2):
		HLLC_Riemann_Solver();
		break;
	case(3):
		Rusanov_Riemann_Solver();
		break;
	case(4):
		FFT_Solver();              //Existing problems!
		break;
	case(5):
		FST_Solver();
		break;
	case(6):
		GaussSeidelSolver();
		break;
	case(7):
		ConjugateGradientSolver();
		break;
	case(8):
		LDC_Solver();
		break;
	case(9):
		FittingFunction();
		break;
	case(10):
		HOM1DHelmholtz();
		break;

	default:
		cout << "Invalid input" << endl;
		break;
	}

	return 0;
}
