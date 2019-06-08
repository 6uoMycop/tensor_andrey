#include <iostream>
#include <fstream>
#include "Tensor.h"


using namespace std;

int main(int argc, char* argv[])
{
    double dTensor[3][3];
    double dCharacteristicValues[3];
    double dCharacteristicVectors[3][3] = { 0.0 };
    Tensor prog;


    if (argc == 1)
    {
        cout << "No file name. Terminate" << endl;
        return -1;
    }
    ifstream F(argv[1]);
    if (!F.is_open())
    {
        cout << "File was not opened. Terminate" << endl;
        return -1;
    }

    int k = 0, j = 0;
    for (int i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            F >> dTensor[i][j];
        }
        k++;
    }

    cout << "Input:" << endl;
    for (int m = 0; m < 3; m++)
    {
        for (int n = 0; n < 3; n++)
        {
            cout << dTensor[m][n] << "\t";
        }
        cout << endl;
    }
    cout << endl;
    F.close();


    prog.calculate(dTensor, dCharacteristicValues, dCharacteristicVectors);


    cout << "Characteristic values:" << endl;
    for (int m = 0; m < 3; m++)
    {
        cout << "lamda" << m << " = " << dCharacteristicValues[m] << endl;
    }
    cout << endl;

    for (int m = 0; m < 3; m++)
    {
        cout << "Characteristic vector for lamda" << m << ":" << endl;
        for (int n = 0; n < 3; n++)
        {
            cout << dCharacteristicVectors[m][n] << endl;
        }
    }


    system("pause");
    return 0;
}