#include <iostream>
#include "poly34.h"

using namespace std;


using namespace std;

// Ищем коэффициенты при lamda^2, lamda^1 и lamda^0
//                       c[0]     c[1]      c[2]
void coefs(double *c, double dTensor[3][3])
{
    c[0] = dTensor[0][0] + dTensor[1][1] + dTensor[2][2];

    c[1] =
        - dTensor[0][0] * dTensor[2][2]
        - dTensor[1][1] * dTensor[2][2]
        - dTensor[0][0] * dTensor[1][1]
        + dTensor[0][2] * dTensor[2][0]
        + dTensor[1][2] * dTensor[2][1];

    c[2] =
          dTensor[0][0] * dTensor[1][1] * dTensor[2][2]
        + dTensor[0][1] * dTensor[1][2] * dTensor[2][0]
        + dTensor[1][0] * dTensor[2][1] * dTensor[0][2]
        - dTensor[0][2] * dTensor[2][0] * dTensor[1][1]
        - dTensor[1][2] * dTensor[2][1] * dTensor[0][0];

    c[0] /= -1;
    c[1] /= -1;
    c[2] /= -1;
}

int main()
{
    double dTensor[3][3];
    double dCoefficients[3];
    double dCharacteristicValues[3];

    cout << "Print coordinates:" << endl;

    int k = 0, j = 0;
    for (int i = 0; i < 3; i++)
    {
        for (j = k; j < 3; j++)
        {
            cout << "A " << i + 1 << " " << j + 1 << endl << "> ";
            cin >> dTensor[i][j];
            dTensor[j][i] = dTensor[i][j];
        }
        k++;
    }
    
    for (int m = 0; m < 3; m++)
    {
        for (int n = 0; n < 3; n++)
        {
            cout << dTensor[m][n] << " ";
        }
        cout << endl;
    }
    
    coefs(dCoefficients, dTensor);

    //for (int m = 0; m < 3; m++)
    //{
    //    cout << dCoefficients[m] << endl;
    //}

    SolveP3(dCharacteristicValues, dCoefficients[0], dCoefficients[1], dCoefficients[2]);

    cout << "Characteristic values:" << endl;
    for (int m = 0; m < 3; m++)
    {
        cout << dCharacteristicValues[m] << endl;
    }



    return 0;
}