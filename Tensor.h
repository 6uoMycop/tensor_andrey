#pragma once
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

class Tensor
{
public:
    Tensor();
    ~Tensor();

    //
    // Общее вычисление
    //
    int calculate(double dTensor[3][3], double* dCharacteristicValues, double dCharacteristicVectors[3][3]);

private:
    //
    // Ищем коэффициенты при lamda^2, lamda^1 и lamda^0
    //                       c[0]     c[1]      c[2]
    void getCoefs(double* c, double dTensor[3][3]);

    // Проверка на то, остался ли единственный элемент с ненулевым коэффициентом в уравнении
    // Возвращает 1 в случае неудачи, в случае успеха - 0
    // iIndex - индекс элемента
    int check(double* dArray, int* pIndex);

    // 
    // Нахождение нормированного собственного вектора
    // 
    void getCharacteristicVector(double dTensorIn[3][3], double dCharacteristicNumber, double* dCharacteristicVector);

    //
    // Solve cubic equation x^3 + a*x^2 + b*x + c = 0
    // x - array of size 3
    // return 3: 3 real roots x[0], x[1], x[2]
    // return 1: 1 real root x[0] and pair of complex roots
    int SolveP3(double* x, double a, double b, double c);
};

