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
    // ����� ����������
    //
    int calculate(double dTensor[3][3], double* dCharacteristicValues, double dCharacteristicVectors[3][3]);

private:
    //
    // ���� ������������ ��� lamda^2, lamda^1 � lamda^0
    //                       c[0]     c[1]      c[2]
    void getCoefs(double* c, double dTensor[3][3]);

    // �������� �� ��, ������� �� ������������ ������� � ��������� ������������� � ���������
    // ���������� 1 � ������ �������, � ������ ������ - 0
    // iIndex - ������ ��������
    int check(double* dArray, int* pIndex);

    // 
    // ���������� �������������� ������������ �������
    // 
    void getCharacteristicVector(double dTensorIn[3][3], double dCharacteristicNumber, double* dCharacteristicVector);

    //
    // Solve cubic equation x^3 + a*x^2 + b*x + c = 0
    // x - array of size 3
    // return 3: 3 real roots x[0], x[1], x[2]
    // return 1: 1 real root x[0] and pair of complex roots
    int SolveP3(double* x, double a, double b, double c);
};

