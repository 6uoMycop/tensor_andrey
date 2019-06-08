#include "Tensor.h"

#define	TwoPi  6.28318530717958648
const double eps = 1e-14;



Tensor::Tensor()
{
}


Tensor::~Tensor()
{
}

//
// Общее вычисление
//
int Tensor::calculate(double dTensor[3][3], double* dCharacteristicValues, double dCharacteristicVectors[3][3])
{
    double dCoefficients[3];

    getCoefs(dCoefficients, dTensor);

    SolveP3(dCharacteristicValues, dCoefficients[0], dCoefficients[1], dCoefficients[2]);

    for (int m = 0; m < 3; m++)
    {
        getCharacteristicVector(dTensor, dCharacteristicValues[m], dCharacteristicVectors[m]);
    }

    return 0;
}

//
// Ищем коэффициенты при lamda^2, lamda^1 и lamda^0
//                       c[0]     c[1]      c[2]
//
void Tensor::getCoefs(double* c, double dTensor[3][3])
{
    c[0] = dTensor[0][0] + dTensor[1][1] + dTensor[2][2];

    c[1] =
         - dTensor[0][0] * dTensor[1][1]
         - dTensor[0][0] * dTensor[2][2]
         - dTensor[1][1] * dTensor[2][2]
         + dTensor[0][2] * dTensor[2][0]
         + dTensor[1][2] * dTensor[2][1]
         + dTensor[0][1] * dTensor[1][0];

    c[2] =
           dTensor[0][1] * dTensor[1][2] * dTensor[2][0]
         + dTensor[1][0] * dTensor[2][1] * dTensor[0][2]
         - dTensor[0][2] * dTensor[2][0] * dTensor[1][1]
         - dTensor[1][2] * dTensor[2][1] * dTensor[0][0]
         - dTensor[0][1] * dTensor[1][0] * dTensor[2][2]
         + dTensor[0][0] * dTensor[1][1] * dTensor[2][2];

    c[0] *= -1;
    c[1] *= -1;
    c[2] *= -1;
}

//
// Проверка на то, остался ли единственный элемент с ненулевым коэффициентом в уравнении
// Возвращает 1 в случае неудачи, в случае успеха - 0
// iIndex - индекс элемента
//
int Tensor::check(double* dArray, int* pIndex)
{
    int flag = 0;
    for (int i = 0; i < 3; i++)
    {
        if (dArray[i] != 0.0)
        {
            flag++;
            if (pIndex != NULL)
            {
                *pIndex = i;
            }
        }
    }
    if (flag == 1)
    {
        return 0;
    }
    return 1;
}

// 
// Нахождение нормированного собственного вектора
// 
void Tensor::getCharacteristicVector(double dTensorIn[3][3], double dCharacteristicNumber, double* dCharacteristicVector)
{
    double dTensor[3][3];

    dCharacteristicVector[0] = 1.0; // проверить, норм ли единицы
    dCharacteristicVector[1] = 1.0;
    dCharacteristicVector[2] = 1.0;

    for (int i = 0; i < 3; i++)
    {
        memcpy(dTensor[i], dTensorIn[i], 3 * sizeof(double));
    }

    dTensor[0][0] -= dCharacteristicNumber;
    dTensor[1][1] -= dCharacteristicNumber;
    dTensor[2][2] -= dCharacteristicNumber;

    double dBuffer[3][6] = { 0 };
    bool bReady[3] = { false };

    int flag0 = 0;
    int iIndex = 0;
    int iIndex1 = 0;
    int iLine = 0;

    const double dZeroArray[6] = { 0 };

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            dBuffer[i][j] = dTensor[i][j];
        }
    }

    while (1)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (abs(dBuffer[i][j]) < eps)
                {
                    dBuffer[i][j] = 0.0;
                }
            }
        }

        for (int i = 0; i < 3; i++)
        {
            flag0 = 0;
            iIndex = 0;
            if (check(dBuffer[i], &iIndex) == 0 && memcmp(&dBuffer[i][3], dZeroArray, 3 * sizeof(double)) == 0) // Exists! j: in a_j*x_j a_j != 0  =>  a_j = 0
            {
                if (!bReady[iIndex])
                {
                    dCharacteristicVector[iIndex] = 0.0;
                    for (int k = 0; k < 3; k++)
                    {
                        dBuffer[k][iIndex] = 0.0;
                        dBuffer[k][iIndex + 3] = 0.0;
                    }
                    bReady[iIndex] = 1;
                }
            }
        }


        flag0 = 0;
        iIndex = 0;
        for (int i = 0; i < 3; i++)
        {
            if (check(dBuffer[i], NULL) != 0) // если элемент не единственный (есть, что перенести)
            {
                for (int j = 0; j < 3; j++)
                {

                    if (dBuffer[i][j] != 0.0) // переносим вправо
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (k != j)
                            {
                                dBuffer[i][k + 3] -= dBuffer[i][k] / dBuffer[i][j];
                                dBuffer[i][k] = 0.0;
                            }
                        }
                        dBuffer[i][j] = 1;
                        flag0 = 1;
                        iLine = i; // номер уравнения, из которого выражали
                        break;
                    }

                }
            }

            if (flag0) // что-то выразили, можем подставить
            {
                break;
            }
        }

        // если сократились одинаковые переменные, то записываем 1 в результат 
        for (int i = 0; i < 3; i++)
        {
            if (check(dBuffer[i], &iIndex) == 0 && check(&dBuffer[i][3], &iIndex1) == 0)
            {
                if (iIndex == iIndex1 && dBuffer[i][iIndex] == dBuffer[i][iIndex1 + 3])
                {
                    dCharacteristicVector[iIndex] = 1.0;
                    bReady[iIndex] = true;
                    dBuffer[i][iIndex] = 0.0;
                    dBuffer[i][iIndex + 3] = 0.0;
                    continue;
                }
            }
        }

        //for (int i = 0; i < 3; i++) // подставляем то, что выразили
        //{
        if (check(dBuffer[iLine], &iIndex) == 0)
        {
            for (int j = 0; j < 3; j++)
            {

                //if (j != i)
                if (j != iLine)
                {

                    if (dBuffer[j][iIndex] != 0.0) // слева
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            dBuffer[j][k] += dBuffer[j][iIndex] * dBuffer[iLine][k + 3];
                        }
                        dBuffer[j][iIndex] = 0.0;
                    }
                    if (dBuffer[j][iIndex + 3] != 0.0) // справа
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            dBuffer[j][k + 3] += dBuffer[j][iIndex + 3] * dBuffer[iLine][k + 3];
                        }
                        dBuffer[j][iIndex + 3] = 0.0;
                    }
                }

            }
            //break;
        }
        //}

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                if (abs(dBuffer[i][j]) < eps)
                {
                    dBuffer[i][j] = 0.0;
                }
            }
        }
        // проверка
        flag0 = 0;
        for (int i = 0; i < 3; i++)
        {
            if ((memcmp(dZeroArray, dBuffer[i], 6 * sizeof(double)) == 0) ||
                ((dBuffer[i][0] == 1 || dBuffer[i][1] == 1 || dBuffer[i][2] == 1) &&
                (check(dBuffer[i], NULL) == 0 && check(&dBuffer[i][3], NULL) == 0)))
            {
                flag0++;
            }
        }
        if (flag0 == 3)
        {
            break;
        }

    }

    // запись результата
    for (int i = 0; i < 3; i++)
    {
        check(dBuffer[i], &iIndex1);
        if (!bReady[iIndex1])
        {
            if (check(&dBuffer[i][3], &iIndex) == 0)
            {
                dCharacteristicVector[iIndex1] = dBuffer[i][iIndex + 3];
                bReady[iIndex1] = true;
            }
            else
            {
                dCharacteristicVector[iIndex1] = 0;
                bReady[iIndex1] = true;//?
            }
        }
    }

    // нормирование
    double dAbsX = sqrt(
        dCharacteristicVector[0] * dCharacteristicVector[0] +
        dCharacteristicVector[1] * dCharacteristicVector[1] +
        dCharacteristicVector[2] * dCharacteristicVector[2]
    );

    if (dAbsX != 0.0)
    {
        dCharacteristicVector[0] /= dAbsX;
        dCharacteristicVector[1] /= dAbsX;
        dCharacteristicVector[2] /= dAbsX;
    }

    return;
}




// solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
// Thanks to Alexandr Rakhmanin <rakhmanin (at) gmail.com>
// public domain
//

//=============================================================================
// _root3, root3 from http://prografix.narod.ru
//=============================================================================
static double _root3(double x)
{
    double s = 1.;
    while (x < 1.)
    {
        x *= 8.;
        s *= 0.5;
    }
    while (x > 8.)
    {
        x *= 0.125;
        s *= 2.;
    }
    double r = 1.5;
    r -= 1. / 3. * (r - x / (r * r));
    r -= 1. / 3. * (r - x / (r * r));
    r -= 1. / 3. * (r - x / (r * r));
    r -= 1. / 3. * (r - x / (r * r));
    r -= 1. / 3. * (r - x / (r * r));
    r -= 1. / 3. * (r - x / (r * r));
    return r * s;
}

double root3(double x)
{
    if (x > 0) return _root3(x); else
        if (x < 0) return-_root3(-x); else
            return 0.;
}


//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
int Tensor::SolveP3(double* x, double a, double b, double c) {	// solve cubic equation x^3 + a*x^2 + b*x + c = 0
    double a2 = a * a;
    double q = (a2 - 3 * b) / 9;
    double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
    // equation x^3 + q*x + r = 0
    double r2 = r * r;
    double q3 = q * q * q;
    double A, B;
    if (r2 <= (q3 + eps)) {//<<-- FIXED!
        double t = r / sqrt(q3);
        if (t < -1) t = -1;
        if (t > 1) t = 1;
        t = acos(t);
        a /= 3; q = -2 * sqrt(q);
        x[0] = q * cos(t / 3) - a;
        x[1] = q * cos((t + TwoPi) / 3) - a;
        x[2] = q * cos((t - TwoPi) / 3) - a;
        return(3);
    }
    else {
        //A =-pow(fabs(r)+sqrt(r2-q3),1./3); 
        A = -root3(fabs(r) + sqrt(r2 - q3));
        if (r < 0) A = -A;
        B = (A == 0 ? 0 : B = q / A);

        a /= 3;
        x[0] = (A + B) - a;
        x[1] = -0.5 * (A + B) - a;
        x[2] = 0.5 * sqrt(3.) * (A - B);
        if (fabs(x[2]) < eps) { x[2] = x[1]; return(2); }
        return(1);
    }
}// SolveP3(double *x,double a,double b,double c)
