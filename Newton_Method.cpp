#include <iostream>
#include <cmath>
#include <math.h>       
#include <array>
#include <iomanip>
#include <algorithm>
#include <numeric>

using namespace std;

//шаблон вывода вектора
template<typename T>
void PrintArray(const T& e) {
    std::cout << std::setw(3) << e << ' ';
}

template<typename T, std::size_t N>
void PrintArray(const std::array<T,N>& A) {

    for (const auto& e : A) 
        PrintArray(e);
    std::cout << '\n';
}

//вычисление нормы вектора
template<typename T, std::size_t N>
double euclid_norm(const std::array<T,N>& A){
    //подсчет нормы через numeric
    return sqrt(inner_product(A.begin(), A.end(), A.begin(), 0.0L));
}

//функция вывода якобиана
array< double, 2 > &jacstep( array< double, 2 > &ksol){
    //подсчет переменных якобиана
    double x00 = 6*pow(ksol[0], 5);
    double x01 = 6*pow(ksol[1], 5);
    double x10 = -exp(ksol[0]);
    double x11 = exp(ksol[1]);

    //вычисление обратного якобиана
    double jacobian[2][2] = { {x00, x01}, {x10, x11}};
    double det = jacobian[0][0]*jacobian[1][1] - jacobian[0][1]*jacobian[1][0];
    
    double resjac[2][2] = { {x11/det, -x01/det}, {-x10/det, x00/det}};
    double y0 = pow(ksol[0],6) + pow(ksol[1],6) - 64;
    double y1 = exp(ksol[1]) - exp(ksol[0]) - 1;    
    double res[2] = {y0*resjac[0][0] + y1*resjac[0][1], y0*resjac[1][0] + y1*resjac[1][1]};
    
    //перегоняем классический массив в контейнер
    static std::array<double, 2> A;
    std::move(std::begin(res), std::end(res), A.begin());    
    return A;
}

// вычитание векторов
array< double, 2 > &substract( array< double, 2 > &ksol, array< double, 2 > &ssol){  
    static std::array<double, 2> A;
    //вычитание через algorithm
    transform(ksol.begin(), ksol.end(), ssol.begin(), A.begin(), minus<>());;        
    return A;
}





int main(){
    int eps = 0.001;
    int iters = 0;
    std::array<double, 2> ksol = {90, 11.7};
//    std::array<double, 2> osol = {1, 1};
//    cout<<euclid_norm(osol)<<endl;
//    std::array<double, 2> sisol = substract(osol, ksol);
    std::array<double, 2> ssol = substract(ksol, jacstep(ksol));

    ksol = ssol; 

    ssol = substract(ksol, jacstep(ksol));

 
    while(euclid_norm(substract(ssol, ksol))>eps){
        cout<<euclid_norm(substract(ssol, ksol))<<endl;

        if(euclid_norm(substract(ssol, ksol))<eps){
            break;
        }
        ksol = ssol;
        PrintArray(ksol);
        ssol = substract(ksol, jacstep(ksol));
        PrintArray(ssol);
        iters++;

       
    }
    
    cout<<iters<<endl;
    PrintArray(ksol);

};