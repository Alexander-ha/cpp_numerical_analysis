#include <iostream>
#include <cmath>
#include <math.h>       
#include <array>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>

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



//вычитание
array< double, 2 > &substract( array< double, 2 > &ksol, array< double, 2 > &ssol){  
    static std::array<double, 2> A;
    //вычитание через algorithm
    transform(ksol.begin(), ksol.end(), ssol.begin(), A.begin(), plus<>());;        
    return A;
}

array< double, 2 > &countstep( array< double, 2 > &ksol, double &tau, array< array<double, 2>, 2> &B){
    double y0 = pow(ksol[0],6) + pow(ksol[1],6) - 64;
    double y1 = exp(ksol[1]) - exp(ksol[0]) - 1;    
    const double mult = tau;
    //умножение на матрицу
    static array< double, 2 > Mrix = {B[0][0]*y0 + B[0][1]*y1, B[1][0]*y0 + B[1][1]*y1};
    //умножение на параметр 
    for_each(Mrix.begin(), Mrix.end(), [mult](double& t) {t*=mult;});

    return Mrix;
}


int main(){
    int eps = 0.001;
    int iters = 0;
    std::array<double, 2> ksol = {-1, 1};
    double tau = 3;
    std::array<std::array<double, 2>, 2> B = { { {0,1}, {0,1} } };

    
  
    
    std::array<double, 2> ssol = substract(ksol, countstep(ksol, tau, B));
    while(euclid_norm(substract(ssol, ksol))>eps){

        ksol = ssol;
        if(euclid_norm(substract(ssol, ksol))<eps){
            break;
        }
    
        ssol = substract(ksol, countstep(ksol, tau, B));
    
        iters++;
    }
    cout<<iters<<endl;
    PrintArray(ksol);
    PrintArray(ssol);


}


