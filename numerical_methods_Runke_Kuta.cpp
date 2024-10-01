#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace std;


//численное решение
struct Solarr{ 
    double num_sol[2];
};

//y' = f(x,y)
Solarr func_resh(double tmp, double y1, double y2){
    Solarr num_sol;
    num_sol.num_sol[0] = -2 * y1 - 3 * y2 + 4 * exp(-tmp);
    num_sol.num_sol[1] = 3 * y1 - 2 * y2 + 4 * exp(-tmp);
    return num_sol;
}


//аналитическое решение
double exact_sol(double tmp){
    return exp(-tmp);
}

//подсчет коэффициентов и численного решения для step+1
void Runge_Kutta(int step, double tmp, double h, double exact_sol(double), Solarr num_resh, double* y1, double* y2, double y_hp1, double y_hp2, double matrix[4][4]){
    double k1[3];
    double k2[3];

    double y1_start, y2_start;
    y1_start = y_hp1;
    y2_start = y_hp2;

  //  func_resh(tmp, *y1, *y2);
    k1[0] =   -2 * y1_start - 3 * y2_start + 4 * exp(-tmp);
    k2[0] =  3 * y1_start - 2 * y2_start - 2 * exp(-tmp) ;
    for (int i = 0; i <= step;  i++){
      for (int j = 0; j <= step-1;  j++){
      if (i<3){
        if(j <=2){
        y1_start += matrix[i][j]*h*k1[j];
        y2_start += matrix[i][j]*h*k2[j];
        }
        else{
        y1_start += 0;
        y2_start += 0;
        }
       // std::cout<<k1[i]<<std::endl;
       // func_resh(tmp + h*matrix[i][0], y1_start, y2_start);
        k1[i] =  -2 * y1_start - 3 * y2_start + 4 * exp(-tmp + h*matrix[i][0]);
        k2[i] =  3 * y1_start - 2 * y2_start - 2 * exp(-tmp + h*matrix[i][0]);
       // std::cout<<k1[i]<<std::endl;

      }
    }

    }
//std::cout<<'_'<<matrix[1][1]<<'_'<<matrix[1][2]<<'_'<<matrix[1][3]<<std::endl;

for (int i = 0; i<=2; i++){
    //std::cout<<k1[i]<<'_'<<matrix[3][i+1]<<std::endl;
     *y1 += h*k1[i]*matrix[3][i+1];
     *y2 += h*k2[i]*matrix[3][i+1];
    }
    std::cout<<*y1<<std::endl;
 

}

//вычисление нормы ошибки
double max_norm(double exact_sol(double), Solarr num_resh, double tmp, double y1, double y2){
    Solarr recount = func_resh(tmp, y1, y2);
    double* y0 = recount.num_sol; 
    double norm0 = abs(sqrt(exact_sol(tmp)*sqrt(2))  - sqrt(y0[0]*y0[0]+y0[1]*y0[1]));
   // std::cout<<norm0<<std::endl;
    return norm0;
}

int main(){

    std::cout<<"reading matrix!";
    double matrix[4][4] = {{0, 0, 0, 0}, 
                           {0.666666666, 0.666666666, 0, 0 },
                            {0.6666666666, 0, 0.6666666666, 0},
                             {0, 0.25, 0.375, 0.375}};
    double eps = 10e-10;
    Solarr num_resh;
    double t_start = 0, y1 = 0, y2 = 0, t_end = 4, h = 0.001;
    double err = max_norm(exact_sol, num_resh, t_start, y1, y2);
    int iter = 0;
    double N = t_end/h;
    cout<<N;
    double t = t_start;
//    Runge_Kutta(iter, t, h, exact_sol, num_resh, &y1, &y2, y1, y2, reinterpret_cast<double (*)[4]>(matrix));
    //std::cout<<y1<<std::endl;
    //std::cout<<y2<<std::endl;
    while((iter<=N) && (err>=10e-20)){
       Runge_Kutta(iter, t, h, exact_sol, num_resh, &y1, &y2, y1, y2, matrix);
       err = max_norm(exact_sol, num_resh, t, y1, y2);

       t+=h;
       iter+=1;
    std::cout<<err<<std::endl;
    }
    std::cout<<err;
  


}