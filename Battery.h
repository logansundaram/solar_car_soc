#ifndef BATTERY_H
#define BATTERY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>

struct ECM {
    double R0;
    double R1;
    double R2;
    double C1;
    double C2;
};

class Battery{
        double Q_n;
        double C_n;
        double R_th;
        double SOC;
        double Temp;
        double amb;
        std::map<std::pair<double,double>,ECM> Params;
 
        using Coordinate = std::pair<double,double>;
        
        double interpolate(double SOC, double T, int flag);

        double map_params(double SOC, double T, int flag);

    public:
        Battery();

        Battery(double Qn, double Cn, double Rth, double SOC, double battery_temp, double amb, std::map<std::pair<double,double>,ECM> Equiv);

        double get_R0();

        double get_R1();

        double get_R2();

        double get_C1();

        double get_C2();

        double return_OCV(double power);

        double calc_SOC(double prev_SOC, double I_L, double ts);
        
        double Calc_IL(double R0, double OCV, double V1, double V2, double hyst, double P);

        double Calc_V2(double prev_V2, double R2, double I_L, double C_2, double delta_t);

        double Calc_V1(double prev_V1, double R1, double I_L, double C_1, double delta_t);

        double calc_P_diss(double V1, double V2, double IL, double R0, double R1, double R2);

        double calculate_delta_T(double P_diss, double T_prev, double ts);

        double get_hysteresis_voltage(double I_L, double prev_I_L, double ts,double prev_hyst, double prev_s_K, double gamma);

        double s_k(double I_L, double prev_s_K);

        double sgn(double I_L);

        double get_voltage(double V1, double V2, double I_L, double R0, double hyst);

        double get_Qn();

        double get_Cn();

        double get_Rth();

        double get_SOC();

        double get_Temp();

        double get_amb();

        void set_Qn(double new_Qn);

        void set_Cn(double new_Cn);

        void set_Rth(double new_Rth);

        void set_SOC(double new_SOC);

        void set_Temp(double new_Temp);



        
        

};

#endif// BATTERY_H