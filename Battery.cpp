#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>
#include "Battery.h"

/*
Implementation for Battery.h
Scope: Creates a Battery object with the R-2RC model and has getters and setters for all parameter values needed to model the battery.
*/

Battery::Battery() : Q_n(0), C_n(0.0), R_th(0.0),SOC(1.0),Temp(20), Params(){};

Battery::Battery(double Qn, double Cn, double Rth, double SOC, double battery_temp, double amb, std::map<std::pair<double,double>,ECM> Equiv) : Q_n(Qn), C_n(Cn), R_th(Rth), SOC(SOC), Temp(battery_temp), amb(amb), Params(Equiv){};

double Battery::get_R0(){
    int flag = 0;
  
    double R0 = map_params(SOC,Temp,flag);
    if(R0 < 0) throw std::runtime_error("Negative R0");
    return R0;
}

double Battery::get_R1(){
    int flag = 1;
    double R1 = map_params(SOC,Temp,flag);
    if(R1 < 0) throw std::runtime_error("Negative R1");
    return R1;
}

double Battery::get_R2(){
    int flag = 2;
    double R2 = map_params(SOC,Temp,flag);
    if(R2 < 0) R2 = 0.004072542;
    return R2;
}

double Battery::get_C1(){
    int flag = 3;
    double C1 = map_params(SOC,Temp,flag);
    if(C1 < 0) throw std::runtime_error("Negative C1");
    return C1;
}

double Battery::get_C2(){
    int flag = 4;
    double C2 = map_params(SOC,Temp,flag);
    if(C2 < 0) throw std::runtime_error("Negative C2");
    return C2;
}
//double OCV = OCV_charge[0]*pow(SOC,10)+	OCV_charge[1]*pow(SOC,9) + OCV_charge[2]*pow(SOC,8)	+ OCV_charge[3]*pow(SOC,7) + OCV_charge[4]*pow(SOC,6) + OCV_charge[5]*pow(SOC,5) + OCV_charge[6]*pow(SOC,4)	+ OCV_charge[7]*pow(SOC,3)	+ OCV_charge[8]*pow(SOC,2)	+ OCV_charge[9]*SOC	+ OCV_charge[0];
//double OCV = OCV_discharge[0]*pow(SOC,0) + OCV_discharge[1]*pow(SOC,9) + OCV_discharge[2]*pow(SOC,8) + OCV_discharge[3]*pow(SOC,7) + OCV_discharge[4]*pow(SOC,6) + OCV_discharge[5]*pow(SOC,5) + OCV_discharge[6]*pow(SOC,4) + OCV_discharge[7]*pow(SOC,3) + OCV_discharge[8]*pow(SOC,2) + OCV_discharge[9]*SOC	+ OCV_discharge[10];
//-2147.665528*pow(SOC,10)+	12225.18877*pow(SOC,9)	-29649.31849*pow(SOC,8)	+40108.79954*pow(SOC,7)	-33272.1398*pow(SOC,6)	+17542.01074*pow(SOC,5)	-5913.575841*pow(SOC,4)	+1260.295618*pow(SOC,3)	-165.6978598*pow(SOC,2)	+13.49868912*SOC	+2.810697493;
//656.066758235956*pow(SOC,10)	-2212.74515839224*pow(SOC,9)+	2369.13403011035*pow(SOC,8)+	110.185175339433*pow(SOC,7)	-2361.54165234520*pow(SOC,6)	+2256.65970770073*pow(SOC,5)	-1066.60421354360*pow(SOC,4)	+292.595299268461*pow(SOC,3)	-48.4676488348531*pow(SOC,2)+	6.17904316492938*SOC	+2.74072164077091;
double Battery::return_OCV(double power){
    //std::cout << SOC << std::endl;
    if(SOC >= 1){
        SOC = 1;
    }
    else if(SOC < 0){
        SOC = 0;
    }
    if(power > 0){
       double OCV = -4736.75561818155*pow(SOC,10)	+24517.2664738478*pow(SOC,9)	-54173.2191016299*pow(SOC,8)	+66835.4026038897*pow(SOC,7)	-50609.3374108806*pow(SOC,6)	+24397.7218928215*pow(SOC,5)	-7542.51311268701*pow(SOC,4)+	1477.00354052990*pow(SOC,3)	-177.372653840129*pow(SOC,2)	+13.1032330345485*pow(SOC,1)+	2.89419627230423;
        return OCV;
    }  
    else{
        double OCV = 656.066758235956*pow(SOC,10)	-2212.74515839224*pow(SOC,9)+	2369.13403011035*pow(SOC,8)+	110.185175339433*pow(SOC,7)	-2361.54165234520*pow(SOC,6)	+2256.65970770073*pow(SOC,5)	-1066.60421354360*pow(SOC,4)	+292.595299268461*pow(SOC,3)	-48.4676488348531*pow(SOC,2)+	6.17904316492938*SOC	+2.74072164077091;
        return OCV;
    
    }  

}
double Battery::Calc_V1(double prev_V1, double R1, double I_L, double C_1, double delta_t)
{
    double new_V1 = prev_V1 + delta_t * (-1 * prev_V1/R1 + I_L)/C_1;
    return new_V1;
}

double Battery::Calc_V2(double prev_V2, double R2, double I_L, double C_2, double delta_t)
{
    double new_V2 = prev_V2 + delta_t * (-1 * prev_V2/R2 + I_L)/C_2;
    return new_V2;
}

double Battery::Calc_IL(double R0, double OCV, double V1, double V2, double hyst, double P){
    
    if(P == 0) return 0;
    double I_L1 = (-1*(OCV+V1+V2) + sqrt(abs(pow(OCV-hyst+V1+V2,2) - 4 * (R0 * -P))))/(2*R0);
    double I_L2 = (-1*(OCV+V1+V2) - sqrt(abs(pow(OCV-hyst+V1+V2,2) - 4 * (R0 * -P))))/(2*R0);
    if(abs(I_L1) > abs(I_L2)){
        return I_L2;
    }
    else{
        return I_L1;
    }
}

double Battery::calc_P_diss(double V1, double V2, double IL, double R0, double R1, double R2){

    return pow(V1,2)/R1 + pow(V2,2)/R2 + pow(IL,2)*R0;
}

double Battery:: calculate_delta_T(double P_diss, double T_prev, double ts){
    return ts*(-(T_prev)/(R_th) + P_diss)/C_n;
}


double Battery::calc_SOC(double prev_SOC, double I_L, double delta_T){
    SOC = prev_SOC + I_L*delta_T/Q_n/3600;
    return SOC;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Hysteresis AHHHHHHH (Dont show our fearless leader Matthew Anderson)*/

double Battery::get_hysteresis_voltage(double I_L,double prev_I_L,double ts, double prev_hyst, double prev_s_K, double gamma){
    double M= (return_OCV(1)-return_OCV(-1))/2;
    double M_0 = 0.004;
    double eta;

    if(I_L < 0){
        eta = 1.0;
    }
    else{
        eta = 0.999;
    }
    
    
    double hyst_state = exp(-abs(eta*I_L*gamma*ts/(3600*Q_n))) * prev_hyst - (1-exp(-eta*I_L*gamma*ts/(Q_n*3600)))*sgn(I_L);

    return M_0*s_k(I_L,prev_s_K) + M*hyst_state;



}

double Battery:: s_k(double I_L, double prev_s_K){
    if(abs(I_L)>0){
        if(I_L > 0){
            return 1;
        }
        if(I_L < 0){
            return -1;
        }
        else{
            return 0;
        }
    }
    else{
        return prev_s_K;
    }
}

double Battery:: sgn(double I_L){
    if(abs(I_L)>0){
        if(I_L > 0){
            return 1;
        }
        if(I_L < 0){
            return -1;
        }
        else{
            return 0;
        }
    }
    else{
        return 0;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*  Battery Voltage!*/

double Battery:: get_voltage(double V1, double V2, double I_L, double R0, double hyst){
    double OCV = return_OCV(I_L);
    //std::cout << OCV << std::endl;
    

    return OCV + hyst + V1 +V2 +I_L*R0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Getters and Setters for Private Variables
*/

double Battery::get_Qn(){
    return Q_n;

}

double Battery::get_Cn(){
    return C_n;
}

double Battery::get_Rth(){
    return R_th;
}

double Battery::get_amb(){
    return amb;
}
double Battery::get_SOC(){
    return SOC;
}

double Battery::get_Temp(){
    return Temp;
}

void Battery::set_Qn(double new_Qn){
    Q_n = new_Qn;
}

void Battery::set_Cn(double new_Cn){
    C_n = new_Cn;
}

void Battery::set_Rth(double new_Rth){
    R_th = new_Rth;
}

void Battery::set_SOC(double new_SOC){
    SOC = new_SOC;
}

void Battery::set_Temp(double new_Temp){
    Temp = new_Temp;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Private Functions!!!*/
/// @brief Private function for get R0, R1,R2,C1,C2
/// @param SOC 
/// @param Temp 
/// @return a vector


double Battery:: interpolate(double SOC, double Temp, int flag){
    double SOC1 = -INFINITY, SOC2 = INFINITY, temp1 = -INFINITY, temp2 = INFINITY;
    bool S1_found = false, S2_found = false, T1_found = false, T2_found = false;

    for(const auto& point : Params){
        double point_SOC = point.first.second;
        double point_T = point.first.first; 

    

        if(point_SOC <= SOC && (!S1_found || point_SOC > SOC1)){
            SOC1 = point_SOC;
            S1_found = true;
        }
        if(point_SOC >= SOC && (!S2_found || point_SOC < SOC2)){
            SOC2 = point_SOC;
            S2_found = true;
        }
        if(point_T <= Temp && (!T1_found || point_T > temp1)){
            temp1 = point_T;
            T1_found = true;
        }
        if(point_T >= Temp && (!T2_found || point_T < temp2)){
            temp2 = point_T;
            T2_found = true;
        }

    }

    if(SOC1 == -INFINITY || SOC2 == INFINITY || temp1 == -INFINITY|| temp2 == INFINITY){
        for(const auto&point:Params){
            double point_SOC = point.first.second;
            double point_temp = point.first.first;

            if(SOC1 == -INFINITY && point_SOC > SOC2) SOC1 = point_SOC;
            if(SOC2 == INFINITY && point_SOC < SOC1) SOC2 = point_SOC;
            if(temp1 == -INFINITY && point_temp > temp2) temp1= point_temp;
            if(temp2 == INFINITY && point_temp < temp1) temp2= point_temp;
        }

    }
        double val_11 = 0.0;
        double val_21 = 0.0;
        double val_12 = 0.0;
        double val_22 = 0.0;

        double denom = 0.0;
        double val;


    switch(flag){

        case 0:
            
            val_11 = Params.at({temp1,SOC1}).R0;
            val_21 = Params.at({temp1,SOC2}).R0;
            val_12 = Params.at({temp2,SOC1}).R0;
            val_22 = Params.at({temp2,SOC2}).R0;

            denom = (SOC2-SOC1) * (temp2-temp1);
            val = 0.0;

            val+= ((SOC2-SOC)*(temp2-Temp) * val_11)/denom;
            val+= ((SOC-SOC1)*(temp2-Temp) * val_21)/denom;
            val+= ((SOC2-SOC)*(Temp-temp1) * val_12)/denom;
            val+= ((SOC-SOC1)*(Temp-temp1) * val_22)/denom;

            return val;
        break;

        case 1:
            val_11 = Params.at({temp1,SOC1}).R1;
            val_21 = Params.at({temp1,SOC2}).R1;
            val_12 = Params.at({temp2,SOC1}).R1;
            val_22 = Params.at({temp2,SOC2}).R1;

            denom = (SOC2-SOC1) * (temp2-temp1);
            val = 0.0;

            val+= ((SOC2-SOC)*(temp2-Temp) * val_11)/denom;
            val+= ((SOC-SOC1)*(temp2-Temp) * val_21)/denom;
            val+= ((SOC2-SOC)*(Temp-temp1) * val_12)/denom;
            val+= ((SOC-SOC1)*(Temp-temp1) * val_22)/denom;

            return val;
        break;

        case 2:
            val_11 = Params.at({temp1,SOC1}).R2;
            val_21 = Params.at({temp1,SOC2}).R2;
            val_12 = Params.at({temp2,SOC1}).R2;
            val_22 = Params.at({temp2,SOC2}).R2;

            // std::cout << std::endl;
            // std::cout << "SOC1 " << SOC1 << " SOC2 " << SOC2 << std::endl;
            // std::cout << "temp1 " << temp1 << " temp2 " << temp2 <<std::endl;
            // std::cout << val_11 <<" " << val_21 << std::endl;
            

            denom = abs((SOC2-SOC1) * (temp2-temp1));

            // std::cout << denom << std::endl;
            val = 0.0;

            val+= ((SOC2-SOC)*(temp2-Temp) * val_11)/denom;
            val+= ((SOC-SOC1)*(temp2-Temp) * val_21)/denom;
            val+= ((SOC2-SOC)*(Temp-temp1) * val_12)/denom;
            val+= ((SOC-SOC1)*(Temp-temp1) * val_22)/denom;

            return val;
        break;

        case 3:
            val_11 = Params.at({temp1,SOC1}).C1;
            val_21 = Params.at({temp1,SOC2}).C1;
            val_12 = Params.at({temp2,SOC1}).C1;
            val_22 = Params.at({temp2,SOC2}).C1;

            denom = (SOC2-SOC1) * (temp2-temp1);
            val = 0.0;

            val+= ((SOC2-SOC)*(temp2-Temp) * val_11)/denom;
            val+= ((SOC-SOC1)*(temp2-Temp) * val_21)/denom;
            val+= ((SOC2-SOC)*(Temp-temp1) * val_12)/denom;
            val+= ((SOC-SOC1)*(Temp-temp1) * val_22)/denom;

            return val;
        break;

        case 4:
            val_11 = Params.at({temp1,SOC1}).C2;
            val_21 = Params.at({temp1,SOC2}).C2;
            val_12 = Params.at({temp2,SOC1}).C2;
            val_22 = Params.at({temp2,SOC2}).C2;

            denom = (SOC2-SOC1) * (temp2-temp1);
            val = 0.0;

            val+= ((SOC2-SOC)*(temp2-Temp) * val_11)/denom;
            val+= ((SOC-SOC1)*(temp2-Temp) * val_21)/denom;
            val+= ((SOC2-SOC)*(Temp-temp1) * val_12)/denom;
            val+= ((SOC-SOC1)*(Temp-temp1) * val_22)/denom;

            return val;
        break;
    }
    return -INFINITY;


}
double evaluatePoly(const std::vector<double>& coeffs, const std::vector<std::vector<int>>& terms, double x, double y) {
    double result = 0.0;
    
    for (size_t i = 0; i < coeffs.size(); ++i) {
        double termValue = coeffs[i] * pow(x, terms[i][0]) * pow(y, terms[i][1]);
        result += termValue;
    }
    
    return result;
}

//A 2D polyfit. Coefficients are hard coded in, but this can be updated to run based off a csv.
/*

R0 Coeffs: -2.11657326333673	0.000846885579792063	3.91423897182570	-0.000137801794197026	0.00547177834585406	-2.25943315447019	1.33715944364830e-05	-0.000620873275285851	-0.0274226630476610	-4.84508488806610
R1 Coeffs: 7.22590582443639	-0.142200032478815	-2.47169988508999	-0.00236670751609961	0.169550403838385	-4.50174576034747	1.02899329653893e-05	-0.000453436776804045	-0.00209840570620878	-4.18127431467941
R2 Coeffs: -8.84564638807578	0.0411661045930874	12.3581650530000	0.00171763468028220	-0.109585279607469	-6.07027689223863	5.37782199099250e-05	-0.00434350973009193	0.0889079827226382	-4.08912110337676
C1 Coeffs: -25.9042613406855	0.267657607597391	16.0077264556218	0.00986034348323903	-0.617962609040007	6.30046535769818	-4.67989997084294e-05	0.000868275268393775	0.0667360905281632	10.9325434146271
C2 Coeffs: 21.0571460671353	-0.204124090987535	-21.1778768327544	-0.00775628669504216	0.491615610018904	4.30495900657677	5.81665580739370e-05	-0.00245415262941661	0.0167375616597262	7.96804533313503


*/
double Battery::map_params(double SOC, double T, int flag){
    std::vector<double> R0_coeffs = {
        -2.11657326333673, 0.000846885579792063, 3.91423897182570, 
        -0.000137801794197026, 0.00547177834585406, -2.25943315447019, 
        0.0000133715944364830, -0.000620873275285851, -0.0274226630476610, 
        -4.84508488806610
    }; 
    std::vector<double> R1_coeffs = {
        7.22590582443639,	-0.142200032478815,	-2.47169988508999,	
        -0.00236670751609961,	0.169550403838385,	-4.50174576034747,	1.02899329653893e-05,	
        -0.000453436776804045,	-0.00209840570620878,	-4.18127431467941,
    };
    std::vector<double> R2_coeffs = {
        -8.84564638807578,	0.0411661045930874,	12.3581650530000,	
        0.00171763468028220,	-0.109585279607469,	-6.07027689223863,	5.37782199099250e-05,	
        -0.00434350973009193,	0.0889079827226382,	-4.08912110337676,
    };
    std::vector<double> C1_coeffs = {
        -25.9042613406855,	0.267657607597391,	16.0077264556218,	0.00986034348323903,	
        -0.617962609040007,	6.30046535769818,	-4.67989997084294e-05,	0.000868275268393775,	
        0.0667360905281632,	10.9325434146271,
    };
    std::vector<double> C2_coeffs = {
        21.0571460671353,	-0.204124090987535,	-21.1778768327544,	-0.00775628669504216,
        	0.491615610018904,	4.30495900657677,	5.81665580739370e-05,	-0.00245415262941661,
            	0.0167375616597262,	7.96804533313503,
    };

    std::vector<std::vector<int>> terms = {
        {3, 0},  // Constant
        {2, 1},  // x
        {2, 0},  // y
        {1, 2},  // x^2
        {1, 1},  // xy
        {1, 0},  // y^2
        {0, 3},  // x^3
        {0, 2},  // x^2y
        {0, 1},  // xy^2
        {0, 0}   // y^3
    };

    double val = -INFINITY;
    switch(flag){
        case 0:
            val = evaluatePoly(R0_coeffs,terms,SOC,T);
        break;
        case 1:
            val = evaluatePoly(R1_coeffs,terms,SOC,T);
        break;
        case 2:
            val = evaluatePoly(R2_coeffs,terms,SOC,T);
        break;
        case 3:
            val = evaluatePoly(C1_coeffs,terms,SOC,T); 
        break;
        case 4:
            val = evaluatePoly(C2_coeffs,terms,SOC,T);
        break;           
    }
    
    return exp(val);
    
}


    

