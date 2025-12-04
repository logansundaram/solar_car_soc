/*
Old simulator

Still works, but annoying to work with

You have to change all values to match your battery -> can be prone to mistakes 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>


std::vector<std::pair<double,double> > power; //First column is delta t (s), second column is Power (Wh)

double Q_n = 13.2*3600; //
double T_amb = 20;
double Remaining_Charge = Q_n; // Amp-Seconds
double P_diss = 0.0;
double E_diss = 0.0;
double R_th = 4;
double C_n = 213;
/*
PARAMETER FLAGS FOR READING IN CSV
C1 = 0
C2 = 1
R0 = 2
R1 = 3
R2 = 4
OCV = 5
power = 6
power_charge = 7
*/

#define C1_FLAG 0
#define C2_FLAG 1
#define R0_FLAG 2
#define R1_FLAG 3
#define R2_FLAG 4
#define OCV_DISCHARGE_FLAG     5
#define POWER_FLAG   6
#define OCV_CHARGE_FLAG 7

double Calc_SOC(double Current, double delta_t)
{
    Remaining_Charge = Remaining_Charge + Current*delta_t;
    double new_SOC = Remaining_Charge/Q_n;
    return new_SOC;
}

double Calc_V1(double prev_V1, double R1, double I_L, double C_1, double delta_t)
{
    double new_V1 = prev_V1 + delta_t * (-1 * prev_V1/R1 + I_L)/C_1;
    return new_V1;
}

double Calc_V2(double prev_V2, double R2, double I_L, double C_2, double delta_t)
{
    double new_V2 = prev_V2 + delta_t * (-1 * prev_V2/R2 + I_L)/C_2;
    return new_V2;
}

double Calc_IL(double R0, double OCV, double V1, double V2, double P){
    
    if(P == 0) return 0;
    //std::cout << "R0; " << R0 << " v1: " << V1 << " V2 " << V2 << "Power_out: " << P<<  "OCV: " << OCV<< std::endl;
    double I_L1 = (-1*(OCV+V1+V2) + sqrt(abs(pow(OCV+V1+V2,2) - 4 * (R0 * -P))))/(2*R0);
    double I_L2 = (-1*(OCV+V1+V2) - sqrt(abs(pow(OCV+V1+V2,2) - 4 * (R0 * -P))))/(2*R0);
    if(abs(I_L1) > abs(I_L2)){
        return I_L2;
    }
    else{
        return I_L1;
    }
}
double return_OCV(double SOC, double Power){

        //error checking
    if(SOC >= 1){
        SOC = 1;
    }
    else if(SOC < 0){
        SOC = 0;
    }
 

// ocv polyfit charging -4783.78848787996x^10 + 27347.9312861984x^9 +	-67026.6988407755x^8 +
//	92144.6069787111x^7	- 78029.8368377740x^6 + 42070.3825796194x^5 - 14424.1632098688x^4
//	3053.81797893469x^3	- 375.676774229697x^2 + 24.9063595277509x - 2.72423333860835

// ocv polyfit discharging 4313.04847052565x^10 + 24891.2165949856x^9 + -61632.9921177867x^8 + 85715.6445387034x^7
//-73578.8063784470x^6 + 40329.2056130290x^5 -14113.3560968725x^4 + 3067.12971117682x^3	-390.407463955115x^2
//+ 26.9957718142629x + 2.60049866650654
    
    if(Power > 0){
        double OCV = -2147.665528*pow(SOC,10)+	12225.18877*pow(SOC,9)	-29649.31849*pow(SOC,8)	+40108.79954*pow(SOC,7)	-33272.1398*pow(SOC,6)	+17542.01074*pow(SOC,5)	-5913.575841*pow(SOC,4)	+1260.295618*pow(SOC,3)	-165.6978598*pow(SOC,2)	+13.49868912*SOC	+2.810697493;

        return OCV;
    }  
    else{
        double OCV = 656.066758235956*pow(SOC,10)	-2212.74515839224*pow(SOC,9)+	2369.13403011035*pow(SOC,8)+	110.185175339433*pow(SOC,7)	-2361.54165234520*pow(SOC,6)	+2256.65970770073*pow(SOC,5)	-1066.60421354360*pow(SOC,4)	+292.595299268461*pow(SOC,3)	-48.4676488348531*pow(SOC,2)+	6.17904316492938*SOC	+2.74072164077091;


        return OCV;
       
    }  
       
}


double return_Param(int parameter_flag, double SOC)
{
    //error checking
    if(SOC > 1){
        SOC = 1;
    }
    else if(SOC < 0){
        SOC = 0;
    }
    /*
    0.764683051
0.64283833
0.50437338
0.39533691
0.265840507

*/
    switch(parameter_flag){
        case C1_FLAG:
            if(0.784193001000000 < SOC && SOC <= 1.0)
            {
                 return 237232.7729;

            }    
            else if( 0.670603757000000 < SOC && SOC <= 0.784193001000000)
            {
                return  -1732218.587527529802173*(SOC- 0.784193001000000) + 237232.772900000010850;
            }
            else if(0.555858516000000 < SOC && SOC <= 0.670603757000000)
            {
               return  -776948.675370336626656*(SOC- 0.670603757000000) + 433994.172699999995530;
            }
            else if(0.451429289000000 < SOC && SOC <= 0.555858516000000)
            {
               return  2672682.849601098801941*(SOC- 0.555858516000000) + 523145.335699999995995;
            }
            else if(0.350681661000000 < SOC && SOC <=0.451429289000000)
            {
                return  262914.232581237447448*(SOC- 0.451429289000000) + 244039.131699999998091;
            }
            else if (0.238675630000000 < SOC && SOC <= 0.350681661000000)
            {
                return 288390.409977119998075*(SOC- 0.350681661000000) + 217551.146399999997811;
            }
            else if (0.130011937000000 < SOC && SOC <= 0.238675630000000){
                return 480553.046361124317627*(SOC- 0.238675630000000) + 185249.681199999991804;
            }

            else if(0.05560093 < SOC && SOC <= 0.130011937000000){
                return 853722.704357434762642*(SOC- 0.130011937000000) + 133031.012500000011642;
            }
            else if(0.0 < SOC && SOC <= 0.05560093){
                return 853722.704357434762642*(SOC- 0.130011937000000) + 133031.012500000011642;
            }
            
        break;
        case C2_FLAG:
          if(0.784193001000000 < SOC && SOC <= 1.0)
            {
                 return -15981.188553380983649*(SOC- 0.784193001000000) + 6182.167394000000058;
            }    
            else if( 0.670603757000000 < SOC && SOC <= 0.784193001000000)
            {
                return  -15981.188553380983649*(SOC- 0.784193001000000) + 6182.167394000000058;
            }
            else if(0.555858516000000 < SOC && SOC <= 0.670603757000000)
            {
               return  -3356.885450264559040*(SOC- 0.670603757000000) + 7997.458520000000135;
            }
            else if(0.451429289000000 < SOC && SOC <= 0.555858516000000)
            {
               return  4634.261977252790530*(SOC- 0.555858516000000) + 8382.645150000000285;
            }
            else if(0.350681661000000 < SOC && SOC <=0.451429289000000)
            {
                return  3275.680316761402082*(SOC- 0.451429289000000) + 7898.692753999999695;
            }
            else if (0.238675630000000 < SOC && SOC <= 0.350681661000000)
            {
                return 14154.627093250008329*(SOC- 0.350681661000000) + 7568.675731999999698;
            }
            else if (0.130011937000000 < SOC && SOC <= 0.238675630000000){
                return 22405.192597310306155*(SOC- 0.238675630000000) + 5983.272130999999717;
            }

            else if(0.05560093 < SOC && SOC <= 0.130011937000000){
                return 26031.942908661349065*(SOC- 0.130011937000000) + 3548.641161000000011;
            }
            else if(0.0 < SOC && SOC <= 0.05560093){
                return 1611.578075;
;
            }
            
        break;
        case R0_FLAG:
           if(0.784193001000000 < SOC && SOC <= 1.0)
            {
                 return  -0.000885708861659*(SOC - 0.784193001000000) + 0.002300153000000;
            }    
            else if( 0.670603757000000 < SOC && SOC <= 0.784193001000000)
            {
                return   -0.000885708861659*(SOC - 0.784193001000000) + 0.002300153000000;
            }
            else if(0.555858516000000 < SOC && SOC <= 0.670603757000000)
            {
               return  0.000731856060157*(SOC- 0.670603757000000) + 0.002400760000000;
            }
            else if(0.451429289000000 < SOC && SOC <= 0.555858516000000)
            {
               return  0.000797678986937*(SOC- 0.555858516000000) + 0.002316783000000;
            }
            else if(0.350681661000000 < SOC && SOC <=0.451429289000000)
            {
                return  -0.000330131841913*(SOC- 0.451429289000000) + 0.002233482000000;
            }
            else if (0.238675630000000 < SOC && SOC <= 0.350681661000000)
            {
                return -0.001637965369918*(SOC- 0.350681661000000) + 0.002266742000000;
            }
            else if (0.130011937000000 < SOC && SOC <= 0.238675630000000){
                return -0.004600377423212*(SOC- 0.238675630000000) + 0.002450204000000;
            }
            else if(0.05560093 < SOC && SOC <= 0.130011937000000){
                return -0.008522811685642*(SOC- 0.130011937000000) + 0.002950098000000;
            }
            else if(0.0 < SOC && SOC <= 0.05560093){
                return -0.008522811685642*(SOC- 0.130011937000000) + 0.002950098000000;
            }
            
        break;
        case R1_FLAG:
            if(0.784193001000000 < SOC && SOC <= 1.0)
            {
                 return 0.018893760750798*(SOC- 0.784193001000000) + 0.003390619000000;
            }    
            else if( 0.670603757000000 < SOC && SOC <= 0.784193001000000)
            {
                return 0.018893760750798*(SOC- 0.784193001000000) + 0.003390619000000;

            }
            else if(0.555858516000000 < SOC && SOC <= 0.670603757000000)
            {
               return  -0.001720358929744*(SOC- 0.670603757000000) + 0.001244491000000;
            }
            else if(0.451429289000000 < SOC && SOC <= 0.555858516000000)
            {
               return  -0.013927221734582*(SOC- 0.555858516000000) + 0.001441894000000;
            }
            else if(0.350681661000000 < SOC && SOC <=0.451429289000000)
            {
                return  -0.007443321643265*(SOC- 0.451429289000000) + 0.002896303000000;
            }
            else if (0.238675630000000 < SOC && SOC <= 0.350681661000000)
            {
                return -0.005789036485009*(SOC- 0.350681661000000) + 0.003646200000000;
            }
            else if (0.130011937000000 < SOC && SOC <= 0.238675630000000){
                return -0.014902742169825*(SOC- 0.238675630000000) + 0.004294607000000;
            }
            else if(0.05560093 < SOC && SOC <= 0.130011937000000){
                return -0.063817615047193*(SOC- 0.130011937000000) + 0.005913994000000;
            }
            else if(0.0 < SOC && SOC <= 0.05560093){
                return -0.063817615047193*(SOC- 0.130011937000000) + 0.005913994000000;
            }
            
        break;
        case R2_FLAG:
            if(0.784193001000000 < SOC && SOC <= 1.0)
            {
                 return 0.000293680975639*(SOC- 0.784193001000000) + 0.004072542000000;
            }    
            else if( 0.670603757000000 < SOC && SOC <= 0.784193001000000)
            {
                return  0.000293680975639*(SOC- 0.784193001000000) + 0.004072542000000;
            }
            else if(0.555858516000000 < SOC && SOC <= 0.670603757000000)
            {
               return  -0.001122861382983*(SOC- 0.670603757000000) + 0.004039183000000;
            }
            else if(0.451429289000000 < SOC && SOC <= 0.555858516000000)
            {
               return  -0.009131763466946*(SOC- 0.555858516000000) + 0.004168026000000;
            }
            else if(0.350681661000000 < SOC && SOC <=0.451429289000000)
            {
                return  -0.006874196581581*(SOC- 0.451429289000000) + 0.005121649000000;
            }
            else if (0.238675630000000 < SOC && SOC <= 0.350681661000000)
            {
                return -0.008946134338070*(SOC- 0.350681661000000) + 0.005814208000000;
            }
            else if (0.130011937000000 < SOC && SOC <= 0.238675630000000){
                return-0.034495928644722*(SOC- 0.238675630000000) + 0.006816229000000;
            }
            else if(0.05560093 < SOC && SOC <= 0.130011937000000){
                return -0.126711912392208*(SOC- 0.130011937000000) + 0.010564684000000;
            }
            else if(0.0 < SOC && SOC <= 0.05560093){
                return -0.126711912392208*(SOC- 0.130011937000000) + 0.010564684000000;
            }
            
        break;

    }
    return -1;
}

double calc_P_diss(double V1, double V2, double IL, double R0, double R1, double R2){

    return pow(V1,2)/R1 + pow(V2,2)/R2 + pow(IL,2)*R0;
}

double calculate_delta_T(double P_diss, double T_prev, double ts){
    return ts*(-(T_prev)/(R_th) + P_diss)/C_n;
}

void read_power(std::string parameter_csv)
{
    std::ifstream parameter_file(parameter_csv);
    if(!parameter_file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }
    
    std::string line;
    std::string first_val, second_val;

    if(parameter_file.good())
    {
        while(std::getline(parameter_file,line))
        {
            std::stringstream ss(line);
            std::getline(ss,first_val,',');
            std::getline(ss,second_val);

            double first = atof(first_val.c_str());
            double second = atof(second_val.c_str());

            std::pair<double, double> csv_pair;
            csv_pair  = std::make_pair(first,second);

            power.push_back(csv_pair);
            
            
        }
    }

    parameter_file.close();
    if(parameter_file.is_open())
    {
        throw std::runtime_error("File not closed properly");
    }
}

void run_sim(){

    //First time step
    //initial conditions
    double V1 = 0;
    double V2 = 0;
    double p_out = 0;
    double delta_t = 0;
    double SOC = 1;
    double R0 = return_Param(R0_FLAG,SOC);
    double R1 = return_Param(R1_FLAG,SOC);
    double R2 = return_Param(R2_FLAG,SOC);
    double C1 = return_Param(C1_FLAG,SOC);
    double C2 = return_Param(C2_FLAG,SOC);
    double OCV = return_OCV(SOC,p_out);
    //Calculate load current
    double I_L = 0;

    double p_diss = 0.0;
    double e_diss = 0.0;
    double e_del = 0.0;
    int num_cells = 80;
    const double max_timestep = 0.05;
    for(int i = 0; i < power.size(); i++){
        double time_remaining = power[i].first;
        while(time_remaining > 0) {
        double ts;
        if(time_remaining < max_timestep){
            ts = time_remaining;
        }
        else
        {
            ts = max_timestep;
        }
        //Calculate SOC based on previous I_L
        SOC = Calc_SOC(I_L, ts); // Simple Coulomb Counter
        //Use new SOC to find new paramters
        R0 = return_Param(R0_FLAG,SOC);
        R1 = return_Param(R1_FLAG,SOC);
        R2 = return_Param(R2_FLAG,SOC);
        C1 = return_Param(C1_FLAG,SOC);
        C2 = return_Param(C2_FLAG,SOC);
        OCV = return_OCV(SOC,p_out);
        //change units on p_out
        p_out = ((float)power[i].second / (float)num_cells);
        //Calculate V1, V2, and load current
        V1 = Calc_V1(V1,R1,I_L,C1,ts);
        V2 = Calc_V2(V2,R2,I_L,C2,ts);
        I_L = Calc_IL(R0,OCV,V1,V2,p_out);
        //std::cout << "V1: " << V1 << "V2: " << V2 << std::endl;

        //std::cout << "Load Current:  " << I_L << " SOC: " << SOC << std::endl;
        //std::cout << "OCV: " << OCV << std::endl;
        //calculate total load current
        p_diss = calc_P_diss(V1,V2,I_L,R0,R1,R2); // Power Lost to Heat
        e_diss += p_diss*ts/3600; // Wh
        e_del += p_out*ts/3600; // Wh
        //std::cout << i << std::endl;
        time_remaining -= ts;
        if(SOC <= 0.0001) break;
        }
        if(SOC <= 0.0001) break;
        std:: cout << "power: " << p_diss << " energy: " << e_diss << " i: " << i << " SOC: " << SOC << std::endl ;
        if(e_diss == 0.0495186)
        {

        }
    }
    std::cout<<"Total Power Dissipated: " << p_diss << std::endl;
    std::cout<<"Total Energy Dissipated: " << e_diss << std::endl;
    std::cout << "SoC is " << SOC << std::endl;
    std::cout<< "Total Energy Delivered to Load:" << e_del << std::endl;
}

/// @brief Models a C rate discharge from the battery 
/// @param C_rate 
void C_rate_discharge(double C_rate){
    
    //initial conditions
    
    double V1 = 0.0;
    double V2 = 0.0;
    double p_out = 0.0;
    double delta_t = 0.0;
    double SOC = 1.0;
    double R0 = return_Param(R0_FLAG,SOC);
    double R1 = return_Param(R1_FLAG,SOC);
    double R2 = return_Param(R2_FLAG,SOC);
    double C1 = return_Param(C1_FLAG,SOC);
    double C2 = return_Param(C2_FLAG,SOC);
    double OCV = 4.2;
    double batt_voltage = OCV;

    double p_diss = 0;
    double e_diss = 0;


    
    //Calculate load current
    double I_L = C_rate*-13.2;
    double ts = 0.5; //max time step
    
    double total_cap = 0.0;
    double delta_T = 0;

    double time = 0;
    while(SOC > 0){
        
        //Calculate SOC based on previous I_L
        SOC = SOC + I_L*ts/13.2/3600; // Simple Coulomb Counter
        //Use new SOC to find new paramters
        R0 = return_Param(R0_FLAG,SOC);
        R1 = return_Param(R1_FLAG,SOC);
        R2 = return_Param(R2_FLAG,SOC);
        C1 = return_Param(C1_FLAG,SOC);
        C2 = return_Param(C2_FLAG,SOC);
        OCV = return_OCV(SOC,-1);
        //change units on p_out
        //Calculate V1, V2, and load current
        V1 = Calc_V1(V1,R1,I_L,C1,ts);
        V2 = Calc_V2(V2,R2,I_L,C2,ts);
        time += ts;
       
        
        batt_voltage = OCV + V1 + V2 +I_L * R0;

        p_diss = calc_P_diss(V1,V2,I_L,R0,R1,R2); // Calculate power dissipated the impedance of the battery
        delta_T += calculate_delta_T(p_diss,delta_T,ts); // calculate the change in temperature from the thermal circuit

        std:: cout << "power, " << p_diss << ", time, " <<  time << ", SOC, " << SOC << ", Temperature, " << T_amb+delta_T << std::endl ;
    }
    std::cout<<"Total Power Dissipated: " << p_diss << std::endl;
    std::cout<<"Total Energy Dissipated: " << e_diss << std::endl;
    std::cout << "SoC is " << SOC << std::endl;
    
    
}

void C_rate_Charge(double C_rate){
        //initial conditions
    
    double V1 = 0.0;
    double V2 = 0.0;
    double p_out = 0.0;
    double delta_t = 0.0;
    double SOC = 0.00;
    double R0 = return_Param(R0_FLAG,SOC);
    double R1 = return_Param(R1_FLAG,SOC);
    double R2 = return_Param(R2_FLAG,SOC);
    double C1 = return_Param(C1_FLAG,SOC);
    double C2 = return_Param(C2_FLAG,SOC);
    double OCV = 2.75;
    double batt_voltage = OCV;

    double p_diss = 0;
    double e_diss = 0;


    
    //Calculate load current
    double I_L = C_rate*13.2;
    double ts = 0.5; //max time step
    
    double total_cap = 0.0;
    double delta_T = 0;

    double time = 0;
    while(SOC < 1){
        
        //Calculate SOC based on previous I_L
        SOC = SOC + I_L*ts/13.2/3600; // Simple Coulomb Counter
        //Use new SOC to find new paramters
        R0 = return_Param(R0_FLAG,SOC);
        R1 = return_Param(R1_FLAG,SOC);
        R2 = return_Param(R2_FLAG,SOC);
        C1 = return_Param(C1_FLAG,SOC);
        C2 = return_Param(C2_FLAG,SOC);
        OCV = return_OCV(SOC,1);
        //change units on p_out
        //Calculate V1, V2, and load current
        V1 = Calc_V1(V1,R1,I_L,C1,ts);
        V2 = Calc_V2(V2,R2,I_L,C2,ts);
        time += ts;
       
        
        batt_voltage = OCV + V1 + V2 +I_L * R0;

        p_diss = calc_P_diss(V1,V2,I_L,R0,R1,R2); // Calculate power dissipated the impedance of the battery
        delta_T += calculate_delta_T(p_diss,delta_T,ts); // calculate the change in temperature from the thermal circuit

        std::cout << "Power Dissipated: " << p_diss << "Temperature, " << T_amb + delta_T << std::endl;
    }

}

void Power_draw_test(double P){
    
    //First time step
    //initial conditions
    double V1 = 0;
    double V2 = 0;
    double p_out = 0;
    double delta_t = 0;
    double SOC = 1;
    double R0 = return_Param(R0_FLAG,SOC);
    double R1 = return_Param(R1_FLAG,SOC);
    double R2 = return_Param(R2_FLAG,SOC);
    double C1 = return_Param(C1_FLAG,SOC);
    double C2 = return_Param(C2_FLAG,SOC);
    double OCV = return_OCV(SOC,p_out);
    //Calculate load current
    double I_L = 0;

    double p_diss = 0.0;
    double e_diss = 0.0;
    double e_del = 0.0;
    double ts = 0.5;
    double time = 0.0;

    double delta_T = 0;
    while(SOC > 0){
        
        //Calculate SOC based on previous I_L
        SOC = SOC + I_L*ts/13.2/3600; // Simple Coulomb Counter
        //Use new SOC to find new paramters
        R0 = return_Param(R0_FLAG,SOC);
        R1 = return_Param(R1_FLAG,SOC);
        R2 = return_Param(R2_FLAG,SOC);
        C1 = return_Param(C1_FLAG,SOC);
        C2 = return_Param(C2_FLAG,SOC);
        OCV = return_OCV(SOC,p_out);
        //change units on p_out
        //Calculate V1, V2, and load current
        V1 = Calc_V1(V1,R1,I_L,C1,ts);
        V2 = Calc_V2(V2,R2,I_L,C2,ts);
        I_L = Calc_IL(R0,OCV,V1,V2,P);
        //std::cout << "V1: " << V1 << "V2: " << V2 << std::endl;

        //std::cout << "Load Current:  " << I_L << " SOC: " << SOC << std::endl;
        //std::cout << "OCV: " << OCV << std::endl;
        //calculate total load current
        p_diss = calc_P_diss(V1,V2,I_L,R0,R1,R2); // Power Lost to Heat
        delta_T += calculate_delta_T(p_diss,delta_T,ts);
        e_diss += p_diss*ts/3600; // Wh
        e_del += p_out*ts/3600; // Wh
        time += ts;

        //std::cout << i << std::endl;
        
        std:: cout <<"time, " << time, "Temperature: ", T_amb+delta_T;
    }
}






int main(int argc, char** argv)
{
   
    
    //C_rate_discharge(1);
    C_rate_Charge(1);
    //Power_draw_test(-50);

    return 0;
}