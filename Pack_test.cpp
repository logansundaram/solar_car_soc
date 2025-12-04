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
#include "Pack_test.h"
#include "Battery_Test.h"

Pack::Pack():pack(){}


Pack::Pack(std::vector<Battery*> t_pack) : pack(t_pack){}

void Pack::pack_charge(double I_L){




}

std::vector<double> slice_power(double p_draw, Battery *batt,double ts, double V1_prev, double V2_prev,double I_L_prev, double time){
 //First time step
    //initial conditions

    double p_out = 0;
    double delta_t = 0;
    double R0 = batt->get_R0();
    double R1 = batt->get_R1();
    double R2 = batt->get_R2();
    double C1 = batt->get_C1();
    double C2 = batt->get_C2();
    //std::cout << "Params: " << R0 << " " << R1 << " " << R2 << " " << C1 << " " << C2 << std::endl;
    double OCV = batt->return_OCV(p_draw);
    //std::cout<< OCV << std::endl;
    double V1 = batt->Calc_V1(V1_prev,R1,I_L_prev,C1,ts);
    double V2 = batt->Calc_V2(V2_prev,R2,I_L_prev,C2,ts);

    //std::cout << V1 << " " << V2 <<std::endl;

    //Calculate load current
    double I_L = batt->Calc_IL(R0,OCV,V1,V2,0,p_draw);
    //std::cout<< "I_L " << I_L << std::endl;

    
    double p_diss = 0.0;
    double e_diss = 0.0;
    double e_del = 0.0;
    

    double delta_T = 0;
    double hyst=0.0;

    p_diss = batt->calc_P_diss(V1,V2,I_L,R0,R1,R2); // Power Lost to Heat
    delta_T = batt->calculate_delta_T(p_diss,delta_T,ts);

    time += ts;
    //std::cout<<batt->get_SOC() << std::endl;
    batt->set_SOC(batt->get_SOC()+(I_L*ts)/13.2/3600);
    batt->set_Temp(batt->get_amb()+delta_T);

    double V_t = batt->get_voltage(V1,V2,I_L,R0,0);
        
    std::cout<< "SOC " <<batt->get_SOC() << ", time, " << time  << " , " << " Voltage " << V_t << ", Temperature, " <<  batt->get_Temp()<< ", " <<std::endl;
    std::vector<double> stuff_to_return = {V_t,V1,V2,I_L,time};
    
    //std::cout<< V1 << " " << V2 << " " << I_L << " " << time << std::endl;
    return stuff_to_return;
}

void Pack::pack_power_draw(double p_draw){
    double p_drawn = p_draw/pack.size();

    double lowest_SOC = 100;
    int index = 100;
    for(int i = 0; i<pack.size(); i++){
        if(pack[i]->get_SOC() < lowest_SOC){
            lowest_SOC = pack[i]->get_SOC();
            index = i;
        }
    }
    std::vector<std::vector<double>> stuff_to_return= {
        {0,0,0,0,0},

    };
    double smallest_V = 100;
    std::vector<double> Bat_V;
    for(int i = 0; i < pack.size(); i++){
        Bat_V.push_back(0);
    }
    while( smallest_V > 2.75){
        lowest_SOC = pack[index]->get_SOC();
        for(int i = 0; i < pack.size(); i++){
            //std::cout << "Battery " << i << std::endl;
            stuff_to_return[i]=  slice_power(p_drawn,pack[i],0.5,stuff_to_return[i][1],stuff_to_return[i][2],stuff_to_return[i][3],stuff_to_return[i][4]);
            //std::cout << stuff_to_return[i][0] << std::endl;
            Bat_V[i] = stuff_to_return[i][0];
            //std::cout<< " end of first for" << std::endl;
        }
        smallest_V = 100;
        int index = 0;
        for(int i = 0; i < pack.size(); i++){
            if(Bat_V[i]< smallest_V){
                smallest_V = Bat_V[i];
                
                index = i;
            }
        }
        
    }

    for(int i = 0; i < pack.size(); i++){
        std::cout << " Battery No. " << i;
    }
    std::cout << std::endl;
    for(int i = 0; i < pack.size(); i++){
        std::cout << " SOC: " << pack[i]->get_SOC();
    }
    std::cout << std::endl;
    for(int i = 0; i < pack.size(); i++){
        std::cout << " Voltage: " << Bat_V[i];
    }
    std::cout << std::endl;
    for(int i = 0; i < pack.size(); i++){
        std::cout << " Temp: " << pack[i]->get_Temp();
    }
}

void Pack::pack_drive_cycle(std::vector<std::pair<double,double>> dri_cy){
    std::vector<std::vector<double>> stuff_to_return= {
        {0,0,0,0,0},
        // {0,0,0,0,0},
        // {0,0,0,0,0},
        // {0,0,0,0,0},
        // {0,0,0,0,0},
        // {0,0,0,0,0},
        // {0,0,0,0,0},
        // {0,0,0,0,0},

    };

    std::vector<double> Bat_V;
    for(int i = 0; i < pack.size(); i++){
        Bat_V.push_back(0);
    }
    for(int i = 0; i < dri_cy.size(); i++){
        double time = dri_cy[i].second;
        std::cout<<time<<std::endl;
        double power = dri_cy[i].first;
        
        const double max_ts = 0.5;
        double ts = 0.5;
        double t_elapsed = 0;
        

        while(t_elapsed < time){
            
            if(time-t_elapsed < max_ts){
                
                ts = time-t_elapsed;
            }
            else{
                
                ts = max_ts;
            }
            for(int i = 0; i < pack.size(); i++){
                
                //std::cout << "Battery " << i << std::endl;
                stuff_to_return[i]=  slice_power(power,pack[i],0.5,stuff_to_return[i][1],stuff_to_return[i][2],stuff_to_return[i][3],stuff_to_return[i][4]);
                //std::cout << stuff_to_return[i][0] << std::endl;
                Bat_V[i] = stuff_to_return[i][0];
                //std::cout<< " end of first for" << std::endl;
                
            }
            t_elapsed+=ts;
            



        }
        for(int i = 0; i < pack.size(); i++){
            std::cout << " Battery No. " << i;
        }
        std::cout << std::endl;
        for(int i = 0; i < pack.size(); i++){
            std::cout << " SOC: " << pack[i]->get_SOC();
        }
        std::cout << std::endl;
        for(int i = 0; i < pack.size(); i++){
            std::cout << " Voltage: " << Bat_V[i];
        }
        std::cout << std::endl;
        for(int i = 0; i < pack.size(); i++){
            std::cout << " Temp: " << pack[i]->get_Temp();
        }
    }

}

