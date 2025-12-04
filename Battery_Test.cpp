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
#include "Battery_Test.h"

Battery_Test::Battery_Test() : batt(){}

Battery_Test::Battery_Test(Battery battery) : batt(battery){};

void Battery_Test:: C_rate_discharge_SOC( double current,double starting_temp, double SOC_start, double SOC_end,double ts){
    batt.set_Temp(starting_temp);
    batt.set_SOC(SOC_start);
    

    double V1 = 0.0;
    double V2 = 0.0;
    double p_out = 0.0;

    

    double R0 = batt.get_R0();

    
    double R1 = batt.get_R1();
    double R2 = batt.get_R2();
    double C1 = batt.get_C1();
    double C2 = batt.get_C2();
    double OCV = 4.2;
    double batt_voltage = OCV;
    
    

    double p_diss = 0;
    double p_diss_total = 0.0;
    double e_diss = 0;
    double e_diss_total = 0.0;

    

    

    double total_time = 0.0;
    double delta_T = 0.0;
    double s_k = 0.0;
    double hyst_voltage = 0.0;
    double voltage = 0.0;
    

    while(batt.get_SOC()>SOC_end){
        batt.set_SOC(batt.calc_SOC(batt.get_SOC(),current,ts));
        //std::cout << "R0, " << batt.get_R0() << ", R1, " << batt.get_R1() << ", R2, " << batt.get_R2() << ", C1, " << batt.get_C1() << ", C2, " <<batt.get_C2() << std::endl;
        R0 = batt.get_R0();
        R1 = batt.get_R1();
        R2 = batt.get_R2();
        C1 = batt.get_C1();
        C2 = batt.get_C2();

        V1 = batt.Calc_V1(V1,R1,current,C1,ts);
        V2 = batt.Calc_V2(V2,R2,current,C2,ts);

        total_time+=ts;

        p_diss = batt.calc_P_diss(V1,V2,current,R0,R1,R2);
        e_diss = p_diss*ts/3600; // Wh
        p_diss_total+=p_diss;
        e_diss_total += e_diss;
        delta_T += batt.calculate_delta_T(p_diss,delta_T,ts);
  
        //std::cout<<  << "," << batt.get_hysteresis_voltage(current,current,ts,hyst_voltage,s_k,1) <<  "," << batt.get_hysteresis_voltage(current,current,ts,hyst_voltage,s_k,10) <<  "," << batt.get_hysteresis_voltage(current,current,ts,hyst_voltage,s_k,100)<< std::endl;

        std::cout << batt.get_SOC() << ", " << batt.get_Temp()+delta_T << std::endl;

        //std::cout<<" SOC" << batt.get_SOC() <<  " power dissipated " << p_diss << " Energy Dissipated: " << e_diss << " Battery_temperature " << batt.get_Temp() + delta_T << std::endl;

    }

    //std::cout << "Total Power Dissipated: " << p_diss_total << " Watts, Total energy dissipated: " << e_diss_total << " Wh, Final Temperature: " << batt.get_Temp() << " Battery Time: " << total_time;
}

void Battery_Test:: C_rate_charge_SOC(double current,double starting_temp, double SOC_start, double SOC_end,double ts){
 batt.set_Temp(starting_temp);
    batt.set_SOC(SOC_start);

    double V1 = 0.0;
    double V2 = 0.0;
    double p_out = 0.0;

    double R0 = batt.get_R0();
    double R1 = batt.get_R1();
    double R2 = batt.get_R2();
    double C1 = batt.get_C1();
    double C2 = batt.get_C2();
    double OCV = 4.2;
    double batt_voltage = OCV;

    double p_diss = 0;
    double p_diss_total = 0.0;
    double e_diss = 0;
    double e_diss_total = 0.0;

    double total_time = 0.0;
    double delta_T = 0.0;
    double hyst_voltage =0.0;
    double s_k = 0.0;
    while(batt.get_SOC()<SOC_end){
        batt.calc_SOC(batt.get_SOC(),current,ts);
        //std::cout << "R0, " << batt.get_R0() << ", R1, " << batt.get_R1() << ", R2, " << batt.get_R2() << ", C1, " << batt.get_C1() << ", C2, " <<batt.get_C2() << std::endl;
        R0 = batt.get_R0();
        R1 = batt.get_R1();
        R2 = batt.get_R2();
        C1 = batt.get_C1();
        C2 = batt.get_C2();

        V1 = batt.Calc_V1(V1,R1,current,C1,ts);
        V2 = batt.Calc_V2(V2,R2,current,C2,ts);

        total_time+=ts;

        p_diss = batt.calc_P_diss(V1,V2,current,R0,R1,R2);
        e_diss = p_diss*ts/3600; // Wh
        p_diss_total+=p_diss;
        e_diss_total += e_diss;
        
        delta_T += batt.calculate_delta_T(p_diss,delta_T,ts);
        batt.set_Temp(batt.get_amb()+delta_T);
        s_k = batt.s_k(current,s_k);
        hyst_voltage = batt.get_hysteresis_voltage(current,current,ts,hyst_voltage,s_k,1);
        //std::cout<<"get hyst voltage " << hyst_voltage << std::endl;
        //std::cout << batt.get_SOC() << ", " << batt.return_OCV(-1)+hyst_voltage - V1 -V2 - current*R0 << std::endl;
        

        std::cout<<"SOC" << batt.get_SOC() <<  " power dissipated " << p_diss << " Energy Dissipated: " << e_diss << " Battery_temperature " << batt.get_Temp() << std::endl;

    }

    std::cout << "Total Power Dissipated: " << p_diss_total << " Watts, Total energy dissipated: " << e_diss_total << " Wh, Final Temperature: " << batt.get_Temp() << " Battery Time: " << total_time;


}

void Battery_Test:: drive_cycle(double starting_temp,std::vector<std::pair<double,double>> drive_cycle,double ts){
    double V1 = 0;
    double V2 = 0;
    double p_out = 0;
    double delta_t = 0;
    double SOC = 1;
    double R0 = batt.get_R0();
    double R1 = batt.get_R1();
    double R2 = batt.get_R2();
    double C1 = batt.get_C1();
    double C2 = batt.get_C2();
    double OCV = batt.return_OCV(-1.0);
    //Calculate  load current
    double I_L = 0;

    double p_diss = 0.0;
    double e_diss = 0.0;
    double e_del = 0.0;
    int num_cells = 70;
    const double max_timestep = 0.05;
    double hyst = 0.0;
    double delta_T = 0.0;
    double time = 0.0;

    for(int i = 0; i < drive_cycle.size(); i++){
        double time_remaining = drive_cycle[i].first;
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
            SOC =batt. calc_SOC(SOC,I_L, ts); // Simple Coulomb Counter
            //Use new SOC to find new paramters
            R0 = batt.get_R0();
            R1 = batt.get_R1();
            R2 = batt.get_R2();
            C1 = batt.get_C1();
            C2 = batt.get_C2();


            //change units on p_out
            p_out = ((float)drive_cycle[i].second);
            //Calculate V1, V2, and load current
            V1 = batt.Calc_V1(V1,R1,I_L,C1,ts);
            V2 = batt.Calc_V2(V2,R2,I_L,C2,ts);
            I_L = batt.Calc_IL(R0,OCV,V1,V2,hyst, p_out);
            //std::cout << "V1: " << V1 << "V2: " << V2 << std::endl;

            //std::cout << "Load Current:  " << I_L << " SOC: " << SOC << std::endl;
            //std::cout << "OCV: " << OCV << std::endl;
            //calculate total load current
            p_diss = batt.calc_P_diss(V1,V2,I_L,R0,R1,R2); // Power Lost to Heat
            delta_T += batt.calculate_delta_T(p_diss,delta_T,ts);

            batt.set_Temp(batt.get_amb() + delta_T);
            e_diss += p_diss*ts/3600; // Wh
            e_del += p_out*ts/3600; // Wh
            //std::cout << i << std::endl;
            time_remaining -= ts;
            time += ts;
            if(SOC <= 0.0001) break;
            std::cout << SOC << ", " << time << ", " << batt.get_Temp() << "," << batt.calc_P_diss(V1,V2,I_L,R0,R1,R2) << std::endl;
        
        }   
        if(SOC <= 0.0001) break;
        

    }

    std::cout<<"Total Energy Dissipated: " << e_diss << std::endl;
    std::cout << "SoC is " << SOC << std::endl;
    std::cout<< "Total Energy Delivered to Load:" << e_del << std::endl;

}

void Battery_Test:: constant_power_draw( double starting_temp, double P,double ts){
    //First time step
    //initial conditions
    double V1 = 0;
    double V2 = 0;
    double p_out = 0;
    double delta_t = 0;
    batt.set_SOC(1);
    double R0 = batt.get_R0();
    double R1 = batt.get_R1();
    double R2 = batt.get_R2();
    double C1 = batt.get_C1();
    double C2 = batt.get_C2();
    double OCV = batt.return_OCV(p_out);
    //Calculate load current
    double I_L = 0;

    double p_diss = 0.0;
    double e_diss = 0.0;
    double e_del = 0.0;
    
    double time = 0.0;

    double delta_T = 0;
    double hyst=0.0;

    std::cout<<" Time, Temperature" << std::endl;

    while(batt.get_SOC()>0.0){
        
        //Calculate SOC based on previous I_L
        batt.calc_SOC(batt.get_SOC(),I_L,ts); // Simple Coulomb Counter
        //Use new SOC to find new paramters
        R0 = batt.get_R0();
        R1 = batt.get_R1();
        R2 = batt.get_R2();
        C1 = batt.get_C1();
        C2 = batt.get_C2();
        OCV = batt.return_OCV(p_out);
        //change units on p_out
        //Calculate V1, V2, and load current
        V1 = batt.Calc_V1(V1,R1,I_L,C1,ts);
        V2 = batt.Calc_V2(V2,R2,I_L,C2,ts);
        I_L = batt.Calc_IL(R0,OCV,V1,V2,hyst,P);
        //std::cout << "V1: " << V1 << "V2: " << V2 << std::endl;

        //std::cout << "Load Current:  " << I_L << " SOC: " << SOC << std::endl;
        //std::cout << "OCV: " << OCV << std::endl;
        //calculate total load current
        p_diss = batt.calc_P_diss(V1,V2,I_L,R0,R1,R2); // Power Lost to Heat
        delta_T += batt.calculate_delta_T(p_diss,delta_T,ts);
        e_diss += p_diss*ts/3600; // Wh
        e_del += p_out*ts/3600; // Wh
        time += ts;
        //std::cout << i << std::endl;
        
        std::cout<< time  << " , " << batt.get_Temp() +delta_T<< std::endl;
    }
    
    while(time < 15506 + 1200){
        delta_T += batt.calculate_delta_T(0,delta_T,ts);
        time += ts;
        std::cout<< time  << " , " << batt.get_Temp() +delta_T<< std::endl;
    }

}

void Battery_Test:: new_OCV(double current){
    batt.set_SOC(1.0);
    double SOC = batt.get_SOC();
    while(batt.get_SOC() > 0.0){
        if(batt.get_SOC()>0.5){
            std::cout <<batt.get_SOC() <<  ", " << batt.return_OCV(-1) - current*(batt.get_R0()+batt.get_R1()+batt.get_R2()) << std::endl;
        }
        else{
            std::cout << batt.get_SOC() <<  ", " << batt.return_OCV(1) +  current*(batt.get_R0()+batt.get_R1()+batt.get_R2()) << std::endl;
        }
        SOC = batt.calc_SOC(SOC, current,0.5);
        batt.set_SOC(SOC);
    }
}