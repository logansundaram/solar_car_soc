
#include "Battery.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>
#include <map>
#include "Battery_Test.h"
#include "Pack_test.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* 
    Loading Parameters from CSV
    Can be done for all future battery types
*/

/*Data format for CSV File with parameters
    Temp_0 | SOC_0 | R0_00 | R1_00 | ...
    Temp_0 | SOC_1 | R0_01 | R1_01 | ... 
    Temp_0 | SOC_N | R0_0N | R1_0N | ...
    Temp_1 | SOC_0 | R0_10 | R1_10 | ...
    Temp_1 | SOC_1 | R0_11 | R1_11 | ... 
    Temp_1 | SOC_N | R0_1N | R1_1N | ...
    :      |       |       |       | ...
    Temp_N | SOC_0 | R0_N0 | R1_N0 | ...
    Temp_N | SOC_1 | R0_N1 | R1_N1 | ... 
    Temp_N | SOC_N | R0_NN | R1_NN | ...

*/



//format used:
/*
    Map with the key being a pair of SOC and Temperature, the value is a struct ECM object with all of the R-2RC component values
    The whole map will be passed into when creating a battery object.
*/



std::map<std::pair<double,double>,ECM> get_param_from_csv(std::string filename){
std::ifstream data(filename);
    if(!data.is_open()) throw std::runtime_error("Cannot open data");

    std::string line;
    std::map<std::pair<double,double>,ECM> Data;

    while (std::getline(data, line)) {
        std::vector<std::string> row;
        std::istringstream lineStream(line);
        std::string cell;

        std::getline(lineStream, cell, ',');
        double temperature = std::stof(cell);
        double SOC;
        // Parse each cell in the line
        while(std::getline(lineStream,cell,',')){
            float SOC = std::stof(cell);
            ECM rc;

            std::getline(lineStream, cell, ',');
            rc.R0 = std::stof(cell);
            std::getline(lineStream, cell, ',');
            rc.R1 = std::stof(cell);
            std::getline(lineStream, cell, ',');
            rc.R2 = std::stof(cell);
            std::getline(lineStream, cell, ',');
            rc.C1 = std::stof(cell);
            std::getline(lineStream, cell, ',');
            rc.C2 = std::stof(cell);

            Data[{temperature,SOC}] = rc;
        }

        


    }
    
    data.close();
    return Data;
}

void print_ECM(std::map<std::pair<double,double>,ECM> batt){
        for (const auto& [key, value] : batt) {
        std::cout << "Temp: " << key.first << ", SOC: " << key.second
                  << ", R0: " << value.R0 << ", R1: " << value.R1 << ", R2: " << value.R2 
                  << ", C1: " << value.C1 << ", C2: " << value.C2 << std::endl;
    }
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Battery Tests!:
Table of contents:
1. CC Charge
2. CC Discharge
3. Drive Cycle
4. Constant power draw
*/



using dataset = std::map<std::pair<double,double>,ECM>;
int main(int argc, char** argv)
{
    //create a dataset of R-2RC values.
    //The format is pair <double Temperature ,double SOC> as the key (unique values that organize the map), with the value being a struct ECM with the format
    /*
    struct ECM {
        double R0;
        double R1;
        double R2;
        double C1;
        double C2;
    };
    */
    dataset batt = get_param_from_csv("Enpower_params_median.csv");

    //create a Battery object
    //use the dataset to fill the R-2RC values
    //First value: Capacity, second: heat capacity, third: Thermal Resistance, fourth: SOC, fifth: battery temperature, sixth: ambient temperature.
    Battery *Enpower = new Battery(13.2,143.4,4.95,0.75,19.8,19.8,batt);

    //Creates a test object
    //The only input is the battery
    Battery_Test *test1 = new Battery_Test(*Enpower); 

    //print_ECM(batt);

    //Runs a test
    //test1->C_rate_charge_SOC(13.2/3,19.8,0.0,1.0,0.5);
    std::cout << "Hellooooo"<< std::endl;

    //std::cout << Enpower->get_SOC();
    //test1->C_rate_discharge_SOC(-13.2/5,19.8,1.0,0.0,0.5);
    //test1->new_OCV(-13.2/10);
    //test1->constant_power_draw(20,-800/70,0.5);

    Battery *Enpower1 = new Battery(13.2,143.4,4.95,1.0,40,40,batt);
    
    Battery_Test *Hot_lap = new Battery_Test(*Enpower1);

    std::vector<std::pair<double,double>> HL = {
        {180,-6000/70},
        {120,-2000/70},
        {180,-6000/70},
        {120,-2000/70},
        {180,-6000/70},
        {120,-2000/70},
        {180,-6000/70},
        {120,-2000/70},
        {180,-6000/70},
        {120,-2000/70},
        {180,-6000/70},
        {120,-2000/70},

    };


    Hot_lap->drive_cycle(40,HL,0.5);

    Battery *Enpower2 =  new Battery(13.2,143.4,4.95,1.0,40,40,batt);

    Battery_Test *Hill_Climb = new Battery_Test(*Enpower2);
    std::vector<std::pair<double,double>> HC = {
        {480,-7000/70}
    };

    //Hill_Climb->drive_cycle(40,HC,0.5);

    
    return 0;
}
