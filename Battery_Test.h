#ifndef BATTERY_TEST_H
#define BATTERY_TEST_H
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

class Battery_Test{
    Battery batt;

    public:
        Battery_Test();

        Battery_Test(Battery battery);

        void C_rate_discharge_SOC(double current, double starting_temp, double SOC_start, double SOC_end,double ts);

        void C_rate_charge_SOC(double current, double starting_temp, double SOC_start, double SOC_end,double ts);

        void drive_cycle( double starting_temp, std::vector<std::pair<double,double>> drive_cycle,double ts);

        void constant_power_draw( double starting_temp, double P,double ts);

        void new_OCV(double current);


};


#endif
