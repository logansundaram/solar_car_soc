#ifndef PACK_TEST_H
#define PACK_TEST_H
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

class Pack{
    std::vector<Battery*> pack;

    public:
        Pack();

        Pack(std::vector<Battery*> pack);

        void pack_charge(double I_L);

        void pack_power_draw(double p_draw);

        void pack_drive_cycle(std::vector<std::pair<double,double>> dri_cy);
};
#endif 

//PACK_TEST_H