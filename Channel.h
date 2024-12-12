#ifndef CHANNEL_H_INCLUDED
#define CHANNEL_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <ctime>
#include <chrono>
//#include <random>
#include <boost/math/distributions/rayleigh.hpp>
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"

using namespace ns3;

double Get_Incidence_Angle_AP_UE(Ptr<Node> AP,Ptr<Node> UE);

double RadtoDeg(double radian);


double DegtoRad(double degree);

double calculateDistance(Ptr<Node> AP,Ptr<Node> UE);


double Estimate_one_VLC_Channel_Gain(Ptr<Node> VLC_AP,Ptr<Node> UE);

double Get_Incidence_Angle_AP_UE(Ptr<Node> AP,Ptr<Node> UE);

double calculate_strong_user_data_rate(double channel_gain ,double strong_user_power);

double calculate_weak_user_VLC_data_rate(double channel_gain ,double strong_user_power, double weak_user_power);

double calculate_RF_data_rate(double channel_gain, Ptr<Node> user1, Ptr<Node> user2);


#endif // CHANNEL_H_INCLUDED
