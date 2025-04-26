#ifndef GLOBAL_ENVIRONMENT_H_INCLUDED
#define GLOBAL_ENVIRONMENT_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

#include "ns3/core-module.h"
#include "ns3/applications-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
//#include "ns3/point-to-point-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/vlc-channel-helper.h"
#include "ns3/vlc-device-helper.h"
#include "ns3/netanim-module.h"

using namespace ns3;

extern int UE;
extern int pairing_count;

extern double room_size_x;
extern double room_size_y;
extern double room_size_z;
extern double simulation_time;

extern int VLC_AP;
extern long double d_ref;
extern long double central_freq;

extern long double VLC_field_of_view;
extern long double VLC_PHI_half;
extern long double VLC_filter_gain;
extern long double VLC_concentrator_gain;
extern long double VLC_refractive_index;
extern long double VLC_receiver_area;
extern long double VLC_reflect_efficiency;
extern long double VLC_optical_to_electric_factor;
extern long double VLC_electric_to_optical_factor;

extern long double maximum_current;
extern long double minimum_current;
extern long double total_power;

extern long double kappa;
extern long double Nl;
extern long double VLC_AP_Popt;
extern long double VLC_AP_Bandwidth;

extern long double RF_AP_Bandwidth;
extern long double Nw;

extern long double breakpoint;
extern long double dark_saturation_current;
extern long double central_carrier_frequency;
extern long double angle_of_LOS_arrival;
extern long double fill_factor;
extern long double thermal_voltage;

extern int IRS_num;
extern int IRS_per_row;
extern double IRS_coefficient;

extern long double minimum_satisfaction;

extern long double minimum_data_rate_requirement;
extern long double maximum_data_rate_requirement;

extern int power_allocation_formula;
extern int exclude_user_method;
extern int maximum_iteration;
extern int minimum_iteration;
extern int relay_user;


#endif // GLOBAL_ENVIRONMENT_H_INCLUDED
