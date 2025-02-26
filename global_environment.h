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
extern double d_ref;
extern double central_freq;

extern double VLC_field_of_view;
extern double VLC_PHI_half;
extern double VLC_filter_gain;
extern double VLC_concentrator_gain;
extern double VLC_refractive_index;
extern double VLC_receiver_area;
extern double VLC_reflect_efficiency;
extern double VLC_optical_to_electric_factor;
extern double VLC_electric_to_optical_factor;

extern double maximum_current;
extern double minimum_current;
extern double total_power;

extern double kappa;
extern double Nl;
extern double VLC_AP_Popt;
extern double VLC_AP_Bandwidth;

extern double RF_AP_Bandwidth;
extern double Nw;

extern double breakpoint;
extern double dark_saturation_current;
extern double central_carrier_frequency;
extern double angle_of_LOS_arrival;
extern double fill_factor;
extern double thermal_voltage;

extern double minimum_satisfaction;

extern double minimum_data_rate_requirement;
extern double maximum_data_rate_requirement;

extern int power_allocation_formula;
extern int maximum_iteration;
extern int minimum_iteration;
extern int relay_user;


#endif // GLOBAL_ENVIRONMENT_H_INCLUDED
