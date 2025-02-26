#include "global_environment.h"

int UE = 10;
int pairing_count = UE / 2;

double room_size_x = 5;
double room_size_y = 5;
double room_size_z = 3;
double simulation_time = 10.0;

int VLC_AP = 1; // VLC AP number

double d_ref = 10.0; // Reference distance
double central_freq = 5.0e9;

double VLC_field_of_view = 50;
double VLC_PHI_half = 60;
double VLC_filter_gain = 1;
double VLC_concentrator_gain = 1;
double VLC_refractive_index = 1.5;
double VLC_receiver_area = 0.0001;
double VLC_reflect_efficiency = 0.75;
double VLC_optical_to_electric_factor = 0.58;
double VLC_electric_to_optical_factor = 15;

double maximum_current = 0.8; //800mA = 0.8A
double minimum_current = 0.4; //400mA = 0.4A

double total_power = pow((maximum_current - minimum_current) / 2, 2);

double kappa = 0.53;
double Nl = 1e-21;
double VLC_AP_Popt = 1;
double VLC_AP_Bandwidth = 20;

double RF_AP_Bandwidth = 16000000;
double Nw = 3.981071705534986e-21;

double breakpoint = 5;
double dark_saturation_current = 1e-10;
double central_carrier_frequency = 2400000000;
double angle_of_LOS_arrival = 45;
double fill_factor = 0.75;
double thermal_voltage = 0.025;

double minimum_satisfaction = 0.5;


/*
These two variables define the range of maximum data rates.
For example, if a user's minimum satisfaction level is 0.7
and their maximum data rate falls within this range (e.g., 40),
then their minimum satisfied data rate is calculated as:

40 * 0.7 = 28
*/
double minimum_data_rate_requirement = 0;
double maximum_data_rate_requirement = 50;


// 0 = https://ieeexplore.ieee.org/document/9259258 paper propose,
// 1 = chrome-extension://bocbaocobfecmglnmeaeppambideimao/pdf/viewer.html?file=https%3A%2F%2Farxiv.org%2Fpdf%2F2005.09143 paper propose
// 2 = GRPA
// 3 = use formula (22)
int power_allocation_formula = 3;

int maximum_iteration = 100;
int relay_user = 1;

