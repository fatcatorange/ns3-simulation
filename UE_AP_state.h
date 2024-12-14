#ifndef UE_AP_STATE_H_INCLUDED
#define UE_AP_STATE_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <iomanip>


#include "ns3/core-module.h"
#include "ns3/applications-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
//#include "ns3/point-to-point-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/vlc-channel-helper.h"
#include "ns3/vlc-device-helper.h"
#include "ns3/netanim-module.h"
#include "AP_Node.h"
#include "UE_Node.h"

void basic_init();
void initVLC_AP(NodeContainer &VLC_AP_nodes, std::vector<AP_node*> &AP_list);
void initUE(NodeContainer &UE_nodes,std::vector<UE_node*> &UE_node_list);
void printUE(std::vector<UE_node*> &UE_node_list);


void calculate_throughput();
void calculate_Channel_Gain_Matrix();
void print_Channel_gain_matrix();
void calculate_data_rate_matrix();
void calculate_pair_data_rate(int user1, int user2);
void print_data_rate_matrix();
void print_power_allocation_matrix();

extern std::vector<AP_node*> AP_node_list;
extern std::vector<UE_node*> UE_node_list;

extern std::vector<std::vector<double>> Channel_gain_matrix;
extern std::vector<std::vector<int>> pairing_matrix;
extern std::vector<std::vector<double>> data_rate_matrix;
extern std::vector<double> power_allocation_matrix;
extern std::vector<int> link_selection_matrix;




#endif // UE_AP_STATE_H_INCLUDED
