#include <iostream>
#include <fstream>
#include <string>
#include <chrono> //seed
#include <cmath>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "UE_AP_state.h"
#include <boost/math/distributions/rayleigh.hpp>
#include "Channel.h"
#include "global_environment.h"
#include "Ref1_algo.h"

std::vector<AP_node*> AP_node_list;
std::vector<UE_node*> UE_node_list;

std::vector<std::vector<double>> Channel_gain_matrix(UE,std::vector<double>(VLC_AP));
std::vector<std::vector<int>> pairing_matrix(UE, std::vector<int>(UE)); // if i, j == 1 => i and j are pairing
std::vector<std::vector<double>> data_rate_matrix(UE,std::vector<double>(VLC_AP));
std::vector<double> power_allocation_matrix(UE, 0); // power for every users
std::vector<int> link_selection_matrix(UE, 0); // 0 = direct link 1 = relay link


void initVLC_AP(NodeContainer &VLC_AP_nodes, std::vector<AP_node*> &AP_list)
{
    MobilityHelper VLC_AP_mobility;
    Ptr<ListPositionAllocator> VLC_AP_pos_list = CreateObject<ListPositionAllocator>();

    VLC_AP_pos_list->Add(Vector(room_size_x / 2,room_size_y / 2,room_size_z));

    VLC_AP_mobility.SetPositionAllocator(VLC_AP_pos_list);
    VLC_AP_mobility.Install(VLC_AP_nodes);

    std::cout<<std::endl;
    std::cout<<"VLC_AP init"<<std::endl;
    std::cout<<std::left;

    for(int i=0;i<VLC_AP;i++)
    {
        Ptr<MobilityModel> VLC_mobility_model = VLC_AP_nodes.Get(i)->GetObject<MobilityModel>();
        Vector pos = VLC_mobility_model->GetPosition();
        AP_list.push_back(new AP_node(i,pos,VLC_AP_nodes.Get(i)));
        std::cout<<"VLC AP ID:"<<AP_list[i]->AP_ID<<std::endl;
        pos = AP_list[i]->node->GetObject<MobilityModel>()->GetPosition();
        std::cout<<"VLC AP POSITION:"<<pos.x<<" "<<pos.y<<" "<<pos.z<<std::endl;
    }
}

void initUE(NodeContainer &UE_nodes,std::vector<UE_node*> &UE_node_list)
{
    MobilityHelper UE_mobility;
    Ptr<ListPositionAllocator> UE_pos_list = CreateObject<ListPositionAllocator>();

    srand(time(NULL));

    double x = room_size_x;
    double y = room_size_y;

    double minx = 0;
    double maxx = x;

    double miny = 0;
    double maxy = y;

    for(int i=0;i<UE;i++)
    {
        double posX = (maxx - minx) * rand() / (RAND_MAX + 1.0) + minx;
        double posY = (maxy - miny) * rand() / (RAND_MAX + 1.0) + miny;

        UE_node_list.push_back(new UE_node(i,Vector(posX,posY,0),UE_nodes.Get(i)));
        UE_pos_list->Add(Vector(posX,posY,0));
    }


    UE_mobility.SetPositionAllocator(UE_pos_list);
    UE_mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    UE_mobility.Install(UE_nodes);

}

void printUE(std::vector<UE_node*> &UE_node_list)
{
    std::cout<<std::endl;
    std::cout<<"print UE_AP"<<std::endl;
    std::cout<<std::left;

    std::cout<<Simulator::Now().GetSeconds()<<std::endl;
    for(int i=0;i<UE;i++)
    {
        Ptr<MobilityModel> UE_mobility_model = UE_node_list[i]->node->GetObject<MobilityModel>();;
        Vector pos = UE_mobility_model->GetPosition();
        std::cout<<" "<<"UE POSITION:"<<i<<" "<<pos.x<<" "<<pos.y<<" "<<pos.z<<std::endl;
    }

}

void basic_init() {

    NodeContainer VLC_AP_node;
    VLC_AP_node.Create(VLC_AP);
    initVLC_AP(VLC_AP_node,AP_node_list);

    NodeContainer UE_nodes;
    UE_nodes.Create(UE);
    initUE(UE_nodes,UE_node_list);

    printUE(UE_node_list);




}

void calculate_throughput() {
    calculate_Channel_Gain_Matrix();
    print_Channel_gain_matrix();

    ref1_algo();
    print_power_allocation_matrix();
    calculate_data_rate_matrix();
    print_data_rate_matrix();
}


void calculate_Channel_Gain_Matrix()
{
    for(int i=0;i < UE;i++)
    {
        for(int j=0;j < VLC_AP;j++)
        {
            Channel_gain_matrix[i][j] = Estimate_one_VLC_Channel_Gain(AP_node_list[j]->node,UE_node_list[i]->node);
        }
    }
}


void print_Channel_gain_matrix() {
    std::cout<<std::endl;
    std::cout<<"channel_gain_matrix in time:"<<" "<<Simulator::Now().GetSeconds()<<std::endl;
    std::cout<<std::left;
    for(int i=0;i < UE;i++)
    {
        for(int j=0;j<VLC_AP;j++)
        {
            std::cout<<std::setw(15)<<Channel_gain_matrix[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}
void count_relay_user() {
    relay_user = 0;
    for(int i = 0;i < UE;i++) {
        if (link_selection_matrix[i] == 1) {
            relay_user++;
        }
    }
}
void calculate_data_rate_matrix() {
    count_relay_user();
    for(int i = 0;i < UE;i++) {
        for (int j = 0;j < UE;j++) {
            if (pairing_matrix[i][j] == 1) {
                std::cout<<i<<" "<<j<<std::endl;
                calculate_pair_data_rate(i, j);
            }
        }
    }
}

void calculate_pair_data_rate(int user1, int user2) {
    if (Channel_gain_matrix[user1][0] > Channel_gain_matrix[user2][0]) { //user 1 is strong user
        std::cout<<"link"<<link_selection_matrix[user2]<<std::endl;
        data_rate_matrix[user1][0] = calculate_strong_user_data_rate(Channel_gain_matrix[user1][0], power_allocation_matrix[user1]);
        if (link_selection_matrix[user2] == 0) { //direct link
            data_rate_matrix[user2][0] = calculate_weak_user_VLC_data_rate(Channel_gain_matrix[user2][0], power_allocation_matrix[user1], power_allocation_matrix[user2]);
        }
        else { // relay link
            data_rate_matrix[user2][0] = std::min(calculate_RF_data_rate(Channel_gain_matrix[user1][0], UE_node_list[user1]->node, UE_node_list[user2]->node),
                                         calculate_weak_user_VLC_data_rate(Channel_gain_matrix[user1][0], power_allocation_matrix[user1], power_allocation_matrix[user2]));
        }

    }
    else { // user 2 is strong user
        std::cout<<"link"<<link_selection_matrix[user1]<<std::endl;
        if (link_selection_matrix[user1] == 0) {
            data_rate_matrix[user1][0] = calculate_weak_user_VLC_data_rate(Channel_gain_matrix[user1][0], power_allocation_matrix[user2], power_allocation_matrix[user1]);
        }
        else {
             data_rate_matrix[user1][0] = std::min(calculate_RF_data_rate(Channel_gain_matrix[user2][0], UE_node_list[user2]->node, UE_node_list[user1]->node),
                                         calculate_weak_user_VLC_data_rate(Channel_gain_matrix[user2][0], power_allocation_matrix[user2], power_allocation_matrix[user1]));
        }
        data_rate_matrix[user2][0] = calculate_strong_user_data_rate(Channel_gain_matrix[user2][0], power_allocation_matrix[user2]);
    }
}

void print_data_rate_matrix() {
    std::cout<<std::endl;
    std::cout<<"data rate_matrix in time:"<<" "<<Simulator::Now().GetSeconds()<<std::endl;
    std::cout<<std::left;
    for(int i = 0;i < UE;i++) {
        for (int j = 0;j < VLC_AP;j++) {
            std::cout<<std::setw(15)<<data_rate_matrix[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
}

void print_power_allocation_matrix() {
    std::cout<<std::endl;
    std::cout<<"power allocation_matrix in time:"<<" "<<Simulator::Now().GetSeconds()<<std::endl;
    std::cout<<std::left;
    for(int i = 0;i < UE;i++) {
        std::cout<<power_allocation_matrix[i]<<std::endl;
    }
}
