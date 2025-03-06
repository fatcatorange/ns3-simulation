#include <iostream>
#include <fstream>
#include <string>
#include <chrono> //seed
#include <cmath>
#include <random>
#include <unordered_set>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "UE_AP_state.h"
#include <boost/math/distributions/rayleigh.hpp>
#include "Channel.h"
#include "global_environment.h"
#include "Ref1_algo.h"
#include "algo2.h"

std::vector<AP_node*> AP_node_list;
std::vector<UE_node*> UE_node_list;

std::vector<std::vector<double>> Channel_gain_matrix(UE,std::vector<double>(VLC_AP));
std::vector<std::vector<int>> pairing_matrix(UE, std::vector<int>(UE)); // if i, j == 1 => i and j are pairing
std::vector<std::vector<double>> data_rate_matrix(UE,std::vector<double>(VLC_AP));
std::vector<double> power_allocation_matrix(UE, (total_power / UE)); // power for every users
std::vector<int> link_selection_matrix(UE, 0); // 0 = direct link 1 = relay link
std::vector<double> minimum_satisfaction_matrix(UE, 0);
std::vector<double> maximum_requirement_matrix(UE, 0);


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

    unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    srand(seed);

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

void generate_user_requirement() {
    std::random_device rd;   // 真正的隨機種子
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> maximum_requirement_distribution(minimum_data_rate_requirement, maximum_data_rate_requirement);
    std::uniform_real_distribution<double> minimum_satisfaction_distribution(minimum_satisfaction, 1.0);

    //generate maximum requirement
    for(int i = 0;i < UE;i++) {
        double random_value = maximum_requirement_distribution(gen);
        maximum_requirement_matrix[i] = random_value;
    }

    for(int i = 0;i < UE;i++) {
        double random_value = minimum_satisfaction_distribution(gen);
        minimum_satisfaction_matrix[i] = random_value;
    }
}

void print_user_requirement() {
    for(int i = 0;i < UE;i++) {
        std::cout<<"user "<<i<<" maximum requirement:"<<" "<<maximum_requirement_matrix[i]<<" minimum satisfaction: "<<minimum_satisfaction_matrix[i]<<std::endl;
    }
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
    generate_user_requirement();
    print_user_requirement();

    /*
    init_ref1_algo();
    algorithm2();
    */


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
                calculate_pair_data_rate(i, j);
            }
        }
    }
}

void calculate_pair_data_rate(int user1, int user2) {
    if (Channel_gain_matrix[user1][0] > Channel_gain_matrix[user2][0]) { //user 1 is strong user
        //std::cout<<"link"<<link_selection_matrix[user2]<<std::endl;
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
        //std::cout<<"link"<<link_selection_matrix[user1]<<std::endl;
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

double calculate_sum_rate() {
    double sum_rate = 0;
    for(int i = 0;i < data_rate_matrix.size();i++) {
        for(int j = 0;j <data_rate_matrix[i].size();j++) {
            sum_rate+=data_rate_matrix[i][j];
        }

    }

    return sum_rate;
}

double calculate_fairness() {
    double jann_fairness = 0;
    double top = 0;
    double down = 0;
    for(int i = 0;i < UE;i++) {
        double now_rate = data_rate_matrix[i][0];
        top+=now_rate;
        down+=pow(now_rate, 2);
    }

    top = pow(top, 2);
    down = 2 * pairing_count * down;

    return top / down;
}

double calculate_pair_fairness() {
    double jann_fairness = 0;
    double top = 0;
    double down = 0;
    std::unordered_set<int> visited;
    for(int i = 0;i < pairing_matrix.size();i++) {
        for(int j = 0;j < pairing_matrix[i].size();j++) {
            if (pairing_matrix[i][j] == 1 && visited.count(i) == 0 && visited.count(j) == 0) {
                double now_rate = data_rate_matrix[i][0] + data_rate_matrix[j][0];
                //std::cout<<"now rate:"<<i<<" "<<j<<" "<<now_rate<<std::endl;
                top+=now_rate;
                down+=pow(now_rate, 2);
                visited.insert(i);
                visited.insert(j);
            }
        }
    }

    top = pow(top, 2);
    down = pairing_count * down;

    return top / down;
}

double throughput_write_file() {
    std::fstream outFile;
    outFile.open("/home/jimmy/repos/ns-3-allinone/ns-3.25/scratch/thesis/output.csv",std::ios::out|std::ios::app);
    if(outFile.is_open() == false)
    {
        std::cout<<"file not open"<<std::endl;
    }
    else
    {
        //outFile<<std::endl;
        outFile<<calculate_sum_rate()<<','<<calculate_pair_fairness()<<','<<calculate_fairness()<<std::endl;//throughput
        //std::cout<<calculate_sum_rate()<<','<<calculate_pair_fairness()<<','<<calculate_fairness()<<std::endl;
        //outFile<<std::endl;
    }
    outFile.close();
}




void pair_fairness_write_file() {
    //std::cout<<"fairness:"<<calculate_pair_fairness()<<std::endl;
}
