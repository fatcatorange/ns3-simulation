#include <iostream>
#include <fstream>
#include <string>
#include <chrono> //seed
#include <cmath>
#include <unordered_set>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "Channel.h"
#include "global_environment.h"
#include <boost/math/distributions/rayleigh.hpp>

#include <stdlib.h>
#include <cfloat> // for DBL_MAX
#include <cmath>  // for fabs()
#include "algo2.h"
#include "Hungarian.h"
#include "UE_AP_state.h"

using namespace std;

std::vector<std::pair<double, int>> algo2_sorted_ue_list; //list after sorted by channel gain
std::vector<int> algo2_match(UE / 2, -1);
std::vector<double> algo2_average_rate_matrix(UE, 0);
std::unordered_set<int> excluded_user;
const double INF = std::numeric_limits<double>::max();

int now_user = UE;
int now_pair_count = UE / 2;

double calculate_weak_user_require_power(long double strong_user_power, int weak_user) {
    std::cout<<strong_user_power<<std::endl;
    //strong_user_power *= 1e3;
    long double weak_user_channel_gain = Channel_gain_matrix[weak_user][0];
    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];
    long double weak_user_power = pow(2.0L, (2.0L * now_pair_count * weak_user_minimum_rate / VLC_AP_Bandwidth)) - 1.0;

    long double c = 1.0 / (2.0L * M_PI);
    long double v = pow(VLC_optical_to_electric_factor, 2.0L);
    long double ro = pow(VLC_electric_to_optical_factor, 2.0L);

    long double denominator = (c * powl(v, 2.0L) * powl(ro, 2.0L) * powl(weak_user_channel_gain, 2.0L));
    long double numerator = ((VLC_AP_Bandwidth * Nl * 1e6L) / now_pair_count) + (denominator * strong_user_power);
    weak_user_power *= numerator / denominator;


    std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate<<std::endl;
    return weak_user_power;

}

double algo2_calculate_pair_cost(int weak_user, int strong_user) {
    double strong_user_channel_gain = Channel_gain_matrix[strong_user][0];
    double strong_user_minimum_rate = maximum_requirement_matrix[strong_user] * minimum_satisfaction_matrix[strong_user];
    double strong_user_power = VLC_AP_Bandwidth * Nl * pow(10, 6);

    double c = 1 / (2.0 * M_PI);
    double v = pow(VLC_optical_to_electric_factor, 2);
    double ro = pow(VLC_electric_to_optical_factor, 2);

    strong_user_power/=(now_pair_count * c * v * ro * pow(strong_user_channel_gain, 2) );
    strong_user_power*=(pow(2, (2 * now_pair_count * strong_user_minimum_rate) / VLC_AP_Bandwidth) - 1);


    std::cout<<"check:"<<calculate_strong_user_data_rate(strong_user_channel_gain, strong_user_power)<<" "<<strong_user_minimum_rate<<std::endl;
    double weak_user_power = calculate_weak_user_require_power(strong_user_power, weak_user);
    return strong_user_power + weak_user_power;
}


void algo2_construce_cost_matrix(vector<vector<double>> &cost_matrix) {
    double maxi = 0;
    for(int i = 0; i < now_user / 2;i++) {
        for (int j = 0;j < now_user / 2;j++) {
            cost_matrix[i][j] = algo2_calculate_pair_cost(algo2_sorted_ue_list[i].second, algo2_sorted_ue_list[(now_user / 2) + j].second); // weak, strong
            //std::cout<<cost_matrix[i][j]<<std::endl;
        }
    }
    //print_cost_matrix(cost_matrix);
    //transfer

}

void algo2_print_cost_matrix(std::vector<std::vector<double>> &cost_matrix) {
    std::cout << "cost_matrix" << std::endl;
    for (int i = 0; i < cost_matrix.size(); i++) {
        for (int j = 0; j < cost_matrix[i].size(); j++) {
            //std::cout << std::setw(10) << std::fixed << std::setprecision(4) << cost_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void algo2_user_pairing() {
    for (int i = 0;i < UE;i++) {
        for (int j = 0;j < UE;j++) {
            pairing_matrix[i][j] = 0;
        }
    }
    pairing_matrix.resize(UE, std::vector<int>(UE, 0));
    std::vector<std::vector<double>> cost_matrix(UE / 2, std::vector<double>(UE / 2, 0));
    algo2_construce_cost_matrix(cost_matrix);
    algo2_print_cost_matrix(cost_matrix);
    HungarianAlgorithm HungAlgo;

    HungAlgo.Solve(cost_matrix, algo2_match);

    for (unsigned int x = 0; x < cost_matrix.size(); x++)
		std::cout << x << "," << algo2_match[x] << "\t";
    std::cout<<std::endl;
    for (int i = 0;i < algo2_match.size();i++) {
        int wu = algo2_sorted_ue_list[i].second;
        int su = algo2_sorted_ue_list[algo2_match[i] + (now_user / 2)].second;
        std::cout<<i<<" "<<Channel_gain_matrix[wu][0]<<" "<<algo2_match[i] + (now_user / 2)<<" "<<Channel_gain_matrix[su][0]<<std::endl;
        pairing_matrix[wu][su] = 1;
    }

}

void algo2_list_user() {
    algo2_sorted_ue_list.clear();
    for(int i = 0;i < UE;i++) {
        if (excluded_user.count(i) == 0)
            algo2_sorted_ue_list.push_back({Channel_gain_matrix[i][0], i});
    }
    sort(algo2_sorted_ue_list.begin(), algo2_sorted_ue_list.end());
}

void algorithm2() {
    algo2_list_user();
    algo2_user_pairing();
}
