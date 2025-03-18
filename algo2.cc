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

double correction_weak_user_power(long double strong_user_power, int weak_user, double origin_power) {
    long double left = origin_power;
    long double weak_user_channel_gain = Channel_gain_matrix[weak_user][0];
    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];

    long double maximum_rate = calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, total_power);

    if (maximum_rate < weak_user_minimum_rate) {
        return total_power;
    }

    long double weak_user_power = pow(2.0L, (2.0L * now_pair_count * maximum_rate / VLC_AP_Bandwidth));


    long double c = 1.0 / (2.0L * M_PI);
    long double v = pow(VLC_optical_to_electric_factor, 2.0L);
    long double ro = pow(VLC_electric_to_optical_factor, 2.0L);
    long double denominator = (c * powl(v, 2.0L) * powl(ro, 2.0L) * powl(weak_user_channel_gain, 2.0L)) ;
    long double numerator = ((VLC_AP_Bandwidth * Nl * 1e6L) / (long double)now_pair_count ) + (denominator * strong_user_power);
    //std::cout<<"nu:"<<numerator<<" "<<denominator<<std::endl;
    weak_user_power *= numerator;
    weak_user_power /= denominator;
    weak_user_power -= (numerator / denominator);

    long double right = weak_user_power;
    while(right - left > 1e-15) {
        long double mid = (right + left) / 2;
        long double now_rate = calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, mid);
        if (now_rate > weak_user_minimum_rate) {
            weak_user_power = mid;
            right = mid;
        }
        else {
            left = mid;
        }

    }
    //std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate <<std::endl;

    return weak_user_power;
}

double correction_weak_user_RF_power(long double strong_user_power, int weak_user, int strong_user, double origin_power) {
    long double left = origin_power;
    long double strong_user_channel_gain = Channel_gain_matrix[strong_user][0];
    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];

    long double maximum_rate = calculate_RF_data_rate(strong_user_power, UE_node_list[strong_user]->node, UE_node_list[weak_user]->node);;

    if (maximum_rate < weak_user_minimum_rate) {
        return total_power;
    }

    long double weak_user_power = pow(2.0L, (2.0L * now_pair_count * maximum_rate / VLC_AP_Bandwidth));


    long double c = 1.0 / (2.0L * M_PI);
    long double v = pow(VLC_optical_to_electric_factor, 2.0L);
    long double ro = pow(VLC_electric_to_optical_factor, 2.0L);
    long double denominator = (c * powl(v, 2.0L) * powl(ro, 2.0L) * powl(strong_user_channel_gain, 2.0L)) ;
    long double numerator = ((VLC_AP_Bandwidth * Nl * 1e6L) / (long double)now_pair_count ) + (denominator * strong_user_power);
    //std::cout<<"nu:"<<numerator<<" "<<denominator<<std::endl;
    weak_user_power *= numerator;
    weak_user_power /= denominator;
    weak_user_power -= (numerator / denominator);

    long double right = weak_user_power;
    while(right - left > 1e-15) {
        long double mid = (right + left) / 2;
        long double now_rate = calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, mid);
        if (now_rate > weak_user_minimum_rate) {
            weak_user_power = mid;
            right = mid;
        }
        else {
            left = mid;
        }

    }
    //std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate <<std::endl;

    return weak_user_power;
}

double calculate_weak_user_require_power(long double strong_user_power, int weak_user) {
    now_pair_count *= 20;
    long double weak_user_channel_gain = Channel_gain_matrix[weak_user][0];
    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];
    long double weak_user_power = pow(2.0L, (2.0L * now_pair_count * weak_user_minimum_rate / VLC_AP_Bandwidth));
    //std::cout<<(2.0L * now_pair_count * weak_user_minimum_rate / VLC_AP_Bandwidth)<<std::endl;
    long double c = 1.0 / (2.0L * M_PI);
    long double v = pow(VLC_optical_to_electric_factor, 2.0L);
    long double ro = pow(VLC_electric_to_optical_factor, 2.0L);
    //std::cout<<"wc:"<<weak_user_channel_gain<<" "<<pow(weak_user_channel_gain, 2)<<std::endl;
    long double denominator = (c * powl(v, 2.0L) * powl(ro, 2.0L) * powl(weak_user_channel_gain, 2.0L)) ;
    long double numerator = ((VLC_AP_Bandwidth * Nl * 1e6L) / (long double)now_pair_count / 20) + (denominator * strong_user_power);
    //std::cout<<"nu:"<<numerator<<" "<<denominator<<std::endl;
    weak_user_power *= numerator;
    weak_user_power /= pow(2.0L, (2.0L * (long double)now_pair_count * weak_user_minimum_rate * 0.95L / VLC_AP_Bandwidth));
    weak_user_power /= denominator;
    weak_user_power -= (numerator / denominator);

    //std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate / 2 <<std::endl;
    now_pair_count = now_user / 2;

    return correction_weak_user_power(strong_user_power, weak_user, weak_user_power);

}

double calculate_weak_user_RF_require_power(long double strong_user_power, int weak_user, int strong_user, int hybrid_user_count) {
    long double strong_user_channel_gain = Channel_gain_matrix[strong_user][0];
    long double RF_upper_bound = calculate_RF_data_rate(strong_user_power, UE_node_list[strong_user]->node, UE_node_list[weak_user]->node);

    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];
    //std::cout<<"RF upper bound: "<<relay_user<<" "<<RF_upper_bound<<std::endl;
    // can't satisfy weak user
    if (weak_user_minimum_rate > RF_upper_bound) {

        return -1;
    }

    long double weak_user_power = pow(2.0L, (2.0L * now_pair_count * weak_user_minimum_rate / VLC_AP_Bandwidth)) - 1.0;

    long double c = 1.0 / (2.0L * M_PI);
    long double v = pow(VLC_optical_to_electric_factor, 2.0L);
    long double ro = pow(VLC_electric_to_optical_factor, 2.0L);

    long double denominator = (c * powl(v, 2.0L) * powl(ro, 2.0L) * powl(strong_user_channel_gain, 2.0L));
    long double numerator = ((VLC_AP_Bandwidth * Nl * 1e6L) / now_pair_count) + (denominator * strong_user_power);
    weak_user_power *= numerator / denominator;



    //std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate<<std::endl;
    return correction_weak_user_RF_power(strong_user_power, weak_user, strong_user, weak_user_power);

}

double calculate_strong_user_require_power(int strong_user) {
    double strong_user_channel_gain = Channel_gain_matrix[strong_user][0];
    double strong_user_minimum_rate = maximum_requirement_matrix[strong_user] * minimum_satisfaction_matrix[strong_user];
    double strong_user_power = VLC_AP_Bandwidth * Nl * pow(10, 6);

    double c = 1 / (2.0 * M_PI);
    double v = pow(VLC_optical_to_electric_factor, 2);
    double ro = pow(VLC_electric_to_optical_factor, 2);

    strong_user_power/=(now_pair_count * c * v * ro * pow(strong_user_channel_gain, 2) );
    strong_user_power*=(pow(2, (2 * now_pair_count * strong_user_minimum_rate) / VLC_AP_Bandwidth) - 1);

    return strong_user_power;
}

double algo2_calculate_pair_cost(int weak_user, int strong_user) {

    double strong_user_power = calculate_strong_user_require_power(strong_user);

    //std::cout<<"check:"<<calculate_strong_user_data_rate(strong_user_channel_gain, strong_user_power)<<" "<<strong_user_minimum_rate<<std::endl;
    double weak_user_power = calculate_weak_user_require_power(strong_user_power, weak_user);
    return strong_user_power + weak_user_power;
}


void algo2_construce_cost_matrix(vector<vector<double>> &cost_matrix) {
    double maxi = 0;
    for(int i = 0; i < now_user / 2;i++) {
        relay_user = i;
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
        //std::cout<<i<<" "<<Channel_gain_matrix[wu][0]<<" "<<algo2_match[i] + (now_user / 2)<<" "<<Channel_gain_matrix[su][0]<<std::endl;
        pairing_matrix[wu][su] = 1;
    }

}

double calculate_save_power(int weak_user, int strong_user, int hybrid_user_count) {
    double strong_user_power = calculate_strong_user_require_power(strong_user);
    double weak_user_power_using_VLC = calculate_weak_user_require_power(strong_user_power, weak_user);
    double weak_user_power_using_RF = calculate_weak_user_RF_require_power(strong_user_power, weak_user, strong_user, hybrid_user_count);

    if (weak_user_power_using_RF > weak_user_power_using_VLC || weak_user_power_using_RF < 0) {
        return -1;
    }
    //std::cout<<"VLC: "<<weak_user_power_using_VLC<<" RF: "<<weak_user_power_using_RF<<std::endl;
    return weak_user_power_using_VLC - weak_user_power_using_RF;
}

void algo2_link_selection() {
    vector<vector<int>> algo2_link_selection_table(now_user / 2 + 1, vector<int>(now_user / 2, 0));
    for(int hybrid_user = 0 ;hybrid_user < algo2_link_selection_table.size();hybrid_user++) {
        priority_queue<pair<double, int>> pq;
        relay_user = hybrid_user;
        for(int i = 0;i < now_user / 2;i++) {
            int wu = algo2_sorted_ue_list[i].second;
            int su = algo2_sorted_ue_list[algo2_match[i] + (now_user / 2)].second;
            double saved_power = calculate_save_power(wu, su, hybrid_user);
            //std::cout<<"saved:"<<saved_power<<std::endl;
            pq.push({saved_power,i});
        }

        for(int i = 0;i < hybrid_user;i++) {
            int user = pq.top().second;

            algo2_link_selection_table[hybrid_user][user] = (pq.top().first > 0) ? 1 : -1;
            pq.pop();
        }
    }

    double maximum_saved_power = 0;
    int maximum_saved_index = 0;
    for(int i = 1;i < now_user / 2 + 1;i++) {
        relay_user = i;
        double now_saved_power = 0;
        for(int j = 0;j < algo2_link_selection_table[i].size();j++) {
            if (algo2_link_selection_table[i][j] < 0) {
                now_saved_power = -1;
                break;
            }
            else if (algo2_link_selection_table[i][j] == 1) {
                int wu = algo2_sorted_ue_list[j].second;
                int su = algo2_sorted_ue_list[algo2_match[j] + (now_user / 2)].second;
                double saved_power = calculate_save_power(wu, su, j);
                if (saved_power < 0) {
                    now_saved_power = -1;
                    break;
                }
                now_saved_power+=saved_power;
            }
        }

        //std::cout<<"saved power: "<<i<<" "<<now_saved_power<<std::endl;
        if (now_saved_power > maximum_saved_power) {
            maximum_saved_index = i;
            maximum_saved_power = now_saved_power;
        }
    }

    //std::cout<<"maximum index"<<maximum_saved_index<<std::endl;
    relay_user = maximum_saved_index;

    for (int i = 0;i < algo2_link_selection_table[maximum_saved_index].size();i++) {
        if (algo2_link_selection_table[maximum_saved_index][i] == 1) {
            link_selection_matrix[algo2_sorted_ue_list[i].second] = 1;
        }
    }

    for(int i = 0;i < link_selection_matrix.size();i++) {
        //cout<<link_selection_matrix[i]<<" ";
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

//return: satisfied user count
int algo2_power_allocation(double &now_power) {

    for(int i = 0;i < algo2_match.size();i++) {
        int wu = algo2_sorted_ue_list[i].second;
        int su = algo2_sorted_ue_list[algo2_match[i] + (now_user / 2)].second;
        double strong_user_power = calculate_strong_user_require_power(su);
        long double weak_user_power = link_selection_matrix[i] == 0 ? calculate_weak_user_require_power(strong_user_power, wu) :
            calculate_weak_user_RF_require_power(strong_user_power, wu, su, relay_user);
        long double require_power = strong_user_power + weak_user_power;
        //std::cout<<require_power<<" "<<total_power<<std::endl;
        if (now_power >= require_power) {
            now_power-=require_power;
            power_allocation_matrix[wu] = weak_user_power;
            power_allocation_matrix[su] = strong_user_power;

            /*
            long double weak_user_minimum_rate = maximum_requirement_matrix[wu] * minimum_satisfaction_matrix[wu];
            if (link_selection_matrix[i] == 0)
                std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(Channel_gain_matrix[wu][0], strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate<<std::endl;
            else
                std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(Channel_gain_matrix[su][0], strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate<<std::endl;
            */

        }
        else {
            return i + 1;
        }
    }
    return algo2_match.size();
}

void IRS_assignment(double now_power, int satisfied_user) {

}

void algorithm2() {
    algo2_list_user();
    algo2_user_pairing();
    algo2_link_selection();

    for(int i = 0;i < link_selection_matrix.size();i++) {
        std::cout<<link_selection_matrix[i]<<" ";
    }
    std::cout<<std::endl;

    double now_power = total_power;
    int satisfied_user = algo2_power_allocation(now_power);


}
