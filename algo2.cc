#include <iostream>
#include <fstream>
#include <string>
#include <chrono> //seed
#include <cmath>
#include <unordered_set>
#include <tuple>

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
std::ofstream logFile("log.txt");

bool all_satisfied = false;
int now_user = UE;
int now_pair_count = UE / 2;
bool have_psuedo_user = UE % 2 == 1;

double correction_weak_user_power(long double strong_user_power, int weak_user, double origin_power) {
    long double left = origin_power;
    long double weak_user_channel_gain = Channel_gain_matrix[weak_user][0] + IRS_channel_gain_matrix[weak_user];
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

    long double right = total_power;
    while(right - left > 1e-20) {
        long double mid = (right + left) / 2;
        long double now_rate = calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, mid);
        //logFile<<"now rate:"<<now_rate<<std::endl;
        if (now_rate > weak_user_minimum_rate) {
            weak_user_power = mid;
            right = mid;
        }
        else {
            left = mid;
        }

    }
    //logFile<<"strong power:"<<strong_user_power<<" weak user power:"<<weak_user_power<<" weak user channel gain:"<<weak_user_channel_gain<<std::endl;
    //logFile<<"check weak:" <<calculate_weak_user_VLC_data_rate(weak_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate <<std::endl;

    return weak_user_power;
}

double correction_weak_user_RF_power(long double strong_user_power, int weak_user, int strong_user, double origin_power) {
    long double left = origin_power;
    long double strong_user_channel_gain = Channel_gain_matrix[strong_user][0] + IRS_channel_gain_matrix[strong_user];
    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];

    long double maximum_rate = calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, total_power);
    //std::cout<<maximum_rate<<std::endl;
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
    while(right - left > 1e-20) {
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

    //std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate <<std::endl;


    return weak_user_power;
}

double calculate_weak_user_require_power(long double strong_user_power, int weak_user) {
    if (weak_user == -1) {
        return 0;
    }
    now_pair_count *= 20;
    long double weak_user_channel_gain = Channel_gain_matrix[weak_user][0] + IRS_channel_gain_matrix[weak_user];
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
    if (weak_user == -1) {
        return 0;
    }
    long double strong_user_channel_gain = Channel_gain_matrix[strong_user][0] + IRS_channel_gain_matrix[strong_user];
    long double RF_upper_bound = calculate_RF_data_rate(strong_user_power, UE_node_list[strong_user]->node, UE_node_list[weak_user]->node);

    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];
    //std::cout<<"RF upper bound: "<<relay_user<<" "<<RF_upper_bound<<std::endl;
    // can't satisfy weak user
    if (weak_user_minimum_rate > RF_upper_bound) {

        return total_power;
    }

    long double weak_user_power = pow(2.0L, (2.0L * now_pair_count * weak_user_minimum_rate / VLC_AP_Bandwidth)) - 1.0;

    long double c = 1.0 / (2.0L * M_PI);
    long double v = pow(VLC_optical_to_electric_factor, 2.0L);
    long double ro = pow(VLC_electric_to_optical_factor, 2.0L);

    long double denominator = (c * powl(v, 2.0L) * powl(ro, 2.0L) * powl(strong_user_channel_gain, 2.0L));
    long double numerator = ((VLC_AP_Bandwidth * Nl * 1e6L) / now_pair_count) + (denominator * strong_user_power);
    weak_user_power *= numerator / denominator;



   // std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate<<std::endl;
    return correction_weak_user_RF_power(strong_user_power, weak_user, strong_user, weak_user_power);

}


double bsrpa_correction_weak_user_RF_power(long double strong_user_power, int weak_user, int strong_user, double origin_power) {
    long double left = origin_power;
    long double strong_user_channel_gain = Channel_gain_matrix[strong_user][0] + IRS_channel_gain_matrix[strong_user];
    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];
    long double RF_upper_bound = calculate_RF_data_rate(strong_user_power, UE_node_list[strong_user]->node, UE_node_list[weak_user]->node);
    long double maximum_rate = calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, total_power);
    //std::cout<<maximum_rate<<std::endl;
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

    long double right = total_power;
    while(right - left > 1e-20) {
        long double mid = (right + left) / 2;
        long double now_rate = calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, mid);
        if (now_rate > std::min(weak_user_minimum_rate, RF_upper_bound)) {
            weak_user_power = mid;
            right = mid;
        }
        else {
            left = mid;
        }

    }

    //std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate <<std::endl;


    return weak_user_power;
}

double bsrpa_calculate_weak_user_RF_require_power(long double strong_user_power, int weak_user, int strong_user, int hybrid_user_count) {
    if (weak_user == -1) {
        return 0;
    }
    long double strong_user_channel_gain = Channel_gain_matrix[strong_user][0] + IRS_channel_gain_matrix[strong_user];
    long double RF_upper_bound = calculate_RF_data_rate(strong_user_power, UE_node_list[strong_user]->node, UE_node_list[weak_user]->node);

    long double weak_user_minimum_rate = maximum_requirement_matrix[weak_user] * minimum_satisfaction_matrix[weak_user];
    //std::cout<<"RF upper bound: "<<relay_user<<" "<<RF_upper_bound<<std::endl;

    long double weak_user_power = pow(2.0L, (2.0L * now_pair_count * weak_user_minimum_rate / VLC_AP_Bandwidth)) - 1.0;

    long double c = 1.0 / (2.0L * M_PI);
    long double v = pow(VLC_optical_to_electric_factor, 2.0L);
    long double ro = pow(VLC_electric_to_optical_factor, 2.0L);

    long double denominator = (c * powl(v, 2.0L) * powl(ro, 2.0L) * powl(strong_user_channel_gain, 2.0L));
    long double numerator = ((VLC_AP_Bandwidth * Nl * 1e6L) / now_pair_count) + (denominator * strong_user_power);
    weak_user_power *= numerator / denominator;



   // std::cout<<"check weak:" <<calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate<<std::endl;
    return bsrpa_correction_weak_user_RF_power(strong_user_power, weak_user, strong_user, weak_user_power);

}

double calculate_strong_user_require_power(int strong_user) {
    double strong_user_channel_gain = Channel_gain_matrix[strong_user][0] + IRS_channel_gain_matrix[strong_user];
    double strong_user_minimum_rate = maximum_requirement_matrix[strong_user] * minimum_satisfaction_matrix[strong_user];
    double strong_user_power = VLC_AP_Bandwidth * Nl * pow(10, 6);

    double c = 1 / (2.0 * M_PI);
    double v = pow(VLC_optical_to_electric_factor, 2);
    double ro = pow(VLC_electric_to_optical_factor, 2);

    strong_user_power/=(now_pair_count * c * v * ro * pow(strong_user_channel_gain, 2) );
    strong_user_power*=(pow(2, (2 * now_pair_count * strong_user_minimum_rate) / VLC_AP_Bandwidth) - 1);

   // std::cout<<"check strong:"<<calculate_strong_user_data_rate(strong_user_channel_gain, strong_user_power)<<std::endl;

    return strong_user_power;
}

double algo2_calculate_pair_cost(int weak_user, int strong_user) {

    double strong_user_power = calculate_strong_user_require_power(strong_user);

    //std::cout<<"check:"<<calculate_strong_user_data_rate(strong_user_channel_gain, strong_user_power)<<" "<<strong_user_minimum_rate<<std::endl;
    double weak_user_power = link_selection_matrix[weak_user] == 0 ? calculate_weak_user_require_power(strong_user_power, strong_user) :
            calculate_weak_user_RF_require_power(strong_user_power, weak_user, strong_user, relay_user);
    //std::cout<<weak_user<<" "<<strong_user<<std::endl;
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
    //std::cout << "cost_matrix" << std::endl;
    for (int i = 0; i < cost_matrix.size(); i++) {
        for (int j = 0; j < cost_matrix[i].size(); j++) {
            //std::cout << std::setw(10) << std::fixed << std::setprecision(4) << cost_matrix[i][j] << " ";
        }
    }
}

void algo2_user_pairing() {
    for (int i = 0;i < UE;i++) {
        for (int j = 0;j < UE;j++) {
            pairing_matrix[i][j] = 0;
        }
    }
    pairing_matrix.clear();
    pairing_matrix.resize(UE, std::vector<int>(UE, 0));
    std::vector<std::vector<double>> cost_matrix(now_user / 2, std::vector<double>(now_user / 2, 0));
    algo2_construce_cost_matrix(cost_matrix);

    algo2_print_cost_matrix(cost_matrix);
    HungarianAlgorithm HungAlgo;

    HungAlgo.Solve(cost_matrix, algo2_match);

    for (int i = 0;i < algo2_match.size();i++) {
        int wu = algo2_sorted_ue_list[i].second;
        int su = algo2_sorted_ue_list[algo2_match[i] + (now_user / 2)].second;
        //std::cout<<i<<" "<<Channel_gain_matrix[wu][0]<<" "<<algo2_match[i] + (now_user / 2)<<" "<<Channel_gain_matrix[su][0]<<std::endl;
        if (wu != -1)
            pairing_matrix[wu][su] = 1;
        else {
            pairing_matrix[su][su] = 1;
        }
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
            //std::cout<<i<<std::endl;
            int wu = algo2_sorted_ue_list[i].second;
            int su = algo2_sorted_ue_list[algo2_match[i] + (now_user / 2)].second;
            //std::cout<<wu<<" "<<su<<std::endl;
            double saved_power = calculate_save_power(wu, su, hybrid_user);
            //std::cout<<"saved:"<<saved_power<<std::endl;
            pq.push({saved_power,i});
        }

        for(int i = 0;i < hybrid_user;i++) {
            int user = pq.top().second;
            if (user != -1)
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
        if (algo2_link_selection_table[maximum_saved_index][i] == 1 && algo2_sorted_ue_list[i].second != -1) {
            link_selection_matrix[algo2_sorted_ue_list[i].second] = 1;
        }
    }


}

void algo2_list_user() {
    algo2_sorted_ue_list.clear();
    now_user = 0;
    for(int i = 0;i < UE;i++) {
        if (excluded_user.count(i) == 0) {
            algo2_sorted_ue_list.push_back({Channel_gain_matrix[i][0], i});
            //std::cout<<"add:"<<i<<std::endl;
            now_user++;
        }

    }
    //-1 is pseudo user
    have_psuedo_user = false;
    if (now_user % 2 == 1) {
        have_psuedo_user = true;
        now_user++;
        algo2_sorted_ue_list.push_back({0, -1});
    }

    sort(algo2_sorted_ue_list.begin(), algo2_sorted_ue_list.end());
    now_pair_count = (now_user) / 2;
    pairing_count = now_pair_count;
    //std::cout<<"now user:"<<now_user<<std::endl;
}

//return: satisfied user pair count
//sorted_user pair: require power, weak user, strong user (sorted by pair require power)
int algo2_power_allocation(double &now_power, std::vector<std::tuple<double, int, int>> &sorted_user_pair) {

    sorted_user_pair.clear();
    for(int i = 0;i < algo2_match.size();i++) {
        int wu = algo2_sorted_ue_list[i].second;
        int su = algo2_sorted_ue_list[algo2_match[i] + (now_user / 2)].second;
        //std::cout<<"pa user:"<<wu<<" "<<su<<std::endl;
        double strong_user_power = calculate_strong_user_require_power(su);
        //std::cout<<"strong user:"<<strong_user_power<<std::endl;
        long double weak_user_power = link_selection_matrix[wu] == 0 ? calculate_weak_user_require_power(strong_user_power, wu) :
            calculate_weak_user_RF_require_power(strong_user_power, wu, su, relay_user);
        long double require_power = strong_user_power + weak_user_power;
        if (wu != -1)
            power_allocation_matrix[wu] = weak_user_power;
        power_allocation_matrix[su] = strong_user_power;
        //std::cout<<require_power<<" "<<total_power<<std::endl;

        sorted_user_pair.push_back(tuple<double, int, int>(require_power, wu, su));
        //std::get<0>(sorted_user_pair[i]);

    }

    sort(sorted_user_pair.begin(), sorted_user_pair.end());

    for(int i = 0;i < sorted_user_pair.size();i++) {
        double now_user_require_power = std::get<0>(sorted_user_pair[i]);
        if (now_power >= now_user_require_power) {
            now_power-= now_user_require_power;
            //std::cout<<now_power<<" "<<now_user_require_power<<std::endl;
        }
        else {
            for(int j = i;j < sorted_user_pair.size();j++) {
                int wu = std::get<1>(sorted_user_pair[j]);
                int su = std::get<2>(sorted_user_pair[j]);
                if (wu != -1)
                    power_allocation_matrix[wu] = 0;
                power_allocation_matrix[su] = 0;
                //std::cout<<wu<<" "<<su<<std::endl;
            }
            return i;
        }

    }

    all_satisfied = true;
    return algo2_match.size();
}

//check IRS should assign to weak or strong user, and return extra power;
double assign_IRS_element(int wu, int su, int IRS_number) {
    if (wu == -1) {
        double strong_user_IRS_bonus = Estimate_IRS_channel_gain(AP_node_list[0]->node, UE_node_list[su]->node, IRS_nodes.Get(IRS_number));
        IRS_channel_gain_matrix[su]+= strong_user_IRS_bonus;
        double strong_user_power = calculate_strong_user_require_power(su);
        double extra_power = power_allocation_matrix[su] - strong_user_power;
        power_allocation_matrix[su] = strong_user_power;
        return extra_power;
    }

    double weak_user_total_channel_gain = Channel_gain_matrix[wu][0] + IRS_channel_gain_matrix[wu];
    double strong_user_total_channel_gain = Channel_gain_matrix[su][0] + IRS_channel_gain_matrix[su];

    double weak_user_IRS_bonus = Estimate_IRS_channel_gain(AP_node_list[0]->node, UE_node_list[wu]->node, IRS_nodes.Get(IRS_number));
    double strong_user_IRS_bonus = Estimate_IRS_channel_gain(AP_node_list[0]->node, UE_node_list[su]->node, IRS_nodes.Get(IRS_number));

    //std::cout<<"origin:"<<weak_user_total_channel_gain<<" "<<strong_user_total_channel_gain<<std::endl;
    //std::cout<<"bonus: "<<weak_user_IRS_bonus<<" "<<strong_user_IRS_bonus<<std::endl;

    IRS_channel_gain_matrix[wu]+= weak_user_IRS_bonus;

    double new_weak_user_power = algo2_calculate_pair_cost(wu, su);

    IRS_channel_gain_matrix[wu]-= weak_user_IRS_bonus;

    IRS_channel_gain_matrix[su]+= strong_user_IRS_bonus;

    double new_strong_user_power = algo2_calculate_pair_cost(wu, su);

    IRS_channel_gain_matrix[su]-= strong_user_IRS_bonus;

    double origin_power = algo2_calculate_pair_cost(wu, su);
    //std::cout<<"origin:"<<origin_power<<" "<<new_weak_user_power<<" "<<new_strong_user_power<<std::endl;

    //order will not rearrange
    if (weak_user_total_channel_gain + weak_user_IRS_bonus < strong_user_total_channel_gain && new_weak_user_power < new_strong_user_power) {
        //std::cout<<"give weak!"<<std::endl;
        IRS_channel_gain_matrix[wu]+= weak_user_IRS_bonus;
    }
    else {
        //std::cout<<"give strong!"<<std::endl;
        IRS_channel_gain_matrix[su]+= strong_user_IRS_bonus;
    }

    double strong_user_power = calculate_strong_user_require_power(su);
    double weak_user_power = link_selection_matrix[wu] == 0 ? calculate_weak_user_require_power(strong_user_power, wu) :
            calculate_weak_user_RF_require_power(strong_user_power, wu, su, relay_user);
    double extra_power = (power_allocation_matrix[su] - strong_user_power) + (power_allocation_matrix[wu] - weak_user_power);
    power_allocation_matrix[su] = strong_user_power;
    power_allocation_matrix[wu] = weak_user_power;
    double strong_user_channel_gain = Channel_gain_matrix[wu][0] + IRS_channel_gain_matrix[wu];
    double weak_user_minimum_rate = maximum_requirement_matrix[wu] * minimum_satisfaction_matrix[wu];
    //std::cout<<"check normal weak:" <<calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, weak_user_power)<<" "<<weak_user_minimum_rate <<std::endl;
    //std::cout<<"extra power:"<<extra_power<<std::endl;
    return extra_power;
}

void IRS_assignment(double &now_power, int satisfied_pair, std::vector<std::tuple<double, int, int>> &sorted_user_pair) {
    //0 = all satisfied 1 = strong user, 2 = weak user
    int next_user = satisfied_pair == now_pair_count ? 0 : 1;
    int service_pair = satisfied_pair;
    int first_su = next_user == 0 ? 0 : std::get<2>(sorted_user_pair[satisfied_pair]);
    double next_user_requirement = calculate_strong_user_require_power(first_su);

    if (service_pair == 0) {
        return;
    }
    for(int i = 0;i < IRS_num;i++) {
        int now_assign_pair = i % service_pair;

        int wu = std::get<1>(sorted_user_pair[now_assign_pair]);
        int su = std::get<2>(sorted_user_pair[now_assign_pair]);


        //std::cout<<"user:"<<wu<<" "<<su<<std::endl;

        if (next_user != 2 || now_assign_pair != service_pair - 1) {
            now_power+=assign_IRS_element(wu, su, i);
        }
        else {
            now_power+=assign_IRS_element(-1, su, i);
        }

        if (next_user == 1) {
            if (now_power > next_user_requirement) {
                now_power-=next_user_requirement;
                int now_strong_user = std::get<2>(sorted_user_pair[service_pair ]);
                int next_weak_user = std::get<1>(sorted_user_pair[service_pair]);
                service_pair++;
                power_allocation_matrix[now_strong_user] = next_user_requirement;
                if (next_weak_user == -1) {
                    satisfied_pair++;
                    if (satisfied_pair != now_pair_count) {
                        int next_strong_user = std::get<1>(sorted_user_pair[service_pair]);
                        next_user_requirement = calculate_strong_user_require_power(next_strong_user);
                    }
                }
                else {
                    next_user_requirement = link_selection_matrix[next_weak_user] == 0 ? calculate_weak_user_require_power(next_user_requirement, next_weak_user) :
                            calculate_weak_user_RF_require_power(next_user_requirement, next_weak_user, now_strong_user, relay_user);
                    next_user = 2;
                }

            }
        }
        else if (next_user == 2) {
            if (now_power > next_user_requirement) {
                now_power-=next_user_requirement;
                int now_weak_user = std::get<1>(sorted_user_pair[service_pair - 1]);
                power_allocation_matrix[now_weak_user] = next_user_requirement;

                int now_strong_user = std::get<2>(sorted_user_pair[service_pair - 1]);
                double strong_user_power = power_allocation_matrix[now_strong_user];
                double strong_user_channel_gain = Channel_gain_matrix[now_strong_user][0] + IRS_channel_gain_matrix[now_strong_user];
                double weak_user_minimum_rate = maximum_requirement_matrix[now_weak_user] * minimum_satisfaction_matrix[now_weak_user];
                //std::cout<<"check add weak:" <<calculate_weak_user_VLC_data_rate(strong_user_channel_gain, strong_user_power, next_user_requirement)<<" "<<weak_user_minimum_rate <<std::endl;


                next_user = 1;
                satisfied_pair++;
                if (satisfied_pair != now_pair_count) {
                    int next_strong_user = std::get<1>(sorted_user_pair[service_pair]);
                    next_user_requirement = calculate_strong_user_require_power(next_strong_user);
                }

            }
        }

        //all satisfied
        if (satisfied_pair == now_pair_count) {
            next_user = 0;
        }

    }



}

void random_exclude_user() {
    while(1) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distrib(0, algo2_sorted_ue_list.size() - 1);
        int random_number = distrib(gen);
        int delete_user = algo2_sorted_ue_list[random_number].second;
        if (excluded_user.count(delete_user) == 0) {
            excluded_user.insert(delete_user);
            //std::cout<<"delete user:"<<delete_user<<std::endl;
            break;
        }
    }
}

void exclude_largest_requirement_user() {
    double max_power = 0;
    int max_power_user = -1;
    for(int i = 0;i < algo2_match.size();i++) {
        int wu = algo2_sorted_ue_list[i].second;
        int su = algo2_sorted_ue_list[algo2_match[i] + (now_user / 2)].second;
        double strong_user_power = calculate_strong_user_require_power(su);
        //std::cout<<"strong user:"<<strong_user_power<<std::endl;
        long double weak_user_power = link_selection_matrix[i] == 0 ? calculate_weak_user_require_power(strong_user_power, wu) :
            calculate_weak_user_RF_require_power(strong_user_power, wu, su, relay_user);

        if (strong_user_power > max_power) {
            max_power = strong_user_power;
            max_power_user = su;
        }
        if (weak_user_power > max_power) {
            max_power = weak_user_power;
            max_power_user = wu;
        }
        //std::get<0>(sorted_user_pair[i]);

    }

    excluded_user.insert(max_power_user);
    //std::cout<<"delete user:"<<max_power_user<<std::endl;
}

void exclude_user() {
    //std::cout<<"delete user!"<<std::endl;
    //random
    if (exclude_user_method == 0) {
        random_exclude_user();
    }
    else if (exclude_user_method == 1) { // worst channel gain
        int delete_user = algo2_sorted_ue_list[0].second;
        excluded_user.insert(delete_user);
        //std::cout<<"delete user:"<<delete_user<<std::endl;
    }
    else { //largest requirement user
        exclude_largest_requirement_user();
    }
}

//use all residual power or not
bool binary_search_residual_power_allocation(int satisfied_user, double &residual_power, std::vector<std::tuple<double, int, int>> &sorted_user_pair) {
    double accuracy = 1e-10;
    bool max_satisfied = true;
    std::vector<bool> mark_user(now_user, false);
    int next_mark = -1;
    double accurate_satisfication = 0.0;
    double round_final_satisfaction = 0.0;
    for(int i = 0;i < now_user;i++) {
        calculate_user_satisfaction();
        std::vector<double> residual_power_matrix(UE, 0);
        double left = 0;
        double right = 1;
        for(int i = 0;i < satisfied_user;i++) {
            int wu = std::get<1>(sorted_user_pair[i]);
            int su = std::get<2>(sorted_user_pair[i]);
            if (mark_user[su] == false && right > (1.0 - user_satisfaction_matrix[su])) {
                right = 1.0 - user_satisfaction_matrix[su];
                next_mark = su;
            }
            if (wu != -1) {
                if (mark_user[wu] == false && right > (1.0 - user_satisfaction_matrix[wu])) {
                    right = 1.0 - user_satisfaction_matrix[wu];
                    next_mark = wu;
                }
            }
        }
        //std::cout<<"satisfaction:"<<right<<std::endl;
        while(right - left > accuracy) {
            double mid = (right + left) / 2.0;
            //std::cout<<"mid:"<<mid<<std::endl;
            double sum_power = 0.0;
            std::vector<double> tmp_power_matrix(UE, 0);
            for(int i = 0;i < satisfied_user;i++) {
                int wu = std::get<1>(sorted_user_pair[i]);
                int su = std::get<2>(sorted_user_pair[i]);
                double su_temp_satisfaction = minimum_satisfaction_matrix[su];
                if (mark_user[su] == false) {
                    minimum_satisfaction_matrix[su]+=(mid + accurate_satisfication);
                }
                else {
                    minimum_satisfaction_matrix[su] = 1.0;
                }

                double strong_user_require_power = calculate_strong_user_require_power(su);
                strong_user_require_power -= power_allocation_matrix[su];
                minimum_satisfaction_matrix[su] = su_temp_satisfaction;
                tmp_power_matrix[su] = strong_user_require_power;


                double weak_user_require_power = 0;

                //std::cout<<"wu:"<<wu<<"su:"<<su<<std::endl;
                if (wu != -1) {
                    double wu_temp_satisfaction = minimum_satisfaction_matrix[wu];
                    if (mark_user[wu] == false) {
                        minimum_satisfaction_matrix[wu]+=(mid + accurate_satisfication);
                    }
                    else {
                        minimum_satisfaction_matrix[wu] = 1.0;
                    }

                    double strong_user_power = strong_user_require_power + power_allocation_matrix[su];
                    weak_user_require_power = link_selection_matrix[wu] == 0 ? calculate_weak_user_require_power(strong_user_power, wu) :
                                    bsrpa_calculate_weak_user_RF_require_power(strong_user_power, wu, su, relay_user);
                    weak_user_require_power -= power_allocation_matrix[wu];
                    minimum_satisfaction_matrix[wu] = wu_temp_satisfaction;
                    tmp_power_matrix[wu] = weak_user_require_power;
                }

                sum_power+= (strong_user_require_power + weak_user_require_power);

            }
            //std::cout<<"residual power:"<<residual_power<<" sum power:"<<sum_power<<std::endl;
            if (sum_power > residual_power) {
                max_satisfied = false;
                right = mid;
            }
            else {
                round_final_satisfaction = mid;
                residual_power_matrix = tmp_power_matrix;
                left = mid;
            }
        }
        for(int i = 0;i < residual_power_matrix.size();i++) {
            power_allocation_matrix[i] += residual_power_matrix[i];
            residual_power-=residual_power_matrix[i];
        }

        if (max_satisfied == false) {
            return false;
        }
        else { // mark max satisfied user
            mark_user[next_mark] = true;
            accurate_satisfication+=round_final_satisfaction;
        }


    }

    return true;


}

void weak_user_average_residual_power_allocation(int satisfied_user, double residual_power, std::vector<std::tuple<double, int, int>> &sorted_user_pair) {

    for(int i = 0;i < satisfied_user;i++) {
        int wu = std::get<1>(sorted_user_pair[i]);
        if (wu != -1)
            power_allocation_matrix[wu]+=(residual_power / (satisfied_user - have_psuedo_user));
    }
}

void check_all_satisfied() {
    calculate_user_satisfaction();
    for(int i = 0;i < UE;i++) {
        if (excluded_user.count(i) == 0 && user_satisfaction_matrix[i] < minimum_satisfaction_matrix[i] * 0.99) {
            print_user_requirement();
            print_user_satisfaction();

        }
    }
}

void algorithm2() {
    while(1) {
        iteration_count++;
        clear_power_allocation_matrix();
        all_satisfied = false;
        algo2_list_user();

        for(int i = 0;i < UE;i++) {
            power_allocation_matrix[i] = 0.0;
            IRS_channel_gain_matrix[i] = 0.0;
            link_selection_matrix[i] = 0;

        }
        //print_user_requirement();
        calculate_user_satisfaction();
        //print_user_satisfaction();

        algo2_user_pairing();
        //std::cout<<"!"<<std::endl;
        //
        algo2_link_selection();



        double now_power = total_power;
        std::vector<std::tuple<double, int, int>> sorted_user_pair;
        int satisfied_pair = algo2_power_allocation(now_power, sorted_user_pair);

        calculate_user_satisfaction();

        //std::cout<<"IRS assignment start"<<std::endl;

        IRS_assignment(now_power,satisfied_pair, sorted_user_pair);

        //std::cout<<"IRS assignment over"<<std::endl;
        //print_user_requirement();
        calculate_user_satisfaction();
        //print_user_satisfaction();

        /*
        std::cout<<"all satisfied:" <<all_satisfied<<std::endl;
        std::cout<<"residual power: "<<now_power<<std::endl;
        */



        //average_residual_power_allocation(satisfied_pair, now_power, sorted_user_pair);

        if (all_satisfied == false) {
            exclude_user();
        }
        else {
            //check_all_satisfied();
            bool all_max_satisfied = binary_search_residual_power_allocation(satisfied_pair, now_power, sorted_user_pair);
            //std::cout<<"unused power:"<<now_power<<std::endl;
            if (now_power > 0.0)
                weak_user_average_residual_power_allocation(satisfied_pair, now_power, sorted_user_pair);
            break;

        }

    }
}
