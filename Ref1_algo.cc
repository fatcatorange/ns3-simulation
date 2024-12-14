#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <unordered_set>


#include "Channel.h"
#include "global_environment.h"
#include "UE_AP_state.h"

#include "Ref1_algo.h"

using namespace std;

std::vector<std::pair<double, int>> sorted_ue_list; //list after sorted by channel gain
std::vector<double> weight_matrix(UE, 1);

/*
initially, strongest user pairing with the weakest user,
and every weak user using direct link,
*/

void init_ref1_algo() {
    for(int i = 0;i < UE;i++) {
        sorted_ue_list.push_back({Channel_gain_matrix[i][0], i});
    }

    std::sort(sorted_ue_list.begin(), sorted_ue_list.end());

    for(int i = 0;i < UE / 2;i++) { //strongest user pairing with weakest user
        pairing_matrix[sorted_ue_list[i].second][sorted_ue_list[UE - i -1].second] = 1;
        pairing_matrix[sorted_ue_list[UE - i -1].second][sorted_ue_list[i].second] = 1;
    }

}


double calculate_lambda(std::vector<std::pair<int, int>> &tmp_pairing) {
    long double left = 1e-6;
    long double right = 1e9;
    long double convergens = total_power * 0.01;
    long double lambda = 1;

    while(left <= right ) {
        std::cout<<left<<" "<<right<<std::endl;
        long double now_lambda = (left + right) / 2;
        long double now_power = 0;
        for (int i = 0;i < tmp_pairing.size();i++) {
            int now_strong_user = tmp_pairing[i].first;
            int now_weak_user = tmp_pairing[i].second;
            double term1 = 1;
            if (link_selection_matrix[now_weak_user] == 0) { //direct link
                term1 = weight_matrix[now_weak_user] * VLC_AP_Bandwidth / (2 * pairing_count * now_lambda);
            }
            else {
                term1 = weight_matrix[now_strong_user] * VLC_AP_Bandwidth / (2 * pairing_count * now_lambda);
            }
            double noise = VLC_AP_Bandwidth / (2 * pairing_count) * Nl * pow(10, 6);
            double c = 1 / (2 * M_PI);
            double signal = c;
            signal *= pow(VLC_optical_to_electric_factor, 2);
            signal *= pow(VLC_electric_to_optical_factor, 2);
            signal *= pow(Channel_gain_matrix[now_strong_user][0], 2);
            signal /= noise;

            double term2 = 1 / signal;

            now_power+=std::max(0.0, term1 - term2);
        }
        //std::cout<<total_power<<" "<<now_power<<std::endl;
        if (now_power <= total_power) {
            lambda = now_lambda;
            if (total_power - now_power <= convergens) {
                break;
            }
            right = now_lambda ;
        }
        else {
            left = now_lambda;
        }

    }

    return lambda;
}

void pair_power_allocation(int strong_user, int weak_user, double lambda) {
    double term1 = 1;
    if (link_selection_matrix[weak_user] == 0) { //direct link
        term1 = weight_matrix[weak_user] * VLC_AP_Bandwidth / (2 * pairing_count * lambda);
    }
    else {
        term1 = weight_matrix[strong_user] * VLC_AP_Bandwidth / (2 * pairing_count * lambda);
    }
    double noise = VLC_AP_Bandwidth / (2 * pairing_count) * Nl * pow(10, 6);
    double c = 1 / (2 * M_PI);
    double psi_j = c;
    psi_j *= pow(VLC_optical_to_electric_factor, 2);
    psi_j *= pow(VLC_electric_to_optical_factor, 2);
    psi_j *= pow(Channel_gain_matrix[strong_user][0], 2);
    psi_j /= noise;

    double term2 = 1 / psi_j;
    double pair_total_power = std::max(0.0, term1 - term2);
    //std::cout<<"pair power:"<<pair_total_power<<std::endl;
    link_selection_matrix[weak_user] = 1;

    if (link_selection_matrix[weak_user] == 1) { //case 1, formula (20)
        double ro1 = -1 + sqrt(1 + (pair_total_power * psi_j));
        ro1 /= psi_j;

        double A = 2 * calculate_RF_data_rate(Channel_gain_matrix[strong_user][0], UE_node_list[strong_user]->node, UE_node_list[weak_user]->node) * pairing_count / (VLC_AP_Bandwidth * 1e6);
        A = pow(2, A);

        double ro2 = pair_total_power * psi_j + 1 - A;
        //std::cout<<"ro:"<<ro1<<" "<<ro2<<std::endl;
        power_allocation_matrix[strong_user] = std::min(ro1, ro2);
        power_allocation_matrix[weak_user] = pair_total_power - power_allocation_matrix[strong_user];
    }



    //std::cout<<total_power<<std::endl;
}

void ref1_power_allocation(std::vector<std::pair<int, int>> &tmp_pairing) {
    double lambda = calculate_lambda(tmp_pairing);
    for(int i = 0;i < tmp_pairing.size();i++) {
        std::cout<<"st:"<<tmp_pairing[i].first<<" "<<tmp_pairing[i].second<<std::endl;
        pair_power_allocation(tmp_pairing[i].first, tmp_pairing[i].second, lambda);
    }
}

void check_pairing(std::vector<std::pair<int, int>> &tmp_pairing) {
    std::unordered_set<int> visited;
    for(int i = 0;i < UE;i++) {
        if (visited.count(i) == 1) {
            continue;
        }
        for(int j = 0;j < UE;j++) {
            if (pairing_matrix[i][j] == 1) {
                if (Channel_gain_matrix[i][0] > Channel_gain_matrix[j][0]) {
                    tmp_pairing.push_back({i, j});
                }
                else {
                    tmp_pairing.push_back({j, i});
                }
                visited.insert(i);
                visited.insert(j);
                break;
            }
        }
    }
}

void ref1_algo() {
    init_ref1_algo();
    int round = 0;

    //while(round <= maximum_iteration) {
        std::vector<std::pair<int, int>> tmp_pairing; // [strong user, weak user]
        check_pairing(tmp_pairing);
        ref1_power_allocation(tmp_pairing);
        //ref1_user_pairing();
        //ref1_link_selection();
        //clear_pairing_matrix();
        //update_weight();
        round++;

    //}


}
