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
std::vector<int> match(UE / 2, -1);
const double INF = std::numeric_limits<double>::max();

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
    /*
    The paper does not explain what happens if condition
are not satisfied. Therefore, in such cases, I use Case 1 to allocate power.
    */
    //if (weak_user % 2 == 0)
    //link_selection_matrix[weak_user] = 1;
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

    double psi_i = c;
    psi_i *= pow(VLC_optical_to_electric_factor, 2);
    psi_i *= pow(VLC_electric_to_optical_factor, 2);
    psi_i *= pow(Channel_gain_matrix[weak_user][0], 2);
    psi_i /= noise;

    double term2 = 1 / psi_j;
    double pair_total_power = std::max(0.0, term1 - term2);
    //std::cout<<"pair power:"<<pair_total_power<<std::endl;

    //std::cout<<"psi: "<<weight_matrix[weak_user] <<" "<< weight_matrix[strong_user] <<" "<< (psi_j / psi_i)<<std::endl;
    if (link_selection_matrix[weak_user] == 0 && (weight_matrix[weak_user] / weight_matrix[strong_user]) < (psi_j / psi_i)) { //case 2, formula (26)
        double epsilon = 1e-15;
        if (abs(weight_matrix[strong_user] - weight_matrix[weak_user]) < epsilon) {
            weight_matrix[weak_user]+=epsilon;
        }
        double sigma = (weight_matrix[weak_user] * psi_i) - (weight_matrix[strong_user] * psi_j);
        sigma /= psi_j * psi_i * (weight_matrix[strong_user] - weight_matrix[weak_user]);
        std::cout<<"sigma:"<<sigma<<std::endl;
        if (pair_total_power > 2 * sigma) {
            power_allocation_matrix[strong_user] = sigma;
            power_allocation_matrix[weak_user] = pair_total_power - power_allocation_matrix[strong_user];
        }
        else {
            double ro1 = -1 + sqrt(1 + (pair_total_power * psi_j));
            ro1 /= psi_j;

            double A = 2 * calculate_RF_data_rate(Channel_gain_matrix[strong_user][0], UE_node_list[strong_user]->node, UE_node_list[weak_user]->node) * pairing_count / (VLC_AP_Bandwidth * 1e6);
            A = pow(2, A);

            double ro2 = pair_total_power * psi_j + 1 - A;
            //std::cout<<"ro:"<<ro1<<" "<<ro2<<std::endl;
            power_allocation_matrix[strong_user] = std::min(ro1, ro2);
            power_allocation_matrix[weak_user] = pair_total_power - power_allocation_matrix[strong_user];
        }


    }
    else { //case 2, formula (20)
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
        //std::cout<<"st:"<<tmp_pairing[i].first<<" "<<tmp_pairing[i].second<<std::endl;
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

/*
hungarian alogorithm
*/

double calculate_pair_cost(int user1, int user2) { //weak , strong
    double total_rate = 0;

    total_rate = calculate_strong_user_data_rate(Channel_gain_matrix[user2][0], power_allocation_matrix[user2]);
    if (link_selection_matrix[user1] == 1) {
        total_rate += weight_matrix[user1] * std::min(calculate_RF_data_rate(Channel_gain_matrix[user2][0], UE_node_list[user2]->node, UE_node_list[user1]->node),
                                    calculate_weak_user_VLC_data_rate(Channel_gain_matrix[user2][0], power_allocation_matrix[user2], power_allocation_matrix[user1]));
    }
    else {
        total_rate += weight_matrix[user2] * calculate_weak_user_VLC_data_rate(Channel_gain_matrix[user1][0], power_allocation_matrix[user2], power_allocation_matrix[user1]);
    }
    //std::cout<<user1<<" "<<user2<<" "<<total_rate<<std::endl;
    return total_rate;
}

void print_cost_matrix(std::vector<std::vector<double>> &cost_matrix) {
    for(int i = 0; i < UE / 2;i++) {
        for (int j = 0;j < UE / 2;j++) {
            std::cout<<" "<<cost_matrix[i][j];
        }
        std::cout<<std::endl;
    }
}

void calculate_cost_matrix(std::vector<std::vector<double>> &cost_matrix) {
    double maxi = 0;
    for(int i = 0; i < UE / 2;i++) {
        for (int j = 0;j < UE / 2;j++) {
            cost_matrix[i][j] = calculate_pair_cost(sorted_ue_list[i].second, sorted_ue_list[(UE / 2) + j].second); // weak, strong
            //std::cout<<cost_matrix[i][j]<<std::endl;
            maxi = std::max(maxi, cost_matrix[i][j]);
        }
    }
    print_cost_matrix(cost_matrix);
    //transfer

    for(int i = 0; i < UE / 2;i++) {
        for (int j = 0;j < UE / 2;j++) {
            cost_matrix[i][j] = maxi - cost_matrix[i][j];
        }
    }
}





// Helper function for DFS
bool dfs(int x, std::vector<double> &lx, std::vector<double> &ly, std::vector<int> &match, std::vector<bool> &s, std::vector<bool> &t, const std::vector<std::vector<double>> &cost) {
    int n = cost.size();
    s[x] = true;
    for (int y = 0; y < n; ++y) {
        if (t[y]) continue;
        double slack = lx[x] + ly[y] - cost[x][y];
        if (slack == 0) {
            t[y] = true;
            if (match[y] == -1 || dfs(match[y], lx, ly, match, s, t, cost)) {
                match[y] = x;
                return true;
            }
        }
    }
    return false;
}

// Update labels to reduce slack
void update_labels(std::vector<double> &lx, std::vector<double> &ly, const std::vector<bool> &s, const std::vector<bool> &t, const std::vector<std::vector<double>> &cost) {
    int n = cost.size();
    double delta = INF;
    for (int x = 0; x < n; ++x) {
        if (s[x]) {
            for (int y = 0; y < n; ++y) {
                if (!t[y]) {
                    delta = std::min(delta, lx[x] + ly[y] - cost[x][y]);
                }
            }
        }
    }
    for (int x = 0; x < n; ++x) {
        if (s[x]) lx[x] -= delta;
    }
    for (int y = 0; y < n; ++y) {
        if (t[y]) ly[y] += delta;
    }
}

// Hungarian algorithm implementation
std::vector<int> hungarian(const std::vector<std::vector<double>> &cost) {
    int n = cost.size();
    std::vector<double> lx(n, 0), ly(n, 0);
    std::vector<int> match(n, -1);

    // Initialize label lx to the maximum value in each row
    for (int i = 0; i < n; ++i) {
        lx[i] = *std::max_element(cost[i].begin(), cost[i].end());
    }

    // Find a perfect matching
    for (int i = 0; i < n; ++i) {
        while (true) {
            std::vector<bool> s(n, false), t(n, false);
            if (dfs(i, lx, ly, match, s, t, cost)) break;
            update_labels(lx, ly, s, t, cost);
        }
    }

    return match;
}

void ref1_user_pairing() {
    for (int i = 0;i < UE;i++) {
        for (int j = 0;j < UE;j++) {
            pairing_matrix[i][j] = 0;
        }
    }
    pairing_matrix.resize(UE, std::vector<int>(UE, 0));
    std::vector<std::vector<double>> cost_matrix(UE / 2, std::vector<double>(UE / 2, 0));
    calculate_cost_matrix(cost_matrix);
    print_cost_matrix(cost_matrix);
    match = hungarian(cost_matrix);
    for (int i = 0;i < match.size();i++) {
        //std::cout<<"weak user "<<sorted_ue_list[i].second<<" pairing with strong user"<<sorted_ue_list[i + (UE / 2)].second<<std::endl;
        pairing_matrix[sorted_ue_list[i].second][sorted_ue_list[i + (UE / 2)].second] = 1;
    }



}

double calculate_hybrid_link_bonus(int strong_user, int weak_user) {

    double VLC_rate = std::min(calculate_weak_user_VLC_data_rate(Channel_gain_matrix[strong_user][0], power_allocation_matrix[strong_user], power_allocation_matrix[weak_user]),
                                     calculate_RF_data_rate(Channel_gain_matrix[strong_user][0], UE_node_list[strong_user]->node, UE_node_list[weak_user]->node));
    double dir_rate = calculate_weak_user_VLC_data_rate(Channel_gain_matrix[weak_user][0], power_allocation_matrix[strong_user], power_allocation_matrix[weak_user]);
    return VLC_rate - dir_rate;
}

double calculate_link_table_pair_rate(int strong_user, int weak_user, int pair_link) {
    double rate = calculate_strong_user_data_rate(Channel_gain_matrix[strong_user][0], power_allocation_matrix[strong_user]);
    if (pair_link == 0) {
        rate+=std::min(calculate_weak_user_VLC_data_rate(Channel_gain_matrix[strong_user][0], power_allocation_matrix[strong_user], power_allocation_matrix[weak_user]),
                        calculate_RF_data_rate(Channel_gain_matrix[strong_user][0], UE_node_list[strong_user]->node, UE_node_list[weak_user]->node));
    }
    else {
        rate+=calculate_weak_user_VLC_data_rate(Channel_gain_matrix[weak_user][0], power_allocation_matrix[strong_user], power_allocation_matrix[weak_user]);
    }

    std::cout<<"dir:"<<calculate_weak_user_VLC_data_rate(Channel_gain_matrix[weak_user][0], power_allocation_matrix[strong_user], power_allocation_matrix[weak_user])<<std::endl;
    std::cout<<"hybrid:"<<std::min(calculate_weak_user_VLC_data_rate(Channel_gain_matrix[strong_user][0], power_allocation_matrix[strong_user], power_allocation_matrix[weak_user]),
                        calculate_RF_data_rate(Channel_gain_matrix[strong_user][0], UE_node_list[strong_user]->node, UE_node_list[weak_user]->node))<<std::endl;

    return rate;
}

void ref1_link_selection(std::vector<std::pair<int, int>> &tmp_pairing) {
    std::vector<std::vector<int>> link_selection_table(UE / 2 + 1, std::vector<int>(UE / 2, 0));
    for(int hybrid_user = 1;hybrid_user <= tmp_pairing.size() ;hybrid_user++) {
        std::priority_queue<pair<double, int>> rate_pq;
        for(int j = 0;j < tmp_pairing.size(); j++) {
            relay_user = hybrid_user;
            rate_pq.push({calculate_hybrid_link_bonus(tmp_pairing[j].first, tmp_pairing[j].second), j});
        }
        for(int i = 0;i < hybrid_user;i++) {
            int user = rate_pq.top().second;
            rate_pq.pop();
            link_selection_table[hybrid_user][user] = 1;
        }
    }
    double maximum_rate = 0;
    int maximum_index = 0;
    for(int i = 0;i < link_selection_table.size();i++) {
        double now_rate = 0;
        for(int j = 0;j < link_selection_table[i].size();j++) {
            now_rate+=calculate_link_table_pair_rate(tmp_pairing[j].first, tmp_pairing[j].second, link_selection_table[i][j]);
        }
        std::cout<<"now rate: "<<now_rate<<std::endl;
        if (now_rate > maximum_rate) {
            maximum_rate = now_rate;
            maximum_index = i;
        }

    }
    std::cout<<maximum_index<<std::endl;
    relay_user = maximum_index;
    link_selection_matrix.resize(UE, 0);
    for (int i = 0;i < link_selection_table[maximum_index].size();i++) {
        if (link_selection_table[maximum_index][i] == 1) {
            link_selection_matrix[tmp_pairing[i].second] = 1;
        }
    }

}



void ref1_algo() {
    init_ref1_algo();
    int round = 0;
    double last = 0;
    double now_rate = 0;

    while(round <= maximum_iteration) {
        std::vector<std::pair<int, int>> tmp_pairing; // [strong user, weak user]
        check_pairing(tmp_pairing);
        ref1_power_allocation(tmp_pairing);
        ref1_user_pairing();
        ref1_link_selection(tmp_pairing);
        //clear_pairing_matrix();
        //update_weight();
        calculate_data_rate_matrix();
        last = now_rate;
        now_rate = calculate_sum_rate();
        std::cout<<now_rate<<" "<<last<<std::endl;
        if (now_rate - last <= (now_rate / 1000)) {
            break;
        }
        round++;

    }


}
