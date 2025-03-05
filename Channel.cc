#include <iostream>
#include <fstream>
#include <string>
#include <chrono> //seed
#include <cmath>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/mobility-module.h"
#include "Channel.h"
#include "global_environment.h"
#include <boost/math/distributions/rayleigh.hpp>

#define e 2.7182818284

std::normal_distribution<double> Gaussian (0.0,10);    //normal distribution 即 Gaussian distribution
boost::math::rayleigh_distribution<double> rayleigh(0.8); // standard rayleigh distribution
std::uniform_real_distribution<double> random_p(0.0, 1.0);// uniform random variable between 0.0 and 1.0 for inverse transform sampling
std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
int random_count = 100000;
double X=0.0,H=0.0,p;

double RadtoDeg(double radian)
{
    return radian * 180 / M_PI;
}

double DegtoRad(double degree){
    return degree * M_PI / 180;
}

double calculateDistance(Ptr<Node> AP,Ptr<Node> UE)
{
    Ptr<MobilityModel> AP_mobility = AP->GetObject<MobilityModel>();
    Ptr<MobilityModel> UE_mobility = UE->GetObject<MobilityModel>();

    Vector AP_pos = AP_mobility->GetPosition();
    Vector UE_pos = UE_mobility->GetPosition();

    return AP_mobility->GetDistanceFrom(UE_mobility);
}


double Estimate_one_VLC_Channel_Gain(Ptr<Node> VLC_AP,Ptr<Node> UE){

    double incidence_angle = Get_Incidence_Angle_AP_UE(VLC_AP,UE);

    if(RadtoDeg(incidence_angle) >= VLC_field_of_view)
    {
        return 0;
    }

    double lambertian_coefficient = (-1) / (log(cos(DegtoRad(VLC_PHI_half))));

    double irradiance_angle = incidence_angle;

    Ptr<MobilityModel> VLC_AP_MobilityModel = VLC_AP->GetObject<MobilityModel>();
    Vector VLC_AP_Pos = VLC_AP_MobilityModel->GetPosition();

    Ptr<MobilityModel> UE_MobilityModel = UE->GetObject<MobilityModel>();
    Vector UE_Pos = UE_MobilityModel->GetPosition();

    double height_diff = VLC_AP_Pos.z - UE_Pos.z;

    double distance = calculateDistance(VLC_AP,UE);

    double channel_gain = VLC_receiver_area * (lambertian_coefficient + 1)/(2*M_PI*pow(distance,2));

    channel_gain = channel_gain * pow(cos(irradiance_angle),lambertian_coefficient);

    channel_gain = channel_gain * VLC_filter_gain;

    channel_gain = channel_gain * pow(VLC_refractive_index , 2) / pow(sin(DegtoRad(VLC_field_of_view)) , 2);

    channel_gain = channel_gain * cos(incidence_angle);

    return channel_gain;


}

double Get_Incidence_Angle_AP_UE(Ptr<Node> AP,Ptr<Node> UE){

    Ptr<MobilityModel> AP_MobilityModel = AP->GetObject<MobilityModel>();
    Vector AP_Pos = AP_MobilityModel->GetPosition();

    Ptr<MobilityModel> UE_MobilityModel = UE->GetObject<MobilityModel>();
    Vector UE_Pos = UE_MobilityModel->GetPosition();

    double height_diff = AP_Pos.z - UE_Pos.z;

    double dx = AP_Pos.x - UE_Pos.x;
    double dy = AP_Pos.y - UE_Pos.y;

    double plane_diff = sqrt(pow(dx,2) + pow(dy,2));

    double hypoteuse = sqrt(pow(height_diff,2) + pow(plane_diff,2));

    const double angle = acos((pow(height_diff,2) + pow(hypoteuse,2) - pow(plane_diff,2)) / (2*height_diff*hypoteuse));

    return angle;
}

double calculate_strong_user_data_rate(double channel_gain ,double strong_user_power) {
    //std::cout<<"rate:"<<channel_gain<<" "<<strong_user_power<<std::endl;
    double data_rate = VLC_AP_Bandwidth / (2.0 * pairing_count);
    double noise = VLC_AP_Bandwidth / (pairing_count) * Nl * pow(10, 6); // A^2/Hz -> A^2/MHz
    double c = 1 / (2.0 * M_PI);
    double signal = c;


    signal *= pow(VLC_optical_to_electric_factor, 2.0);
    signal *= pow(VLC_electric_to_optical_factor, 2.0);
    signal *= pow(channel_gain, 2.0);
    signal *= strong_user_power;
    double tmp = signal / noise;
    data_rate*= log2(1 + (signal / noise));

    return data_rate;
}

double calculate_weak_user_VLC_data_rate(long double channel_gain ,long double strong_user_power, long double weak_user_power) {
    long double data_rate = VLC_AP_Bandwidth / (2.0 * (long double)pairing_count);
    long double noise = VLC_AP_Bandwidth / (long double)(pairing_count) * Nl * 1e6; // A^2/Hz -> A^2/MHz
    long double c = 1.0 / (2.0 * M_PI);
    long double signal = c;
    long double interference = c;

    interference *= pow(VLC_optical_to_electric_factor, 2.0);
    interference *= pow(VLC_electric_to_optical_factor, 2.0);
    interference *= pow(channel_gain, 2.0);
    interference *= strong_user_power;


    signal *= pow(VLC_optical_to_electric_factor, 2.0);
    signal *= pow(VLC_electric_to_optical_factor, 2.0);
    signal *= pow(channel_gain, 2.0);
    signal *= weak_user_power;


    data_rate*= log2(1.0 + (signal / (noise + interference)));


    return data_rate;
}

std::complex<double> generateComplexGaussian(double mu, double sigma) {
    // define random generator
    std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());
    // define normal distribution
    std::normal_distribution<double> distribution(mu, sigma);

    // 生成實部和虛部的高斯隨機變量
    double realPart = distribution(re);
    double imagPart = distribution(re);

    return std::complex<double>(realPart, imagPart);
}

double get_RF_channel_gain(Ptr<Node> user1, Ptr<Node> user2){
    std::complex<double> rf_channel(0,0);
    double K = 0;
    //std::uniform_real_distribution<double> unif(0, 1);

    //(realPart, imaginaryPart)
    std::complex<double> x1 = generateComplexGaussian(0.0, 1.0); //X1 is a complex Gaussian random variable with zero mean and unit variance
    double c = cos((angle_of_LOS_arrival * M_PI) / 180);
    double s = sin((angle_of_LOS_arrival * M_PI) / 180);
    std::complex<double> euler(c,s); //e^(jθ) = cos(θ) + j*sin(θ)
    double distance =calculateDistance(user1, user2);
    if(distance < breakpoint) K = 1;
    rf_channel = sqrt(K/(K+1));
    rf_channel *= euler;
    rf_channel += std::sqrt(1/(K+1)) * x1;
    double channel = std::abs(rf_channel) * abs(rf_channel);

    double pass_loss;
    double L_FS = 20*log10(distance);
    L_FS += 20*log10(central_carrier_frequency);
    L_FS -= 147.5;
    double sigma = 5;
    if(distance <= breakpoint){
        sigma = 3;
    }

    // define random generator
    std::default_random_engine re(std::chrono::system_clock::now().time_since_epoch().count());
    // define normal distribution
    std::normal_distribution<double> distribution(0, sigma);

    double X_sigma = distribution(re);
    if(distance <= breakpoint){
        pass_loss = L_FS + X_sigma;
    }
    else{
        pass_loss = L_FS + 35 * log(distance / breakpoint) + X_sigma;
    }

    double channel_gain = channel;
    channel_gain *= pow(10,(-pass_loss/10));
    //std::cout<<"RF channel gain: "<<channel_gain<<std::endl;
    return channel_gain;
}

double get_energy_harvest(double channel_gain){
    double dc_bias = (maximum_current + minimum_current) / 2;
    double short_circuit_DC_current = VLC_electric_to_optical_factor * VLC_optical_to_electric_factor;
    short_circuit_DC_current *= channel_gain;
    short_circuit_DC_current *= dc_bias;

    double voltage = short_circuit_DC_current / dark_saturation_current;
    voltage = log(1+voltage) / log(e);
    voltage *= thermal_voltage;
    double energy = fill_factor;
    energy *= short_circuit_DC_current * voltage;

    return energy;
}

double calculate_RF_data_rate(double channel_gain, Ptr<Node> user1, Ptr<Node> user2) {
    double RF_power = 0;
    RF_power = get_energy_harvest(channel_gain);
    double channel = get_RF_channel_gain(user1, user2);
    double SNR = RF_power * channel;
    SNR = SNR / (RF_AP_Bandwidth * Nw);
    double achievable_rate;
    achievable_rate = log(1+SNR);
    achievable_rate *= 0.5;
    achievable_rate *= RF_AP_Bandwidth;
    achievable_rate /= std::max(relay_user, 1);

    //std::cout<<"data rate:"<<achievable_rate / 1e6<<std::endl;
    return achievable_rate / 1e6; //bps to mbps
}

