#ifndef REF1_ALGO_H_INCLUDED
#define REF1_ALGO_H_INCLUDED

#include "ns3/core-module.h"
#include "ns3/applications-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
//#include "ns3/point-to-point-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/vlc-channel-helper.h"
#include "ns3/vlc-device-helper.h"
#include "ns3/netanim-module.h"

/*
simulation for https://ieeexplore.ieee.org/document/9259258

*/

void init_ref1_algo();

void ref1_algo();

#endif // REF1_ALGO_H_INCLUDED
