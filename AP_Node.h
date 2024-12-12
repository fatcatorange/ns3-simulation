#ifndef AP_NODE_H_INCLUDED
#define AP_NODE_H_INCLUDED

#include <ctype.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

#include "ns3/core-module.h"
#include "ns3/applications-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
//#include "ns3/point-to-point-module.h"
#include "ns3/ipv4-global-routing-helper.h"
#include "ns3/vlc-channel-helper.h"
#include "ns3/vlc-device-helper.h"
#include "ns3/netanim-module.h"

using namespace ns3;

class AP_node{

public:
    int AP_ID; //the id < VLC_node are VLC_AP else are RF_AP
    Vector position;
    Ptr<Node> node;

    AP_node(int id,Vector pos,Ptr<Node> _node){
        AP_ID = id;
        position = pos;
        node = _node;
    }

};

#endif // AP_NODE_H_INCLUDED
