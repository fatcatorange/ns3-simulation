/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

//
// Network topology
//
//           10Mb/s, 10ms       10Mb/s, 10ms
//       n0-----------------n1-----------------n2
//
//
// - Tracing of queues and packet receptions to file
//   "tcp-large-transfer.tr"
// - pcap traces also generated in the following files
//   "tcp-large-transfer-$n-$i.pcap" where n and i represent node and interface
// numbers respectively
//  Usage (e.g.): ./waf --run tcp-large-transfer
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
#include "Channel.h"
#include "global_environment.h"
#include "UE_AP_state.h"
#include "Ref1_algo.h"

using namespace ns3;
using namespace std;

NS_LOG_COMPONENT_DEFINE ("VisibleLightCommunication");

void PrintPacketData(Ptr<const Packet> p, uint32_t size);

// The number of bytes to send in this simulation.
static const uint32_t totalTxBytes = 80000;
static uint32_t currentTxBytes = 0;

// Perform series of 1040 byte writes (this is a multiple of 26 since
// we want to detect data splicing in the output stream)
static const uint32_t writeSize = 2048;
uint8_t data[writeSize];

// These are for starting the writing process, and handling the sending
// socket's notification upcalls (events).  These two together more or less
// implement a sending "Application", although not a proper ns3::Application
// subclass.

void StartFlow(Ptr<Socket>, Ipv4Address, uint16_t);

void WriteUntilBufferFull(Ptr<Socket>, uint32_t);

std::vector<double> Received(1, 0);
std::vector<double> theTime(1, 0);
//////////////////////////////////////
//Function to generate signals.
std::vector<double>& GenerateSignal(int size, double dutyRatio);

static void RxEnd(Ptr<const Packet> p) { // used for tracing and calculating throughput

	//PrintPacketData(p,p->GetSize());

	Received.push_back(Received.back() + p->GetSize()); // appends on the received packet to the received data up until that packet and adds that total to the end of the vector
	theTime.push_back(Simulator::Now().GetSeconds()); // keeps track of the time during simulation that a packet is received
	//NS_LOG_UNCOND("helooooooooooooooooo RxEnd");
}

static void TxEnd(Ptr<const Packet> p) { // also used as a trace and for calculating throughput

	Received.push_back(Received.back() + p->GetSize()); // same as for the RxEnd trace
	theTime.push_back(Simulator::Now().GetSeconds()); 	//
	//NS_LOG_UNCOND("helooooooooooooooooo TxEnd");
}

static void CwndTracer(uint32_t oldval, uint32_t newval) {
	NS_LOG_INFO("Moving cwnd from " << oldval << " to " << newval);
}

int main(int argc, char *argv[]) {
	// Users may find it convenient to turn on explicit debugging
	// for selected modules; the below lines suggest how to do this
	//  LogComponentEnable("TcpSocketImpl", LOG_LEVEL_ALL);
	//  LogComponentEnable("PacketSink", LOG_LEVEL_ALL);            // uncomment in original example
	//  LogComponentEnable("TcpLargeTransfer", LOG_LEVEL_ALL);

	//parameters:


	Simulator::Schedule(Seconds(0.0), &basic_init);
	Simulator::Schedule(Seconds(0.0), &calculate_throughput);
	//Simulator::Schedule(Seconds(0.0), &ref1_algo);

	Simulator::Stop(Seconds(simulation_time));//??
    Simulator::Run();

    throughput_write_file();

    //pair_fairness_write_file();



}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//begin implementation of sending "Application"
void StartFlow(Ptr<Socket> localSocket, Ipv4Address servAddress,
		uint16_t servPort) {
	//NS_LOG_UNCOND("helooooooooooooooooo StartFlow");
	localSocket->Connect(InetSocketAddress(servAddress, servPort)); //connect

	// tell the tcp implementation to call WriteUntilBufferFull again
	// if we blocked and new tx buffer space becomes available
	localSocket->SetSendCallback(MakeCallback(&WriteUntilBufferFull));
	WriteUntilBufferFull(localSocket, localSocket->GetTxAvailable());
}

void WriteUntilBufferFull(Ptr<Socket> localSocket, uint32_t txSpace) {
	//NS_LOG_UNCOND("helooooooooooooooooo WriteUntilBufferFull");
	while (currentTxBytes < totalTxBytes && localSocket->GetTxAvailable() > 0) {

		uint32_t left = totalTxBytes - currentTxBytes;
		uint32_t dataOffset = currentTxBytes % writeSize;
		uint32_t toWrite = writeSize - dataOffset;
		toWrite = std::min (toWrite, left);
		toWrite = std::min (toWrite, localSocket->GetTxAvailable ());

		Ptr<Packet> p = Create<Packet>(&data[dataOffset], toWrite);
		Ptr<Node> startingNode = localSocket->GetNode();
		Ptr<VlcTxNetDevice> txOne = DynamicCast<VlcTxNetDevice>(startingNode->GetDevice(0) );
		txOne->EnqueueDataPacket(p);

		int amountSent = localSocket->Send (&data[dataOffset], toWrite, 0);
		if(amountSent < 0)
		{
			// we will be called again when new tx space becomes available.
			return;
		}

		currentTxBytes += amountSent;
	}

	localSocket->Close();
}
/*
std::vector<double>& GenerateSignal(int size, double dutyRatio) {
	std::vector<double> *result = new std::vector<double>();
	result->reserve(size);

	double bias = 0;
	double Vmax = 4.5;
	double Vmin = 0.5;

	for (int i = 0; i < size; i++) {
		if (i < size * dutyRatio) {
			result->push_back(Vmax + bias);
		} else {
			result->push_back(Vmin + bias);
		}
	}

	return *result;
}

void PrintPacketData(Ptr<const Packet> p, uint32_t size) {
	uint8_t *data = new uint8_t[size];

	p->CopyData(data, size);

	for (uint32_t i = 0; i < size; i++) {
		std::cout << (int) data[i] << " ";
	}

	std::cout << std::endl;

	delete[] data;

}
*/
