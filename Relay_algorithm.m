%% proposed solution : Relay Algorithm (Ground BS or Satellite)

%% Define Parameters and Configurations
% Get from main.m
Bandwidth = 150;    % for NB-IoT, khz

% initial set relaying ground BS is true
GBS = true;

%% selection criteria
% distance to BS (If the UE device is far from the BS, choose SAT)
dis_to_BS = 50;   % UE to BS less than 100 km
if(dis_to_BS > 100)
    GBS = false;
end

% obstacle between UE & BS (mountain or buildings)
% If there are obstacles between UE & BS, choose SAT)
if(obstacle == true)
    GBS = false;
end

% weather conditions (rainy days choose Ground BS)
if(weather == "rainy")
    GBS = true;
end

% bandwidth constraints (GBS : 180kHz) 
if(weather == "rainy")
    GBS = true;
end

% Latency and Throughput Requirements
% If low latency or high throughput is required, ground-based transmission might be preferred.
% Satellite communication generally has higher latency due to the longer travel distance.


% Energy Consumption
% Satellite communication might require more power due to longer transmission distances.
% Consider the energy consumption requirements for UE devices, especially battery-operated ones.

%% Relay Selection Algorithm
% Implement algorithms for selecting suitable relays based on criteria 
% such as signal strength, link quality, proximity

% Define relay candidates (e.g., LEO satellites, base stations)
relay_candidates = 1;
% Define criteria for relay selection (e.g., signal strength, link quality)
criteria = 1;

% Relay Selection
selected_relay = relaySelection(relay_candidates, criteria);

%% MATLAB function code for Relay Algorithm Implementation

% decide to access ground networks or NTN(satellites) networks

function [selected_relay] = relaySelection(relay_candidates, criteria)
    
    % Example 1. : select the relay with the strongest signal strength
    % Example 2. : select the relay with the nearest networks
    % d > 100 km : satellites relay
    % d < 100 km : ground networks relay
    [~, idx] = max(criteria.signal_strength);
    selected_relay = relay_candidates(idx);
end

function [handover_decision] = performHandover(current_relay, next_relay)
    % Implement handover decision logic
    % Example: perform handover if the next relay provides better signal quality
    if next_relay.signal_quality > current_relay.signal_quality
        handover_decision = true;
    else
        handover_decision = false;
    end
end
