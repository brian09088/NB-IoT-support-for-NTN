%% Satellite access  and handover strategy

%% Algorithm for decision transmit or relay routing path use ISL (Inter-Satellite Links)

%% Access and handover strategy

function [] = access_handover()

    if new UE arrive:
        determine the current geographic coordinates, 
        send location information and access request to GS

        GS sends : satellite coverage matrix to UE

        for i in N : 
            if t(s,i) < t(now) < t(theta,i)
                popup the element i in S
            end
            for m in S:
                select the satellite which has
            end
        end
    end
if X service-cycle is ended:
    
    pop the satellite j which can service to get the Q(ij,t) for each j in [S'' ...
        return j which satellite has the maximum Q(ij,t)
    
    if j != x :
        handover to satellite j
    else:
        service continue
    end

end        