function [PriceBSm, PowerBSmax, Power0B, PriceS] = price_rank( buyersUp, sellersUp, DUp, distU, dist)
% PriceBSm(b,s): the price b proposes to s, =nan if the proposed price
% falls under the value the seller s is requiring
% PowerBSmax(b,s): max power that can be exchanged between b and s
% Power0B(b): the power that can be exchanged between b and utility
% PriceS: the min price the seller s is accepting from the buyers
% The price is defined following a logic that avoids additional losses
%with the grid

%Parameters of our distribution line system   
beta=0.02;                  % transformer loss rate
U0=50*10^3;                 % (kV) voltage of the utility energy 
U1=22*10^3;                 % (kV) voltage of the local electric network
R= 0.2; %*10^(-3);          % (ohm/km) losses in the main grid
Pu=0.5;

PriceBSm=zeros(length(buyersUp), length(sellersUp));
PriceS=zeros(1,length(sellersUp));
        Power0S=zeros(1,length(sellersUp));
        for s=1:length(sellersUp)
            Power0S(s)=abs(DUp(sellersUp(s))) - (DUp(sellersUp(s))^2)*distU(sellersUp(s))*R/(U0^2) - ( beta*abs(DUp(sellersUp(s))) );
            PriceS(s)=Pu*Power0S(s) /  abs(DUp(sellersUp(s)));
        end

        PriceBS=zeros(length(buyersUp), length(sellersUp));
        PowerBS=zeros(length(buyersUp), length(sellersUp));
        PowerBSmax=zeros(length(buyersUp), length(sellersUp));
        Power0B=zeros(1, length(buyersUp));
        
        for b=1:length(buyersUp)
                
            for s=1:length(sellersUp)
                Eqn=[(R*dist(sellersUp(s),buyersUp(b))/(U1^2)) -1 DUp(buyersUp(b))];
                r=roots(Eqn);
                    if isreal(r)
                        r(r<=0)=nan;
                        PowerBS(b,s)=min(r); 
                    else
                       PowerBS(b,s)=U1^2/(2*R*dist(sellersUp(s),buyersUp(b)));  
                    end
                PowerBSmax(b,s)=min(PowerBS(b,s),abs(DUp(sellersUp(s)))); 
                
                %%%%% Ps,min
%                 Power0S(s) = PowerBSmax(b,s) - (PowerBSmax(b,s)^2)*distU(sellersUp(s))*R/(U0^2) - ( beta*PowerBSmax(b,s) );
%                 PriceS(s)=Pu*Power0S(s) /  PowerBSmax(b,s);
                
                %%%%%% Psb 
                PowerBSrec = PowerBSmax(b,s) - ( (R*dist(sellersUp(s),buyersUp(b))/(U1^2))*PowerBSmax(b,s)^2 ) ;
                Eqn0=[distU(buyersUp(b))*R/(U0^2) (beta-1) PowerBSrec];
                r0=roots(Eqn0);
                if isreal(r0)
                    r0(r0<=0)=nan;
                    Power0B(b)=min(r0);
                else
                    Power0B(b) = U0^2/(2*R*distU(buyersUp(b)));
                end
                PriceBS(b,s)=Pu*Power0B(b) / PowerBSmax(b,s);
                if PriceBS(b,s) < PriceS(s)
                    PriceBSm(b,s) = nan;
                else
                    PriceBSm(b,s) = PriceBS(b,s);
                end
            end
        end
end

