%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Cooperative Game Theory                               %
%                Gale Shapley (Matching) Adapted version                %
%                Integration of the Pricing  Preferences                %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [l,lU,e, iter, Profit]  = matchingPrice(dist,distU,D,C,N)
l=0;   %losses in the networked µgrid
lU=0;  %losses due to exchange with the main grid 
iter=0;%nbr of iterations til convergence 

%The profit the participants make by exchanging energy internally
Profit=zeros(length(C),1);  

%Parameters of our distribution line system   
beta=0.02;                  % transformer loss rate
U0=50*10^3;                 % (kV) voltage of the utility energy 
U1=22*10^3;                 % (kV) voltage of the local electric network
R= 0.2;                     % (ohm/km) losses in the main grid
Pu=0.5;                     % Price of utiliy grid 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(C) > 1
    
%Separate buyers from sellers
[ buyers, sellers ] = split_buyers_sellers( D, C );
      
    DUp=D;               %Updated matrix of energy
    sellersUp=sellers;   %Updated list of sellers
    buyersUp=buyers;     %Updated list of buyers
    loss=zeros(N, N);    %matrix of losses during the exchange
    Exch=zeros(N, N);    %matrix of energy exchanged between µgrids
    PriceBSm=zeros(length(buyersUp), length(sellersUp)); %matrix of prices
    
    %%%    Execute the matching until the sum of energy to sell/buy =0 or
    %%%    the internal trade is no more beneficial
    while ~( isempty(sellersUp) || isempty(buyersUp) || all(all(isnan(PriceBSm))) )

        %%%%%%%%%%%%%%%%%%% Matching Game buyer_optimized %%%%%%%%%%%%%%%%%%%%
        [PriceBSm, PowerBSmax, Power0B, PriceS] = price_rank( buyersUp, sellersUp, DUp, distU, dist);
        
        %buyer_free=zeros(1,length(buyersUp));               %%%Line added for all-loop
        seller_partner=zeros(1,length(sellersUp));
        
        %condition=zeros(1,1);                               %%%Line added for all-loop
        %while min(condition)==0                             %%%Line added for all-loop

            for n=1:length(buyersUp)
                best_s_index=find(min(PriceBSm(n,:)));
               % PriceBSm(n,best_s_index) = 0;               %%%Line added for all-loop
                if isempty(best_s_index)
                    continue
                end
                if seller_partner(best_s_index)==0
                    seller_partner(best_s_index)=buyersUp(n);
                    %buyer_free(n)=1;                        %%%Line added for all-loop
                else 
                    % the buyer itself and not its index
                    comp=seller_partner(best_s_index); 
                    %index of competitor in buyersUp
                    competitor= buyersUp==comp; 

                    if PriceBSm(n,best_s_index) > PriceBSm(competitor,best_s_index)
                        seller_partner(best_s_index)=buyersUp(n);
                        %buyer_free(n)=1;                    %%%Line added for all-loop
                        %buyer_free(competitor)=0;           %%%Line added for all-loop
                    end
                end
            end
            %if length(buyersUp) <= length(sellersUp)        %%%Line added for all-loop
                %condition=buyer_free;                       %%%Line added for all-loop
            %else
                %condition=seller_partner;                   %%%Line added for all-loop
            %end
            
       %end
    %%%%%%%%%%%%%%%%%%%%%% Energy Exchange  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sellersUp2=sellersUp;
        buyersUp2=buyersUp;
        for s=1:length(sellersUp)  %the conflict in the indices btw Up & Up2
            b=seller_partner(s);
            if b~=0
                bId=buyersUp==b;
                E = PowerBSmax(bId,s);
                loss(b,sellersUp(s)) = E^2*R*dist(sellersUp(s),b)/(U1^2);
                DUp(sellersUp(s)) = DUp(sellersUp(s)) + E; 
                DUp(b)= DUp(b)- ( E - loss(b,sellersUp(s)) ); 

                if DUp(b)==0
                   bId2=buyersUp2==b;
                   buyersUp2(bId2)=[]; 
                end
                if DUp(sellersUp(s))==0
                   SId= sellersUp2 == sellersUp(s);
                   sellersUp2(SId)=[];
                end
    
                Exch(b,sellersUp(s))=E;
                Pt = Ptrad(PriceBSm, PriceS, bId, s);
                %E0b = E + E^2*distU(b)*R/(U0^2) + ( beta*E );
                Profit(b) = (Pu*Power0B(bId)) - (Pt*E);
                %Es0 = E - E^2*distU(sellersUp(s))*R/(U0^2) - ( beta*E );
                Profit(sellersUp(s)) = (Pt - PriceS(s)) * E;
            end
        end
        sellersUp=sellersUp2;
        buyersUp=buyersUp2;
        iter=iter+1;
    end       %end big while
  
    %losses of energy exchange between the subset and the utility
    for pf=1:length(C)          
        lU = lU + (DUp(C(pf))^2)*distU(C(pf))*R/(U0^2) + ( beta*abs(DUp(C(pf))) );       
    end
    %add internal losses
     l=sum(sum(loss));
    e = 1 - ( sum((abs(DUp))) / sum(abs(D)) ) ;
else 
    lU = ( D(C)^2*distU(C)*R/(U0^2) ) + ( beta*abs(D(C)) );
    e=0;
end