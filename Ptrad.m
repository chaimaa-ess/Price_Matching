function  Ptrad = Ptrad( PrBS, PrS, b, s )
%Ptrad defines the price of trading: the trading price is the second best
%price proposed (if it exixsts) or the average value between the price
%proposed by the seller and the price offered by the buyer

    vect = PrBS(:,s);
    vectUp = vect;
    
    %%% keep only values lower than PrBS(b,s)
    for i = 1:length(vect)
        if vect(i) >= PrBS(b,s)
            ind = vectUp==vect(i);
            vectUp(ind) = [];
        end
    end
    vect = vectUp;
    if isempty(vect) % < PrS(s)
        Ptrad = (PrBS(b,s)+PrS(s))/2 ;
    else
        Ptrad = max(vect) ;
    end

end

