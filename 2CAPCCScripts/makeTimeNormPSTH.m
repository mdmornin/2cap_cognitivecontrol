function [timeWarpedRateVector] = makeTimeNormPSTH(stMtx, Events, epochBins)


        epoch1_bintimes = linspace(Events(1),Events(2),epochBins(1));
        epoch2_bintimes = linspace(Events(2),Events(3),epochBins(2));
        if Events(3) == Events(4)
            epoch3_bintimes = linspace(Events(3),Events(4)+2,epochBins(3));
            epoch4_bintimes = linspace(Events(4)+2,Events(4)+5,epochBins(4));
        elseif Events(3) ~= Events(4)
            epoch3_bintimes = linspace(Events(3),Events(4),epochBins(3));
            epoch4_bintimes = linspace(Events(4),Events(4)+5,epochBins(4));
        end

        bin_width1 = diff(epoch1_bintimes);
        bin_width2 = diff(epoch2_bintimes);
        bin_width3 = diff(epoch3_bintimes);
        bin_width4 = diff(epoch4_bintimes);

        binned_data_epoch1 = histc(stMtx,[linspace(Events(1),Events(2),epochBins(1))]);
        binned_data_epoch2 = histc(stMtx,[linspace(Events(2),Events(3),epochBins(2))]);
        if Events(3) == Events(4)
            binned_data_epoch3 = histc(stMtx,[linspace(Events(3),Events(4)+2,epochBins(3))]);
            binned_data_epoch4 = histc(stMtx,[linspace(Events(4)+2,Events(4)+5,epochBins(4))]);
        elseif Events(3) ~= Events(4)
            binned_data_epoch3 = histc(stMtx,[linspace(Events(3),Events(4),epochBins(3))]);
            binned_data_epoch4 = histc(stMtx,[linspace(Events(4),Events(4)+5,epochBins(4))]);
        end
        

        
        % When "epochBins" are even, the difference between each bin-width
        % is the same. This might not be the case for odds, and need to
        % implement a way to check this. 
        
        epoch1_rates = binned_data_epoch1./bin_width1(1);
        epoch2_rates = binned_data_epoch2./bin_width2(1);
        epoch3_rates = binned_data_epoch3./bin_width3(1);
        epoch4_rates = binned_data_epoch4./bin_width4(1);

        
        timeWarpedRateVector = [epoch1_rates(1:end-1,:); epoch2_rates(1:end-1,:); epoch3_rates(1:end-1,:); epoch4_rates(1:end-1,:)];
        
        
    
end

