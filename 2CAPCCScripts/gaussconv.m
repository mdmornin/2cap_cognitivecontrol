function [fr,binned_data] = gaussconv(stMtx, pBin)
        isi = diff(stMtx);
        mean_isi = nanmean(isi);
        std_isi = nanstd(isi);
        binned_data = histc(stMtx,[min(stMtx(:)):pBin:max(stMtx(:))]);
        %Gaussian Parameters
        fr = zeros(length(binned_data),min(size(binned_data)));
        for i = 1:min(size(stMtx))
            if ~isnan(mean_isi(i))
                gauss_window = 1./pBin; % Constant
                gauss_SD(i) = mean_isi(i)^(1/2) * (1 / std_isi(i)/mean_isi(i));    %Half of Mean ISI

                gk = fspecial('gaussian',gauss_window,gauss_SD(i)); gk = gk./pBin;
                fr(:,i) = conv2(binned_data(:,i),gk,'same');
            end
        end
end