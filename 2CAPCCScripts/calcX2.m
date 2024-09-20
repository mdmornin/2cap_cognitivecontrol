Wdata = load('WProps.mat');
Pdata = load('PProps.mat');



for i = 1:15

    dataConNumP = floor(Pdata.congruent(i,:)*48);
    dataIncNumP = floor(Pdata.incongruent(i,:)*48);
    X = [dataConNumP(:),dataIncNumP(:)];
    
    [hp(i),pp(i),pXSquared(i)] = chi2cont(X);
end
    [hp,~,~,pp] = fdr_bh(pp);

for i = 1:15

    dataConNumW = floor(Wdata.congruent(i,:)*48);
    dataIncNumW = floor(Wdata.incongruent(i,:)*48);
    X = [dataConNumW(:),dataIncNumW(:)];
    
    [hw(i),pw(i),wXSquared(i)] = chi2cont(X);
end
    [hw,~,~,pw] = fdr_bh(pw);
