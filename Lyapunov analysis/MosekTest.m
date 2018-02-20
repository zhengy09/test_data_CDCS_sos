

clc;clear;close all;

warning on
load LyapunovSeDuMiData.mat

TimeT = zeros(length(n),1); 
Cost = zeros(length(n),1); 

%%

for i = 1:length(n)
    
    At = SeDuMiData{i}.At;
    b = SeDuMiData{i}.b;
    c = SeDuMiData{i}.c;
    K = SeDuMiData{i}.K;
    %% solutions using different method
    
    prob = convert_sedumi2mosek(At, b, c, K);
    
    try
        Tmosek = tic;
        [rcode, res] = mosekopt('minimize info', prob);
        TimeT(i,1) = toc(Tmosek);
    catch
        warning('out of memory')
    end

    save QuarticPolyMosek.mat
    
    
    %% statistics
    Cost(i) = [res.sol.itr.pobjval];
    
    save MosekTest
end