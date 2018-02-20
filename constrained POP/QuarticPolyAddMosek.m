
%% Minimize the quartic polynomial over the hypercube
% f(x) = \sum_{1<i<j<n} xixj + xi^2xj-xj^3-xi^2xj^2
% gi(x) = \sum xi^2<1

clear;
load QuarticPoly

TimeT = zeros(length(N),8); % mosek, sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos

TimeT(:,2:8) = TimeTotal;
%TimeMosek = zeros(length(N))

for Index = 1:length(N)
    
    SeDuMiData = data{Index};
    
    % coverted into mosek format
    prob = convert_sedumi2mosek(SeDuMiData.At, SeDuMiData.b, SeDuMiData.c, SeDuMiData.K);
    
    try
        Tmosek = tic;
        [rcode, res] = mosekopt('minimize info', prob);
        TimeT(Index,1) = toc(Tmosek);
    catch
        warning('out of memory')
    end

    save QuarticPolyMosek.mat
end


