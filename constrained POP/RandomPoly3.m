
%% Minimize the random polynomial
% min f(x) 
% gi(x) = \sum xi^2<1
clc;clear

N = round(logspace(log10(10),log10(30),10));
% N = [5,10,15,30];% round(logspace(log10(10),log10(40),10));
d = 2;         % half degree of the polynomial

TimePro = zeros(length(N),1);     % time for problem generation
TimeTotal = zeros(length(N),7);   % sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos
TimeSetup = zeros(length(N),3);   % scs-direct, scs-indirect, cdcs-sos
TimeADMM = zeros(length(N),3);    % scs-direct, scs-indirect, cdcs-sos
TimeAver = zeros(length(N),3);    % scs-direct, scs-indirect, cdcs-sos
Cost = zeros(length(N),7);        % sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos
Iter = zeros(length(N),6);        % sedumi, sdpt3, sdpa, scs-direct, scs-indirect, cdcs-sos
Density = zeros(length(N),4);

data = cell(length(N),1);

TimeTotalV = zeros(length(N),7); 
TimeADMMV  = zeros(length(N),3);
TimeAverV  = zeros(length(N),3);  
TimeSetupV = zeros(length(N),3); 
%%
Maxiter = 2e3;
Tol     = 1e-3;
TolIPM  = 1e-8;

load RandomPolyN5.mat

for Index = 7:length(N)

    n = N(Index);
    
    Nt = 5;
    
    TimeTotalj = zeros(Nt,7);   % sedumi, cdcs(primal), scs(direct)
    TimeSetupj = zeros(Nt,3);   % sosadmm, cdcs(primal), scs(direct)
    TimeADMMj = zeros(Nt,3);
    TimeAverj = zeros(Nt,3);
    Costj = zeros(Nt,7);
    Iterj = zeros(Nt,6);
    
    for j = 1:Nt
        fprintf('testing : N = %i', N(Index))
        %% using GloptiPloy 3 to construct random polynomials
        tic
        n = N(Index);
        %% generating POP
        mpol('x',n,1)
        b  = mmon(x,0,2*d-1);    % degree up tp 2d - 1;
        p0 = rand(1,length(b));%sprandn(1,length(b),0.2);
        p0 = p0/norm(p0);
        p = p0*b - sum(x.^(2*d));
        g = sum(x.^(2));
        K = [g <=1];
        P = msdp(min(p), K, d);
        [A,b,c,K] = msedumi(P);

        %% Problem data in Sedumi form
        [m,n] = size(A);
        Density(Index,:) = [m,n,sum(sum(spones(A)))/m/n,(min(K.s)*(min(K.s)+1))/2];

        TimePro(Index) = toc;
        At = A';
        data{Index}.At = At;
        data{Index}.b = b;
        data{Index}.c = c;
        data{Index}.K = K;
        %% solutions using different method

        % 1 by sedumi
        x1 = zeros(length(c),1);
        x2 = zeros(length(c),1);
        x3 = zeros(length(c),1);
        x4 = zeros(length(c),1);
        if N(Index) < 20
        try
            opts.eps = TolIPM;
            [x1,y1,infoSedumi] = sedumi(At,b,c,K,opts);
        catch
            warning('out of memory')
        end
        end

        if N(Index) < 24
        % 2 by sdpt3
        try
            [blk,At2,C2,b2]      = read_sedumi(At,b,c,K); 
            tsdpt = tic;
            opts.gaptol = TolIPM;
            [obj,X,y,Z,infoSDPT] = sqlp(blk,At2,C2,b2,opts); 
            timeSdpt = toc(tsdpt);
        catch
            warning('out of memory')
        end

        % 3 by sdpa
        if N(Index) < 25  % N = 29, matlab on my computer crashed
        try  % sedumiwrap cannot handle k.q
            if isfield(K,'q') && ~isempty(K.q) && K.q == 0
                K.q = [];  
                tSdpa = tic;
                opts.epsilonStar = TolIPM;
                [x3,y3,infoSdpa] = sedumiwrap(A,b,c,K,[],opts);
                timeSdpa = toc(tSdpa);
            else

            end
        catch
            warning('out of memory')
        end
        end

        % 4 by csdp
        try
            tCsdp    = tic;
            opts.axtol = TolIPM;
            opts.aytol = TolIPM;
            opts.objtol = TolIPM;
            [x4,y4,z4,info_csdp] = csdp(At,b,c,K,opts); 
            timeCsdp = toc(tCsdp);
        catch
            warning('out of memory')
        end
        end

        % 5 by SCS-direct
        params.max_iters = Maxiter;
        params.eps = Tol;
        [x5,y5,cscs5,infoSCSdirect] = solveWithSCSdirect(At,full(b),full(c),K,params);

        % 6 by SCS-indirect
        params.max_iters = Maxiter;
        params.eps = Tol;
        [x6,y6,cscs6,infoSCSindirect] = solveWithSCSindirect(At,full(b),full(c),K,params);


        % 7 by cdcs - sos
        opts.relTol = Tol;
        opts.solver = 'sos';
        opts.maxIter = Maxiter;
        [x7,y7,z7,infoCDCSsos] = cdcs(At,b,c,K,opts);
        
        %%
        TimeTotalj(j,:) = [infoSedumi.wallsec,timeSdpt,timeSdpa,timeCsdp, ...
                          (infoSCSdirect.solveTime+infoSCSdirect.setupTime)/1e3, ...
                          (infoSCSindirect.solveTime+infoSCSindirect.setupTime)/1e3, ...
                          infoCDCSsos.time.total];   
        TimeSetupj(j,:) = [infoSCSdirect.setupTime/1e3,infoSCSindirect.setupTime/1e3,infoCDCSsos.time.setup]; 
        TimeADMMj(j,:) = [infoSCSdirect.solveTime/1e3,infoSCSindirect.solveTime/1e3,infoCDCSsos.time.admm]; 
        Costj(j,:) = [c'*x1,obj(1),c'*x3,c'*x4,cscs5,cscs6,c'*x7];
        Iterj(j,:) = [infoSedumi.iter,infoSDPT.iter,infoSdpa.iteration,infoSCSdirect.iter,infoSCSindirect.iter,infoCDCSsos.iter];
        TimeAverj(j,:)  = TimeADMMj(j,:)./Iterj(j,4:end);

    end
    
    TimeTotal(Index,:) = mean(TimeTotalj);
    TimeSetup(Index,:) = mean(TimeSetupj); 
    TimeADMM(Index,:) = mean(TimeADMMj);
    Cost(Index,:) = mean(Costj);
    Iter(Index,:) = mean(Iterj);
    TimeAver(Index,:)  = mean(TimeAverj);
    
    TimeTotalV(Index,:) = std(TimeTotalj); 
    TimeADMMV(Index,:)  = std(TimeADMMj);
    TimeAverV(Index,:)  = std(TimeAverj);  
    TimeSetupV(Index,:) = std(TimeSetupj); 

    save RandomPolyN5.mat
end

%save RandomPolyN.mat

