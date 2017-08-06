

%% solving the sdps arising from Lyapunov tests

clc;clear;close all;

warning on
load LyapunovSeDuMiData.mat

TimeTotal = zeros(length(n),7);   % sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos
TimeSetup = zeros(length(n),3);   % scs-direct, scs-indirect, cdcs-sos
TimeADMM = zeros(length(n),3);    % scs-direct, scs-indirect, cdcs-sos
TimeAver = zeros(length(n),3);    % scs-direct, scs-indirect, cdcs-sos
Cost = zeros(length(n),7);        % sedumi, sdpt3, sdpa, csdp, scs-direct, scs-indirect, cdcs-sos
Iter = zeros(length(n),6);        % sedumi, sdpt3, sdpa, scs-direct, scs-indirect, cdcs-sos
Density = zeros(length(n),4);
InfoTerm = zeros(length(n),6);

%%
Maxiter = 2e3;
Tol     = 1e-3;
TolIPM  = 1e-8;

for i = 1:9
    
    At = SeDuMiData{i}.At;
    b = SeDuMiData{i}.b;
    c = SeDuMiData{i}.c;
    K = SeDuMiData{i}.K;
    %% solutions using different method
    
    % 1 by sedumi
    x1 = zeros(length(c),1)+NaN;
    x2 = zeros(length(c),1)+NaN;
    x3 = zeros(length(c),1);
    try
        opts.eps = TolIPM;
        [x1,y1,infoSedumi] = sedumi(At,b,c,K,opts);
    catch
        warning('out of memory')
    end
    
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
    if n(i) < 29  % N = 29, matlab on my computer crashed
        try  % sedumiwrap cannot handle k.q
            if ~isfield(K,'q') || isempty(K.q) ||(~isempty(K.q) && K.q == 0)
                K.q = [];  
                tSdpa = tic;
                [x3,y3,infoSdpa] = sedumiwrap(At',b,c,K,[],opts);
                timeSdpa = toc(tSdpa);
            end
        catch
            warning('out of memory')
        end
    end
    
    % 4 by csdp
    try
        opts.axtol = TolIPM;
        opts.aytol = TolIPM;
        opts.objtol = TolIPM;
        if (isfield(K,'f')) %Convert free vars to non-negative LP vars
                n_free = K.f; 
                [A1,b1,c1,K1] = convertf(At,b,c,K); %K.f set to zero
                At1 = A1';
        end
            c1 = full(c1);
            tCsdp    = tic;
            [x4,y4,z4,info_csdp] = csdp(At1,b1,c1,K1,opts);  %JA updated handling of info flag
            timeCsdp = toc(tCsdp);
    catch
        warning('out of memory')
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
    
    
    %% statistics

    TimeTotal(i,:) = [infoSedumi.wallsec,timeSdpt,timeSdpa,timeCsdp, ...
                          (infoSCSdirect.solveTime+infoSCSdirect.setupTime)/1e3, ...
                          (infoSCSindirect.solveTime+infoSCSindirect.setupTime)/1e3, ...
                          infoCDCSsos.time.total];   
    TimeSetup(i,:) = [infoSCSdirect.setupTime/1e3,infoSCSindirect.setupTime/1e3,infoCDCSsos.time.setup]; 
    TimeADMM(i,:) = [infoSCSdirect.solveTime/1e3,infoSCSindirect.solveTime/1e3,infoCDCSsos.time.admm]; 
    Cost(i,:) = [c'*x1,obj(1),c'*x3,c1'*x4,cscs5,cscs6,c'*x7];
    Iter(i,:) = [infoSedumi.iter,infoSDPT.iter,infoSdpa.iteration,infoSCSdirect.iter,infoSCSindirect.iter,infoCDCSsos.iter];
    TimeAver(i,:)  = TimeADMM(i,:)./Iter(i,4:end);
    InfoTerm(i,:) = [infoSedumi.feasratio,infoSDPT.termcode, infoSdpa.dualityGap,infoSCSdirect.relGap,infoSCSindirect.relGap,infoCDCSsos.problem] ; 

    save LyapunovSeDuMiDataR1.mat
end