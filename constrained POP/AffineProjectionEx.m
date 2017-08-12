

clc;clear;close all;

load QuarCase.mat

n = length(data);

TimeTotal  = zeros(n,3);    % blk, ldl, blk-ldl
TimeLinear = zeros(n,3);
TimeConic  = zeros(n,3);
Iter       = zeros(n,3);

opts1.solver = 'sos';
opts1.maxIter = 200;
opts1.KKTfact = 'blk';

opts2.solver = 'sos';
opts2.maxIter = 200;
opts2.KKTfact = 'ldl';

opts3.solver = 'sos';
opts3.maxIter = 200;
opts3.KKTfact = 'blk-ldl';

for i = 1:n
    
    [x,y,z,info2] = cdcs(data{i}.At,data{i}.b,data{i}.c,data{i}.K,opts2);
    [x,y,z,info3] = cdcs(data{i}.At,data{i}.b,data{i}.c,data{i}.K,opts3);
    [x,y,z,info1] = cdcs(data{i}.At,data{i}.b,data{i}.c,data{i}.K,opts1);
    
    
    TimeTotal(i,:) = [info1.time.total,info2.time.total,info3.time.total];
    TimeLinear(i,:) = [info1.time.subiter(1),info2.time.subiter(1),info3.time.subiter(1)];
    TimeConic(i,:) = [info1.time.subiter(2),info2.time.subiter(2),info3.time.subiter(2)];
    Iter(i,:)      = [info1.iter,info2.iter,info3.iter];
    
end

save ResultLinProj




