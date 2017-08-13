
clear;clc;
load('RandomPoly.mat','data')

n = length(data);

TimeTotal  = zeros(n,4);    % blk, ldl, blk-ldl, blk-chol
TimeLinear = zeros(n,4);
TimeConic  = zeros(n,4);
Iter       = zeros(n,4);

opts1.solver = 'sos';
opts1.maxIter = 200;
opts1.KKTfact = 'blk';

opts2.solver = 'sos';
opts2.maxIter = 200;
opts2.KKTfact = 'ldl';

opts3.solver = 'sos';
opts3.maxIter = 200;
opts3.KKTfact = 'blk-ldl';

opts4.solver = 'sos';
opts4.maxIter = 200;
opts4.KKTfact = 'blk-chol';
for i = 1:n
    
    [x,y,z,info2] = cdcs(data{i}.At,data{i}.b,data{i}.c,data{i}.K,opts2);
    [x,y,z,info3] = cdcs(data{i}.At,data{i}.b,data{i}.c,data{i}.K,opts3);
    [x,y,z,info4] = cdcs(data{i}.At,data{i}.b,data{i}.c,data{i}.K,opts4);
    [x,y,z,info1] = cdcs(data{i}.At,data{i}.b,data{i}.c,data{i}.K,opts1);
    
    
    TimeTotal(i,:) = [info1.time.total,info2.time.total,info3.time.total,info4.time.total];
    TimeLinear(i,:) = [info1.time.subiter(1),info2.time.subiter(1),info3.time.subiter(1),info4.time.subiter(1)];
    TimeConic(i,:) = [info1.time.subiter(2),info2.time.subiter(2),info3.time.subiter(2),info4.time.subiter(2)];
    Iter(i,:)      = [info1.iter,info2.iter,info3.iter,info4.iter];
    
end

save ResultRomLinProj