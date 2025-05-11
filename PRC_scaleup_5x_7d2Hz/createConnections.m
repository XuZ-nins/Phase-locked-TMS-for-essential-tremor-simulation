% Create randomized connections among neurons in the 5 CCTC loops
% Author: Xu Zhang @UConn, Apr., 2019

clear;clc;
% rng(479);
% Initialize rng seed
rngSEED = ceil(rand*899+100);
rng(rngSEED);
disp(rngSEED);

%% ION gap junctions
subGroupMat = [0 1 1 0 1 0 0 0;1 0 0 1 0 1 0 0;1 0 0 1 0 0 1 0;0 1 1 0 0 0 0 1;1 0 0 0 0 1 1 0;0 1 0 0 1 0 0 1;0 0 1 0 1 0 0 1;0 0 0 1 0 1 1 0];
subIONgapmat = blkdiag(subGroupMat,subGroupMat,subGroupMat,subGroupMat,subGroupMat);
for i = 11:15
    subIONgapmat(rem(i-1,5)*8+1,rem(i-2,5)*8+2)=1;
    subIONgapmat(rem(i-1,5)*8+2,rem(i,5)*8+1)=1;
    subIONgapmat(rem(i-1,5)*8+3,rem(i-2,5)*8+4)=1;
    subIONgapmat(rem(i-1,5)*8+4,rem(i,5)*8+3)=1;
end
IONgapmat = blkdiag(subIONgapmat,subIONgapmat,subIONgapmat,subIONgapmat,subIONgapmat);
% IONgapmat is the connectivity matrix within ION cells

IONgap_all = [];
for i = 1:size(IONgapmat,1)
    IONgap_all = [IONgap_all;find(IONgapmat(i,:)==1)'];
end

% Set gap junction connectivity
formatSpec = '%d \n';
PYNint_fileID = fopen('params_ION_gapj.txt','w');
fprintf(PYNint_fileID,formatSpec,length(IONgap_all));
fprintf(PYNint_fileID,formatSpec,sum(IONgapmat,2));
fprintf(PYNint_fileID,formatSpec,IONgap_all);
fclose(PYNint_fileID);

% Set gap junction conductance
formatSpec = '%.4f \n';
ION_gapstrength = rand(length(IONgap_all),1)*0.4+0.3;
IONgaps_fileID = fopen('params_ION_gapstrength.txt','w');
fprintf(IONgaps_fileID,formatSpec,ION_gapstrength);
fclose(IONgaps_fileID);

% Set offset currents to ION
formatSpec = '%.4f \n';
ION_oc = 1e-3*(randn(40*5,1)*0.02-1.3);
IONoc_fileID = fopen('params_ION_oc.txt','w');
fprintf(IONoc_fileID,formatSpec,ION_oc);
fclose(IONoc_fileID);

%% PYN connections
modelscale = 5; % 5 loops
width = 2.5; subPYNintmat = tril(triu(true(20*modelscale),-ceil(width)),floor(width))-diag(true(20*modelscale,1));
for i=1:ceil(width)
    subPYNintmat(i,end-ceil(width)+i:end) = true; subPYNintmat(end-floor(width)+i:end,i) = true;
end
PYNintmat = blkdiag(subPYNintmat,subPYNintmat,subPYNintmat,subPYNintmat,subPYNintmat);
% PYNintmat is the connectivity matrix within PYN cells

PYNint_all = [];
for i = 1:length(PYNintmat)
    PYNint_all = [PYNint_all;find(PYNintmat(i,:)==1)'];
end

formatSpec = '%d \n';
PYNint_fileID = fopen('params_PY_int.txt','w');
fprintf(PYNint_fileID,formatSpec,length(PYNint_all));
fprintf(PYNint_fileID,formatSpec,PYNint_all);
fclose(PYNint_fileID);

%% ION-PC connections
% Randomized connections with constraints, see SI Note 3 2)
IONPCmat = zeros(40*5,200*5);
groupRatioIONPC = round(randn(24,1)*4+40);
groupRatioIONPC(25) = 200*5-sum(groupRatioIONPC(1:24));
while min(groupRatioIONPC)<=25 || max(groupRatioIONPC)>=55 || length(find(groupRatioIONPC>40))>=20 || length(find(groupRatioIONPC<40))>=20
    groupRatioIONPC = round(randn(24,1)*4+40);
    groupRatioIONPC(25) = 200*5-sum(groupRatioIONPC(1:24));
end

IONPCcount = 0;
for m = 1:25
    subRatioIONPC = abs(round(randn(7,1)*2+5));
    subRatioIONPC(8) = groupRatioIONPC(m)-sum(subRatioIONPC(1:7));
    while min(subRatioIONPC)<=1 || max(subRatioIONPC)>=10 || length(find(subRatioIONPC>=7))>=6 || length(find(subRatioIONPC<=3))>=6
        subRatioIONPC = abs(round(randn(7,1)*2+5));
    subRatioIONPC(8) = groupRatioIONPC(m)-sum(subRatioIONPC(1:7));
    end
    for n = 1:8
        IONPCmat((m-1)*8+n,(IONPCcount+1):(IONPCcount+subRatioIONPC(n))) = 1;
        IONPCcount = IONPCcount + subRatioIONPC(n);
    end
end

% IONPCmat is the connectivity matrix between ION->PC

IONPC_all = [];
for i = 1:size(IONPCmat,1)
    IONPC_all = [IONPC_all;find(IONPCmat(i,:)==1)'];
end

formatSpec = '%d \n';
PYNint_fileID = fopen('params_ION2PC.txt','w');
fprintf(PYNint_fileID,formatSpec,sum(IONPCmat,2));
fprintf(PYNint_fileID,formatSpec,IONPC_all);
fclose(PYNint_fileID);

%% PC-DCN connections
% Randomized connections with constraints, see SI Note 3 3)
groupPCDCN = round(randn(24,1)*4+40);
groupPCDCN(25) = 200*5-sum(groupPCDCN(1:4));
while min(groupPCDCN)<=25 || max(groupPCDCN)>=55 || length(find(groupPCDCN>40))>=20 || length(find(groupPCDCN<40))>=20
    groupPCDCN = round(randn(24,1)*4+40);
    groupPCDCN(25) = 200*5-sum(groupPCDCN(1:24));
end

PCDCNmat = zeros(200*5,5*5);
PCDCNcount = 0;
for m = 1:25
    PCDCNmat((PCDCNcount+1):(PCDCNcount+groupPCDCN(m)),m) = 1;
    PCDCNcount = PCDCNcount + groupPCDCN(m);
end

% PCDCNmat is the connectivity matrix between PC->DCN

PCDCN_all = [];
for i = 1:size(PCDCNmat,1)
    PCDCN_all = [PCDCN_all;find(PCDCNmat(i,:)==1)'];
end

formatSpec = '%d \n';
PYNint_fileID = fopen('params_PC2DCN.txt','w');
fprintf(PYNint_fileID,formatSpec,sum(PCDCNmat,2));
fprintf(PYNint_fileID,formatSpec,PCDCN_all);
fclose(PYNint_fileID);

%% DCN-ION connections
% Randomized connections with constraints, see SI Note 3 4)
groupRatioDCNION = round(randn(24,1)*4+8);
groupRatioDCNION(25) = 40*5-sum(groupRatioDCNION(1:24));
while min(groupRatioDCNION)<=3 || max(groupRatioDCNION)>=13 || length(find(groupRatioDCNION>8))>=20 || length(find(groupRatioDCNION<8))>=20
    groupRatioDCNION = round(randn(24,1)*4+8);
    groupRatioDCNION(25) = 40*5-sum(groupRatioDCNION(1:24));
end

DCNIONcount = 0;
currentIONlist = 1:40*5;
DCNIONmat = zeros(5*5,40*5);
for m = 1:24
    tmprange = length(currentIONlist);
    tmpIONind = ceil(rand(groupRatioDCNION(m),1)*(tmprange-1))+1;
    while length(unique(tmpIONind))<length(tmpIONind)
        tmpIONind = ceil(rand(groupRatioDCNION(m),1)*(tmprange-1))+1;
    end
    DCNIONmat(m,currentIONlist(tmpIONind)) = 1;
    currentIONlist(tmpIONind)=[];
end

DCNIONmat(25,currentIONlist) = 1;

% DCNIONmat is the connectivity matrix between DCN->ION

DCNION_all = [];
for i = 1:size(DCNIONmat,1)
    DCNION_all = [DCNION_all;find(DCNIONmat(i,:)==1)'];
end

formatSpec = '%d \n';
PYNint_fileID = fopen('params_DCN2ION.txt','w');
fprintf(PYNint_fileID,formatSpec,sum(DCNIONmat,2));
fprintf(PYNint_fileID,formatSpec,DCNION_all);
fclose(PYNint_fileID);