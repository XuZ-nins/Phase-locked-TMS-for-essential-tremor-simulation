clear;clc;
% Initialize synaptic connections between neurons
% createConnections;clear;

if ispc
	copyfile('nrnmech.dll','irTMS_par_15p_fopt/ET_par_hoc');
else
	system('cp -r x86_64 irTMS_par_15p_fopt/ET_par_hoc');
end

rng(1);

Leaplist = [11 6 67 43 45 48 20 78 42 84];
tremor_T = 139;
Initlist = round(linspace(0,tremor_T,11));
Initlist = Initlist(1:10);
par_on = 1;
delete(gcp('nocreate'));

%% Real-Time Phase Adaptation
% Initialize irTMS estimation-related variables
tseries = 1:1e3;
load('iTMS.mat');
dt = 0.0125; 
tstop = 10000;
tvec = (0:dt:tstop);

y_all = [];

for k = 1:length(Leaplist)
    for l = 1:length(Initlist)
        x = [Leaplist(k);Initlist(l)];
        y = [k,l];
        y_all = [y_all;y];
    end
end

% Leapcounter = 1;
ampparam = 1; % Stimulation amplitude (in nA)
stimlag = 0; % Stimulation lag

numCPUs = feature('numCores');
fprintf('Number of CPUs requested = %g\n',numCPUs);
par_dir = pwd;
if numCPUs > 1 && par_on
    if numCPUs > 4 % on cluster
        pc_storage_dir = fullfile(par_dir,'pc_storage',getenv('SLURM_JOB_ID'));
        mkdir(pc_storage_dir);
        pc = parcluster('local');
        pc.JobStorageLocation =  pc_storage_dir;
    else
        pc = parcluster('local'); % use default
    end
    if numCPUs>50
        numCPUs = 50;
    end
    poolobj = parpool(pc,numCPUs);

    parfor mm = 1:size(y_all,1)
    
		if exist(strcat('irTMS_par_15p_fopt/simulation_',num2str(mm),'/recordings_full/Vimap.txt'),'file')
			continue;
		end

        disp(mm);
        y = y_all(mm,:);
        tmpLeap = Leaplist(y(1));
        tmpInit = Initlist(y(2));

    	mkdir(strcat('irTMS_par_15p_fopt/simulation_',num2str(mm),'/recordings_full'));
        IONstim = 2200;
		
		y = y_all(mm,:);
		skipnum = round(rand*100);
		fprintf('tmpLeap: %f\n',tmpLeap);
		fprintf('tmpInit: %f\n',tmpInit);
		p = irTMSset(1,'Leap',tmpLeap);
		irTMSseq = net(p,128);
		
        delparam = tmpInit+IONstim+ceil(tremor_T);
		stimlist = [delparam;zeros(100,1)];
		for ll = 2:101
			stimlist(ll) = round(stimlist(ll-1)+tremor_T*1.1*(1+irTMSseq(ll)-irTMSseq(ll-1)));
		end
		stimlist = stimlist(2:end) + stimlag;
		
		stimlist(stimlist>=tstop) = [];
		iflag = zeros(length(tvec),1);
		
	%     stimlist = stimlist - 2600;
		iflag(stimlist/dt+1) = 1;
		ivec = conv(iflag,irec);     
		if length(tvec) > length(ivec)
			ivec = [ivec;zeros(length(tvec)-length(ivec),1)]; % add trailing zeros
		else
			ivec = ivec(1:length(tvec));
		end
    	writeVectorBin(strcat('irTMS_par_15p_fopt/simulation_',num2str(mm),'/'),tvec,ivec); % write tvec, Evec to params/test

    	system(strcat('cp -r irTMS_par_15p_fopt/ET_par_hoc/*',' irTMS_par_15p_fopt/simulation_',num2str(mm),'/'));
    	system(strcat('cp initialization/currentState.dat irTMS_par_15p_fopt/simulation_',num2str(mm),'/'));

        % Specify paths for storing simulation files
        formatSpec = '%s \n';
    	opt_fileID = fopen(strcat('irTMS_par_15p_fopt/simulation_',num2str(mm),'/trialname_par.txt'),'w');
        f_PCap = strcat('recordings_full/PCap.txt');
        f_ION = strcat('recordings_full/ION.txt');
        f_DCNap = strcat('recordings_full/DCNap.txt');
        f_Vimap = strcat('recordings_full/Vimap.txt');
        f_TMS = strcat('recordings_full/TMS.txt');
        f_15p = strcat('recordings_full/iTMS.txt');
        fprintf(opt_fileID,formatSpec,f_PCap);
        fprintf(opt_fileID,formatSpec,f_ION);
        fprintf(opt_fileID,formatSpec,f_DCNap);
        fprintf(opt_fileID,formatSpec,f_Vimap);
        fprintf(opt_fileID,formatSpec,f_TMS);
        fprintf(opt_fileID,formatSpec,f_15p);
        fclose(opt_fileID);

        formatSpec = '%.8f \n';
    	opt_fileID = fopen(strcat('irTMS_par_15p_fopt/simulation_',num2str(mm),'/IONparams.txt'),'w');
        fprintf(opt_fileID,formatSpec,IONstim);
        fclose(opt_fileID);

        formatSpec = '%d \n';
    	opt_fileID = fopen(strcat('irTMS_par_15p_fopt/simulation_',num2str(mm),'/prcparams.txt'),'w');
        fprintf(opt_fileID,formatSpec,stimlist);
        fclose(opt_fileID);
        system(strcat('cp irTMS_par_15p_fopt/simulation_',num2str(mm),'/prcparams.txt irTMS_par_15p_fopt/simulation_',num2str(mm),'/recordings_full/stimlist.txt'))

    	cd(strcat('irTMS_par_15p_fopt/simulation_',num2str(mm),'/'));
    	[status,cmdout] = system('nrniv -nogui ET_irTMS.hoc');
    	% disp(cmdout);
		cd('../../');

        disp('Done');
    end
end