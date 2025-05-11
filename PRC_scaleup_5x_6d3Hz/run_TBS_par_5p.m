clear;clc;
% Initialize synaptic connections between neurons
% createConnections;clear;

if ispc
	copyfile('nrnmech.dll','TBS_par_5p/ET_par_hoc');
else
	system('cp -r x86_64 TBS_par_5p/ET_par_hoc');
end

rng(1);

% Leapcounter = 1;
ampparam = 1; % Stimulation amplitude (in nA)
stimlag = 0; % Stimulation lag
tremor_T = 158;
Initlist = round(linspace(0,tremor_T,11));
Initlist = Initlist(1:10);
par_on = 1;
delete(gcp('nocreate'));

%% Real-Time Phase Adaptation
% Initialize rTMS estimation-related variables
tseries = 1:1e3;
load('iTMS.mat');
dt = 0.0125;
tstop = 10000;
tvec = (0:dt:tstop);

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

    parfor mm = 1:length(Initlist)
	
		if exist(strcat('TBS_par_5p/simulation_',num2str(mm),'/recordings_full/Vimap.txt'),'file')
			continue;
		end

        disp(mm);
        tmpInit = Initlist(mm);

    	mkdir(strcat('TBS_par_5p/simulation_',num2str(mm),'/recordings_full'));
        fprintf('Stimulation onset: %f ms\n',tmpInit);
        IONstim = 2200;

        delparam = tmpInit+IONstim+ceil(tremor_T);
        stimlist = delparam:200:10000;
		stimlist = [stimlist,stimlist+20,stimlist+40];
        stimlist = stimlist(2:end) + stimlag;
        stimlist(stimlist>=tstop) = [];
        iflag = zeros(length(tvec),1);

        iflag(round(stimlist/dt)+1) = 1;
    	ivec = conv(iflag,irec);
    	if length(tvec) > length(ivec)
    		ivec = [ivec;zeros(length(tvec)-length(ivec),1)]; % add trailing zeros
    	else
    		ivec = ivec(1:length(tvec));
    	end
    	writeVectorBin(strcat('TBS_par_5p/simulation_',num2str(mm),'/'),tvec,ivec); % write tvec, Evec to params/test

    	system(strcat('cp -r TBS_par_5p/ET_par_hoc/*',' TBS_par_5p/simulation_',num2str(mm),'/'));
    	system(strcat('cp initialization/currentState.dat TBS_par_5p/simulation_',num2str(mm),'/'));

        % Specify paths for storing simulation files
        formatSpec = '%s \n';
    	opt_fileID = fopen(strcat('TBS_par_5p/simulation_',num2str(mm),'/trialname_par.txt'),'w');
        f_PCap = strcat('recordings_full/PCap.txt');
        f_ION = strcat('recordings_full/ION.txt');
        f_DCNap = strcat('recordings_full/DCNap.txt');
        f_Vimap = strcat('recordings_full/Vimap.txt');
        f_TMS = strcat('recordings_full/TMS.txt');
        f_5p = strcat('recordings_full/iTMS.txt');
        fprintf(opt_fileID,formatSpec,f_PCap);
        fprintf(opt_fileID,formatSpec,f_ION);
        fprintf(opt_fileID,formatSpec,f_DCNap);
        fprintf(opt_fileID,formatSpec,f_Vimap);
        fprintf(opt_fileID,formatSpec,f_TMS);
        fprintf(opt_fileID,formatSpec,f_5p);
        fclose(opt_fileID);

        formatSpec = '%.8f \n';
    	opt_fileID = fopen(strcat('TBS_par_5p/simulation_',num2str(mm),'/IONparams.txt'),'w');
        fprintf(opt_fileID,formatSpec,IONstim);
        fclose(opt_fileID);

        formatSpec = '%d \n';
    	opt_fileID = fopen(strcat('TBS_par_5p/simulation_',num2str(mm),'/prcparams.txt'),'w');
        fprintf(opt_fileID,formatSpec,stimlist);
        fclose(opt_fileID);
        system(strcat('cp TBS_par_5p/simulation_',num2str(mm),'/prcparams.txt TBS_par_5p/simulation_',num2str(mm),'/recordings_full/stimlist.txt'))

    	cd(strcat('TBS_par_5p/simulation_',num2str(mm),'/'));
    	[status,cmdout] = system('nrniv -nogui ET_TBS.hoc');
    	% disp(cmdout);
		cd('../../');

        disp('Done');
    end
end