clear;clc;
% Initialize synaptic connections between neurons
% createConnections;clear;

if ispc
	copyfile('nrnmech.dll','OL_tACS_par/ET_par_hoc');
else
	system('cp -r x86_64 OL_tACS_par/ET_par_hoc');
end

par_on = 1;  % set to 1 to use parfor loop
delete(gcp('nocreate'));
addpath(pwd);

run_range = 1:50;
ampparam_all = 0.0001:0.0001:0.005; % Stimulation amplitude (in nA)

% Improve the accuracy of phase tracking with a time lag (in ms; optional);
% however, it reduces the responsiveness of phase-locked stimulation
tremorT = 158;

IONstim = 2200;

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

    parfor mm = run_range
		if exist(strcat(par_dir,'/OL_tACS_par/simulation_',num2str(mm)),'dir')
			continue;
		end
		
        ampparam = ampparam_all(mm);
        fprintf('ampparam: %f\n',ampparam);

		cd(par_dir);

        mkdir(strcat(par_dir,'/OL_tACS_par/simulation_',num2str(mm),'/recordings_full'));
        system(strcat('cp -r OL_tACS_par/ET_par_hoc/*',' OL_tACS_par/simulation_',num2str(mm),'/'));
        system(strcat('cp -r initialization/currentState.dat OL_tACS_par/simulation_',num2str(mm),'/'));
		
		cd(strcat(par_dir,'/OL_tACS_par/simulation_',num2str(mm),'/'));
		
        % Specify paths for storing simulation files
        formatSpec = '%s \n';
		opt_fileID = fopen(strcat(par_dir,'/OL_tACS_par/simulation_',num2str(mm),'/trialname_par.txt'),'w');
		f_PCap = 'recordings_full/PCap.txt';
		f_ION = 'recordings_full/ION.txt';
		f_DCNap = 'recordings_full/DCNap.txt';
		f_Vimap = 'recordings_full/Vimap.txt';
		fprintf(opt_fileID,formatSpec,f_PCap);
		fprintf(opt_fileID,formatSpec,f_ION);
		fprintf(opt_fileID,formatSpec,f_DCNap);
		fprintf(opt_fileID,formatSpec,f_Vimap);
		fclose(opt_fileID);
		
        % Update RNG seeds that determine synaptic/membrane noises
		dlmwrite('./IONparams.txt',IONstim);
		dlmwrite('./sinparams.txt',[ampparam;1000/tremorT;(rand-0.5)*2*pi]);

		% Execute NEURON
		[status,cmdout] = system(strcat('nrniv -nogui ET_OL_tACS.hoc'));
		% disp(cmdout);
		
        cd(par_dir);
		disp('Done');
    end
else
    disp('Not enough resources for parallel execution.');
end
