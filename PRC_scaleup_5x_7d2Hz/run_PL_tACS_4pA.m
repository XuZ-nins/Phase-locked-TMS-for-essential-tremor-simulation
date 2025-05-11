clear;clc;
% Initialize synaptic connections between neurons (already generated)
% createConnections; clear;

if ispc
	copyfile('nrnmech.dll','PL_tACS_par_4pA/ET_par_hoc');
else
	system('cp -r x86_64 PL_tACS_par_4pA/ET_par_hoc');
end

delete(gcp('nocreate'));

addpath(pwd);
PRClist = -pi:pi/25:pi-1e-5; % List of phases to be locked onto by TMS/tACS

t_stop = 10000; % Simulation period (ms)
t_step = 0.25; % Time step (ms)

run_range = 1:length(PRClist); % Select which phases to run
ampparam = 0.004; % Stimulation amplitude (in nA)

dt = 0.0125; % dt in NEURON simulation
tvec = (0:dt:t_stop);

%% Real-Time Phase Adaptation
% Initialize phase estimation-related variables
tseries = t_step:t_step:1e3; % Horizon to calculate real-time phase values

% Improve the accuracy of phase tracking with a time lag (in ms; optional);
% however, it reduces the responsiveness of phase-locked stimulation
tremorT = 139; % 7.2 Hz tremor
phalag = 0;

IONstim = 2200;

currentAmp_all = ampparam; % [zeros(2500/t_step,1);ampparam*ones(10000/t_step,1)];
intrinsicPhase = rem(2*pi/tremorT*(t_step:t_step:t_stop)',2*pi);
intrinsicPhase(intrinsicPhase>pi) = intrinsicPhase(intrinsicPhase>pi) - 2*pi;

numCPUs = feature('numCores');
fprintf('Number of CPUs requested = %g\n',numCPUs);
par_dir = pwd;
if numCPUs > 1
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
		if exist(strcat(par_dir,'/PL_tACS_par_4pA/simulation_',num2str(mm)),'dir')
			continue;
		end
		
        delaycounter = 0;
        tmpPRC = PRClist(mm);
        fprintf('tmpPRC: %f\n',tmpPRC);

		cd(par_dir);

        % Create directories and copy initialized states
		mkdir(strcat(par_dir,'/PL_tACS_par_4pA/simulation_',num2str(mm),'/recordings_full'));
        system(strcat('cp -r PL_tACS_par_4pA/ET_par_hoc/*',' PL_tACS_par_4pA/simulation_',num2str(mm),'/'));
        system(strcat('cp -r initialization/currentState.dat PL_tACS_par_4pA/simulation_',num2str(mm),'/'));
        system(strcat('cp -r initialization/* PL_tACS_par_4pA/simulation_',num2str(mm),'/recordings_full','/'));
		
		cd(strcat(par_dir,'/PL_tACS_par_4pA/simulation_',num2str(mm),'/'));
		
        % Specify paths for storing simulation files
        formatSpec = '%s \n';
        opt_fileID = fopen(strcat('./trialname_par.txt'),'w');
        f_PCap = strcat('./recordings_full/PCap.txt');
        f_DCNap = strcat('./recordings_full/DCNap.txt');
        f_Vimap = strcat('./recordings_full/Vimap.txt');
        f_Vimap_RT = strcat('./recordings_full/Vimap_RT.txt');
        f_ION = strcat('./recordings_full','/ION.txt');
        f_PCv = strcat('./recordings_full/v_PC.txt');
        f_sin = strcat('./recordings_full','/i_sin.txt');
        fprintf(opt_fileID,formatSpec,f_PCap);
        fprintf(opt_fileID,formatSpec,f_DCNap);
        fprintf(opt_fileID,formatSpec,f_Vimap);
        fprintf(opt_fileID,formatSpec,f_Vimap_RT);
        fprintf(opt_fileID,formatSpec,f_ION);
        fprintf(opt_fileID,formatSpec,f_PCv);
        fprintf(opt_fileID,formatSpec,f_sin);
        fclose(opt_fileID);
		
        % Update RNG seeds that determine synaptic/membrane noises
		dlmwrite('./IONparams.txt',IONstim);

        stimlist = [];
        tmpphase_pre = 0;
		
		copyfile(f_Vimap,f_Vimap_RT);
		apVim = dlmread(f_Vimap);
		spike_mat_PRC = zeros(length(tseries),25);
		for nn = 1:25
			tmpspkind = find(apVim(:,1)==nn-1);
			tmpVim = round(apVim(tmpspkind,2)/t_step);
			mmm = find(abs(diff(tmpVim))<=1)+1;
			tmpVim(mmm)=[];
			tmpVim = tmpVim - (2000/t_step-1-length(tseries));
			tmpVim(tmpVim<=0) = [];
			spike_mat_PRC(tmpVim,nn) = 1;
		end

        %% Proceed from 2000 ms (First 2000 ms must have been run)
        for i = (2000/t_step+1):(t_stop/t_step) % 10 seconds
            if mod(i,500/t_step)==0
                disp(i);
				system(strcat('cp -r ./currentState.dat ./recordings_full'));
            end
            
            fprintf('Time = %.2f ms\n',i*t_step);
			
			spike_mat_PRC = [spike_mat_PRC(2:end,:);zeros(1,25)];
			VimapInfo = dir(f_Vimap_RT);
			if ~(isempty(VimapInfo) || VimapInfo.bytes == 0)
				apVim = dlmread(f_Vimap_RT);
				for nn = 1:25
					tmpspkind = find(apVim(:,1)==nn-1);
					tmpVim = round(apVim(tmpspkind,2)/t_step);
					mmm = find(abs(diff(tmpVim))<=1)+1;
					tmpVim(mmm)=[];
					tmpVim = tmpVim - (i-1-length(tseries));
					tmpVim(tmpVim<=0) = [];
					spike_mat_PRC(tmpVim,nn) = 1;
				end
			end
			
            % Real-time update of phase-tracking every 1 ms
            if mod(i,1/t_step)==0
				% Real-time estimation of instantaneous oscillation phase in
				% the Vim through ecHT
				tmpVim_all = sum(spike_mat_PRC,2)>0;
				% tmpVim_all(tmpVim_all>1)=1;
				tmp_echt = echt(tmpVim_all, 4, 10, 4e3);
				% test_amp_ec = abs(tmp_echt);
				test_PRC_ec = angle(tmp_echt);
				% Update new phase
				tmpphase = test_PRC_ec(end-phalag);

				fprintf('tmpphase_pre: %f\n',tmpphase_pre)
				fprintf('tmpphase: %f\n',tmpphase)
				
				% Update new phase
				newphase = tmpphase + tmpPRC - intrinsicPhase(i);
				tmpphase_pre = tmpphase;
				dlmwrite('./sinparams.txt',[currentAmp_all;1000/tremorT;newphase]);
            end
            
            dlmwrite('./rngSeeds.txt',ceil(rand(16,1)*300 + 1000));

            % Execute NEURON
            [status,cmdout] = system(strcat('nrniv -nogui ET_PL_tACS.hoc'));
            % disp(cmdout);
        end
        cd(par_dir);
		disp('Done');
    end
else
    disp('Not enough resources for parallel execution.');
end
