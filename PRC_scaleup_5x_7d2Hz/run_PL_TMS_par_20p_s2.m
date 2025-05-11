clear;clc;
% Initialize synaptic connections between neurons (already generated)
% createConnections; clear;

if ispc
	copyfile('nrnmech.dll','PL_TMS_par_20p_s2/ET_par_hoc');
else
	system('cp -r x86_64 PL_TMS_par_20p_s2/ET_par_hoc');
end

delete(gcp('nocreate'));

addpath(pwd);
% Load TMS current shape to be used as IClamp in NEURON simulation
load('iTMS.mat');
PRClist = -pi:pi/25:pi-1e-5; % List of phases to be locked onto by TMS/tACS

t_stop = 10000; % Simulation period (ms)
t_step = 0.25; % Time step (ms)

run_range = 1:length(PRClist); % Select which phases to run

dt = 0.0125; % dt in NEURON simulation
tvec = (0:dt:t_stop);

%% Real-Time Phase Adaptation
% Initialize phase estimation-related variables
tseries = t_step:t_step:1e3; % Horizon to calculate real-time phase values

% Improve the accuracy of phase tracking with a time lag (in ms; optional);
% however, it reduces the responsiveness of phase-locked stimulation
phalag = 0;
tremorT = 139; % 7.2 Hz tremor

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
		if exist(strcat(par_dir,'/PL_TMS_par_20p_s2/simulation_',num2str(mm)),'dir')
			continue;
		end
		
        delaycounter = 0;
        initflg = 2;
        tmpPRC = PRClist(mm);
        fprintf('tmpPRC: %f\n',tmpPRC);

		cd(par_dir);

        % Create directories and copy initialized states
		mkdir(strcat(par_dir,'/PL_TMS_par_20p_s2/simulation_',num2str(mm),'/recordings_full'));
        system(strcat('cp -r PL_TMS_par_20p_s2/ET_par_hoc/*',' PL_TMS_par_20p_s2/simulation_',num2str(mm),'/'));
        system(strcat('cp -r initialization/currentState.dat PL_TMS_par_20p_s2/simulation_',num2str(mm),'/'));
        system(strcat('cp -r initialization/* PL_TMS_par_20p_s2/simulation_',num2str(mm),'/recordings_full','/'));
		
		cd(strcat(par_dir,'/PL_TMS_par_20p_s2/simulation_',num2str(mm),'/'));
		
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
				% Real-time estimation of instantaneous phase from Vim
				tmpVim_all = sum(spike_mat_PRC,2)>0;
				tmp_echt = echt(tmpVim_all, 4, 10, 4e3);
				test_PRC_ec = angle(tmp_echt);
				% Update new phase
				tmpphase = test_PRC_ec(end-phalag);

				fprintf('tmpphase_pre: %f\n',tmpphase_pre)
				fprintf('tmpphase: %f\n',tmpphase)
				if (tmpphase>=tmpPRC && tmpphase_pre<=tmpPRC) || ...
						(tmpphase<tmpphase_pre-pi && tmpphase_pre<=tmpPRC) || ...
						(tmpphase<tmpphase_pre-pi && tmpphase>=tmpPRC)
					disp('In phase now');
					
					if i>2200/t_step && delaycounter<i-2.5*tremorT/t_step
						if initflg>0
							initflg = initflg-1;
							fprintf('Time: %f\n',i*t_step);
							disp(strcat(num2str(initflg),' skips left...'));
						else
							disp('Phase updated.');
							delparam = i*t_step+1;
							fprintf('Time: %f\n',i*t_step);
							if ~isempty(stimlist) && delparam<stimlist(end)+tremorT*2.5
								disp('Too close; skipped');
								continue;
							end
							stimlist = [stimlist;delparam];
							fprintf('Next pulse: %g\n',delparam);
							stimlist(stimlist>=t_stop) = [];
							iflag = zeros(length(tvec),1);

							iflag(stimlist/dt+1) = 1;
							ivec = conv(iflag,irec);
							if length(tvec) > length(ivec)
								ivec = [ivec;zeros(length(tvec)-length(ivec),1)]; % add trailing zeros
							else
								ivec = ivec(1:length(tvec));
							end
							writeVectorBin(strcat('./'),tvec,ivec); % write tvec, Evec to params/test
							delaycounter = i;
						end
					end
				end
				tmpphase_pre = tmpphase;
            end
            
            dlmwrite('./rngSeeds.txt',ceil(rand(16,1)*300 + 1000));

            % Execute NEURON
            [status,cmdout] = system(strcat('nrniv -nogui ET_PL_TMS.hoc'));
            % disp(cmdout);
        end
		
        cd(par_dir);

		disp('Done');
    end
else
    disp('Not enough resources for parallel execution.');
end