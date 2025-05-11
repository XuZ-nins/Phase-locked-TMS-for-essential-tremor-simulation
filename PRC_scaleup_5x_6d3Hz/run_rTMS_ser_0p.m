clear;clc;
% Initialize synaptic connections between neurons
% createConnections;clear;

if ispc
	copyfile('nrnmech.dll','rTMS_ser_0p/ET_ser_hoc');
else
	system('cp -r x86_64 rTMS_ser_0p/ET_ser_hoc');
end

rng(1);

%% Real-Time Phase Adaptation
% Initialize rTMS estimation-related variables

for mm = 1
	disp(mm);
	mkdir(strcat('rTMS_ser_0p/simulation_',num2str(mm),'/recordings_full'));
	IONstim = 2200;

	system(strcat('cp -r rTMS_ser_0p/ET_ser_hoc/*',' rTMS_ser_0p/simulation_',num2str(mm),'/'));
	system(strcat('cp initialization/currentState.dat rTMS_ser_0p/simulation_',num2str(mm),'/'));

	% Specify paths for storing simulation files
	formatSpec = '%s \n';
	opt_fileID = fopen(strcat('rTMS_ser_0p/simulation_',num2str(mm),'/trialname_ser.txt'),'w');
	f_PCap = 'recordings_full/PCap.txt';
	f_ION = 'recordings_full/ION.txt';
	f_DCNap = 'recordings_full/DCNap.txt';
	f_Vimap = 'recordings_full/Vimap.txt';
	fprintf(opt_fileID,formatSpec,f_PCap);
	fprintf(opt_fileID,formatSpec,f_ION);
	fprintf(opt_fileID,formatSpec,f_DCNap);
	fprintf(opt_fileID,formatSpec,f_Vimap);
	fclose(opt_fileID);

	formatSpec = '%.8f \n';
	opt_fileID = fopen(strcat('rTMS_ser_0p/simulation_',num2str(mm),'/IONserams.txt'),'w');
	fprintf(opt_fileID,formatSpec,IONstim);
	fclose(opt_fileID);

	cd(strcat('rTMS_ser_0p/simulation_',num2str(mm),'/'));
	[status,cmdout] = system('nrniv -nogui ET_rTMS.hoc');
	disp(cmdout);
	cd('../../');

	disp('Done');
end