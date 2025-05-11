clear;clc;close all;
addpath(pwd);
tmpdir = pwd;
cd([tmpdir,filesep,'PRC_scaleup_5x_7d2Hz',filesep,'PL_TMS_par_20p_s0']);

for n = [12 46] % Effective vs. Ineffective phases
    %%
    k = round((n+22)/34);
    figure(k);
    subplot(4,1,1);
    apPC = dlmread(strcat('./simulation_',num2str(n),'/recordings_full/PCap.txt'));
    PC_map = zeros(10000,1000);
    hold on;
    for i = 1:1000
        tmpind = find(apPC(:,1)==i-1);
        tmpPC = apPC(tmpind,2);
        yyy = diff(tmpPC);
        mmm = find(abs(yyy)<1)+1;
        tmpPC(mmm)=[];
        tmpPC(tmpPC<1)=[];
        PC_map(round(tmpPC),i)=1;
        plot(tmpPC,ones(length(tmpPC),1)*i,'k.','MarkerSize',2);
    end
    PC_map = PC_map(2001:end,:);
    PC_map_all(:,:,k)=PC_map;
    %
    allStim = abs(dlmread(strcat('./simulation_',num2str(n),'/recordings_full/i_sin.txt')));
    % allStim(1:8000) = []; 
    allStim=allStim>0;
    if length(allStim)<6000
        allStim(end+1:6000)=0;
    end
    yyaxis right
    plot(2000.25:0.25:10000,allStim/max(allStim)*20,'g','LineWidth',2);
    xlim([2500 3500]);ylim([0 100]);
    set(gca,'xtick',[],'ytick',[],'FontSize',14);
    title('PC')

    subplot(4,1,2);

    apDCN = dlmread(strcat('./simulation_',num2str(n),'/recordings_full/DCNap.txt'));
    DCN_map = zeros(10000,25);
    hold on;
    for i = 1:25
        tmpind = find(apDCN(:,1)==i-1);
        tmpDCN = apDCN(tmpind,2);
        yyy = diff(tmpDCN);
        mmm = find(abs(yyy)<1)+1;
        tmpDCN(mmm)=[];
        tmpDCN(tmpDCN<1)=[];
        DCN_map(round(tmpDCN),i)=1;
        plot(tmpDCN,ones(length(tmpDCN),1)*i,'k.','MarkerSize',4);
    end
    DCN_map = DCN_map(2001:end,:);
    DCN_map_all(:,:,k)=DCN_map;
    xlim([2500 3500]);ylim([0 25]);
    title('DCN');ylabel('Number of neurons');
    set(gca,'xtick',[],'FontSize',14);

    subplot(4,1,3);

    apION = dlmread(strcat('./simulation_',num2str(n),'/recordings_full/ION.txt'));
    ION_map = zeros(10000,200);
    hold on;
    for i = 1:200
        tmpind = find(apION(:,1)==i-1);
        tmpION = apION(tmpind,2);
        yyy = diff(tmpION);
        mmm = find(abs(yyy)<1)+1;
        tmpION(mmm)=[];
        tmpION(tmpION<1)=[];
        ION_map(round(tmpION),i)=1;
        plot(tmpION,ones(length(tmpION),1)*i,'k.','MarkerSize',4);
    end
    ION_map = ION_map(2001:end,:);
    ION_map_all(:,:,k)=ION_map;
    xlim([2500 3500]);
    set(gca,'xtick',[],'FontSize',14);
    title('ION')

    subplot(4,1,4);

    apVim = dlmread(strcat('./simulation_',num2str(n),'/recordings_full/Vimap.txt'));
    Vim_map = zeros(10000,25);
    hold on;
    for i = 1:25
        tmpind = find(apVim(:,1)==i-1);
        tmpVim = apVim(tmpind,2);
        yyy = diff(tmpVim);
        mmm = find(abs(yyy)<1)+1;
        tmpVim(mmm)=[];
        tmpVim(tmpVim<1)=[];
        Vim_map(round(tmpVim),i)=1;
        plot(tmpVim,ones(length(tmpVim),1)*i,'k.','MarkerSize',4);
    end
    Vim_map = Vim_map(2001:end,:);
    Vim_map_all(:,:,k)=Vim_map;
    xlim([2500 3500]);ylim([0 25]);
    title('Vim');ylabel('Number of neurons');
    set(gca,'xtick',[],'FontSize',14);

    tmpVim_all = sqrt(sum(Vim_map,2));
    %     tmp_echt = echt(tmpVim_all, 6, 10, 1e3);
    %     test_amp_ec = abs(tmp_echt);
    %     test_PRC_ec = angle(tmp_echt);
    filt_order = 2;
    [b,a] = butter(filt_order, [4 10]/(1e3/2), 'bandpass');
    Xf = hilbert(filtfilt(b,a,tmpVim_all));
    test_amp_ec = abs(Xf);
    test_PRC_ec = angle(Xf);
    yyaxis right;plot(2001:10000,test_amp_ec,'r','LineWidth',2);set(gca,'YColor','r')

    if k==1
        [~,It]=min(abs(test_amp_ec-0.15));
        plot([2000+It,2000+It],[0,1],'g--','Color',[237 177 32]/255,'LineWidth',2);
    end
    xlim([2500 3500]);ylim([0 0.5]);
    set(gca,'XTick',2500:250:3500,'XTickLabelRotation',0,'FontSize',14);
    title('Vim')
    xlabel('Time (ms)');ylabel('a.u.');

    pos = [251.6667+400*(k-1)   82.3333  419.3333  777.3333];
    tmpf = gcf;tmpf.Position=pos;
    set(gcf,'color','w');
end
