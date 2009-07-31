%TESTFITGADAMODEL Test the GADA-JRN and GADA-SMN algorithms
clear all;
disp('# GADA - Genome Alteration Detection Analysis');
disp('# Copyright (C) 2008,2009 Childrens Hospital Los Angeles');
disp('# Author: Roger Pique-Regi rpique@gmail.com');


%% Simulation Data generation parameters
M=1000;        % Size of the chromosome, number of probes on the array. 
N=20;           % Number of array samples.
sigma2=0.2;     % Noise Variance
IntLength=20;   % CNVR1 length
BrLength=0;     % CNVR1 is aligned across multiple samples BrLength=0
BrLength2=20;   % CNVR2 breakpoints are randomly distributed U(-20,20)
IntLength2=100; % CNVR2 length

%% Generate Simlation Data. 
Xmat=zeros(N,M);  % Xmat is the underlying copy number
Ymat=Xmat;        % Ymat is the simulated observed data with microarray
sigma2vec=zeros(N,1);
Ir=cell(1,N);

for n=1:N/2  %First CNVR is roughly placed at M/4
    xini = (round(rand*BrLength)+M/4-IntLength/2-BrLength);
    xend = (round(rand*BrLength)+M/4+IntLength/2);
    Xmat(n,xini:xend) = -1;
    sigma2vec(n) = sigma2;
    Ymat(n,:) = Xmat(n,:)+randn(1,M)*sqrt(sigma2vec(n));
end

for n=n+1:N %Second CNVR is roughly placed at 3*M/4
    xini = (round(rand*BrLength2)+3*M/4-IntLength2/2-BrLength2);
    xend = (round(rand*BrLength2)+3*M/4+IntLength2/2);
    Xmat(n,xini:xend) = 1; %0.58;
    sigma2vec(n) = sigma2;
    Ymat(n,:) = Xmat(n,:)+randn(1,M)*sqrt(sigma2vec(n));
end

% Generating the probe hybridization bias. 
TrueRef=0.5*randn(1,M)+0.5*sin(2*pi*(1:M)/M*10); %Wave spatially uncorrelated + correlated  components

%TrueRho=zeros(N,1); %% No wave bias is added 
%TrueRho=(rand(N,1)*0.6+0.7);  %% Wave with different scale on each sample is added. 
TrueRho=ones(N,1);  %% Wave with same scale on each sample is added. 

for n=1:N
     Ymat(n,:)=Ymat(n,:)+TrueRho(n)*TrueRef; 
end

%% Parameter settings
aAlpha=0.5;
TBackElim=4;
MinSegLen=3;
verbose=2; %Debugging messages.

CommonParameters={...
    'TBackElim',TBackElim,'MinSegLen',MinSegLen,...
    'verbose',verbose};

%  isRemovingOK=0;
%  CommonParameters={...
%      'TBackElim',TBackElim,'MinSegLen',MinSegLen,...
%      'isRemovingOK',isRemovingOK,...
%      'AdjustSampleSigma2',0,'SampleSigma2',sigma2vec...
%      'verbose',verbose};

%% Placing input data in a binary file for memmory mapping

fd=fopen('mytemp.tmp','wb');
fwrite(fd,Ymat','double');
fclose(fd);
Ymat2=memmapfile('mytemp.tmp','format',{'double',[M N], 'x'});

%%
[OldGADA.Iext,OldGADA.Wext,OldGADA.par]=fitGadaModel(Ymat2,...
    'aAlpha',aAlpha,'RefAdjust',0,'RhoAdjust',0,...
    CommonParameters{:});

%% 
[GADA_JRN.Iext,GADA_JRN.Wext,GADA_JRN.par]=fitGadaModel(Ymat2,...
    'aAlpha',aAlpha,'RefAdjust',1,'RhoAdjust',0,...
    CommonParameters{:});

%%
[GADA_SMN.Iext,GADA_SMN.Wext,GADA_SMN.par]=fitGadaModel(Ymat2,...
    'aAlpha',aAlpha,'RefAdjust',2,'RhoAdjust',0,...
    CommonParameters{:});

%%

clear Ymat2

%% Final results Plotting

figure(1);clf;

% Original
subplot(2,3,1);
imagesc(Xmat,[-2.0 2.0]);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('Underlying CNA component $x$','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'A','fontsize',14);

% Observed
%figure(5);
subplot(2,3,2);
imagesc(Ymat,[-2.0 2.0]);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('Observed array data $y_{mn} = x_{mn} + \rho_n r_m + \epsilon_{mn}$','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'B','fontsize',14);

subplot(2,3,4);
YmatN=zeros(N,M);
for n=1:N
    YmatN(n,:)=Ymat(n,:)-GADA_SMN.par.Ref;
end
imagesc(YmatN,[-2.0 2.0]);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('Median normalized data: $\tilde y = y -$ median $(y)$','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'D','fontsize',14);

subplot(2,3,6);
Xrec=zeros(N,M);
for n=1:N
    Xrec(n,:)=ReconstructFromIextWext(GADA_JRN.Iext{n},GADA_JRN.Wext{n});
end
imagesc(Xrec,[-2.0 2.0]);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('GADA with joint reference normalization (GADA-JRN)','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'F','fontsize',14);


subplot(2,3,3);
Xrec2=zeros(N,M);
for n=1:N
    Xrec2(n,:)=ReconstructFromIextWext(OldGADA.Iext{n},OldGADA.Wext{n});
end
imagesc(Xrec2,[-2.0 2.0]);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('GADA without normalization','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'C','fontsize',14);

subplot(2,3,5);
Xrec3=zeros(N,M);
for n=1:N
    Xrec3(n,:)=ReconstructFromIextWext(GADA_SMN.Iext{n},GADA_SMN.Wext{n});
end
imagesc(Xrec3,[-2.0 2.0]);
xlabel('Probe $m$','Interpreter','latex','fontsize',10);
ylabel('Sample $n$','Interpreter','latex','fontsize',10);
title('GADA with separate median normalization (GADA-SMN)','Interpreter','latex','fontsize',11);
text(-M*0.15,-1,'E','fontsize',14);

%%
figure1=gcf
% Create rectangle
annotation(figure1,'rectangle',[0.09564 0.0709 0.5478 0.4212],...
    'FaceColor','flat');

% Create arrow
annotation(figure1,'arrow',[0.406 0.3488],[0.5694 0.4677]);

% Create arrow
annotation(figure1,'arrow',[0.3548 0.3821],[0.272 0.272]);

% Create arrow
annotation(figure1,'arrow',[0.6345 0.6583],[0.7429 0.7439]);

% Create arrow
annotation(figure1,'arrow',[0.6405 0.6875],[0.5662 0.4899]);

% Create ellipse
annotation(figure1,'ellipse',[0.1696 0.7153 0.0266 0.218]);

% Create ellipse
annotation(figure1,'ellipse',[0.2765 0.5726 0.0266 0.218]);

% Create textbox
annotation(figure1,'textbox',[0.1595 0.6837 0.05357 0.0328],...
    'String',{'CNVR-1'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Arial',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',[0.2664 0.7843 0.05357 0.0328],...
    'String',{'CNVR-2'},...
    'FontWeight','bold',...
    'FontSize',12,...
    'FontName','Arial',...
    'LineStyle','none');

delete *.tmp

