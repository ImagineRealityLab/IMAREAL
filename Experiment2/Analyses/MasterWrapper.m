%% Set paths etc.
restoredefaultpath;
clear all;
addpath('Analyses')
dataDir  = 'Data'; 
condDirs = {'imLeft_presNothing','imRight_presNothing'};

plotMainFigure = true;
plotSupplement = true;

%% Collect all data
CT = []; P = []; V = []; CF = []; A = []; V2 = [];
condition      = []; prolificID = [];
c = 1;
for cond = 1:length(condDirs)
    
    fprintf('Loading data for condition %d \n',cond)
    
    dataFiles = str2fullfile(fullfile(dataDir,condDirs{cond}),'*.mat');
    nSubs = length(dataFiles);      
    
    critical_trial = zeros(nSubs,2);
    confidence     = zeros(nSubs,2);
    vividness      = zeros(nSubs,3);
    performance    = zeros(nSubs,3);
    age            = zeros(nSubs,1);
    v2             = zeros(nSubs,10,2);
    for sub = 1:nSubs
        
        load(dataFiles{sub},'data');
        
        % reformat to usable output
        [critical_trial(sub,:),confidence(sub,:),vividness(sub,:),...
            performance(sub,:),~,age(sub,1),~,v2(sub,:,:)] = analyse_data(data);
        
        condition(c) = cond;
        comments{c} = data.comments;
        ima_check{c} = data.ima_check;
        c = c+1;
        
        clear data;
        
    end
    
    CT = [CT; critical_trial];
    P  = [P; performance];
    V  = [V; vividness];
    V2 = cat(1,V2,v2);
    CF = [CF; confidence];
    A  = [A; age]; clear age;
    
end

condition = condition';
comments  = comments';


%% Recode to other - imagined 
R = CT(:,1);
R(condition==2 & CT(:,1)==2,1) = 1;
R(condition==2 & CT(:,1)==1,1) = 2;

%% Exclude participants

% discrimination accuracy
LA = P(:,1)<0.55;
fprintf('Excluded %d pps for for too low accuracy \n',sum(LA))

% imagery check
ima_fail = contains(ima_check,'no','IgnoreCase',1)';
fprintf('Excluded %d pps for failing the imagery_check \n',sum(ima_fail))

% incorrrect response
IR = R(:,1) == 2;
fprintf('Excluded %d pps for incorrect response \n',sum(IR|ima_fail))

% exclude them
EX = ima_fail | LA | IR;
CT(EX,:) = []; V(EX,:) = []; P(EX,:) = []; condition(EX) = [];
A(EX) = []; CF(EX,:) = []; V2(EX,:,:) = []; R(EX,:) = []; 

fprintf('Continuing with %d participants \n',length(A));

fprintf('\t PRESENCE RESPONSES: %d \n',sum(R(:,1)==1))

%% Plot reality judgements and vividness 
if plotMainFigure
f = figure; 
subplot(1,2,2);
bar((sum(R==1)/length(R)),'y'); 
ylabel('p(real)'); ylim([0 0.5]); 
title('Reality monitoring'); 
set(gca,'XTickLabels','no input')

% plot avg vividness per condition
VC = zeros(2,2);
for r = 1:2
    idx = R(:,1)==(r-1);
    VC(r,1) = mean(V(idx,3));
    VC(r,2) = std(V(idx,3))/sqrt(sum(idx));
end

subplot(1,2,1); 
b = barwitherr(squeeze(VC(:,2)),squeeze(VC(:,1))); 
b.FaceColor = 'y'; ylim([0 4]); set(gca,'XTickLabels',{'Imagined','Real'});
ylabel('Vividness'); title('No input') 

f.Position = [466 365 521 252];
set(gcf,'NumberTitle','off')
set(gcf,'Name','Figure 2D')
end

% multinomial regression to test 
fprintf('Running multinomial regression on vividness \n')
Ord_V = ordinal(V(:,3),{'1','2','3','4','5'},[5,4,3,2,1]);
resp = R(:,1);resp(resp==0)= -1;

[~,~,stats] = mnrfit(resp,Ord_V,'model','ordinal','link','probit');

fprintf('\t Congruent: t(1) %.3f, p %.3f \n',stats.t(5),stats.p(5))

%% Vividness over time
if plotSupplement
VT = zeros(2,10,2);
for r = 1:2
    for v = 1:10
        idx = R==(r-1);
        VT(r,:,1) = mean(V2(idx,:,1));
        VT(r,:,2) = std(V2(idx,:,1))./sqrt(sum(idx));
    end
end

f = figure;
errorbar(squeeze(VT(1,:,1))',squeeze(VT(1,:,2))','y--'); hold on
errorbar(squeeze(VT(2,:,1))',squeeze(VT(2,:,2))','y', 'LineWidth',2)
xlim([0.5 10.5]); title('No input'); legend('Imagined','Real'); 
ylabel('Vividness'); xlabel('Trial'); ylim([1.5 3.5]);

f.Position = [349 129 786 263];
set(gcf,'NumberTitle','off')
set(gcf,'Name','Supplementary Figure 2C')
end
