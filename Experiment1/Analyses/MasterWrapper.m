%% Basic set-up
restoredefaultpath;
clear all;
addpath('Analyses')
dataDir  = 'Data'; 
condDirs = {'imLeft_presLeft','imLeft_presRight','imRight_presLeft','imRight_presRight'};

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

%% Recode to other - presented
R = CT(:,1); 
R(condition==2 & CT(:,1)==1,1) = 2;
R(condition==2 & CT(:,1)==2,1) = 1;

R(condition==4 & CT(:,1)==1,1) = 2;
R(condition==4 & CT(:,1)==2,1) = 1;

% recode condition to congruency 
cong = ones(size(condition));
cong(ismember(condition,[2,3]))=2;

% R = reality monitoring (RM) response, V(:,2) = averaged vividness over 
% first 9 trials, V(:,3) = vividness last trial, V(:,1) = VVIQ, P(:,1) =
% discriminiation accuracy second task, P(:,2) = discrimination d' and
% P(:,3) = detection d'. Stats RM were done using a binary logistic
% regression with R as dependent variable and cong, V(:,2) and P(:,3) as
% independent variables (see Methods for more details).

%% Exclude participants

% discrimination accuracy
LA = P(:,1)<0.55;
fprintf('Excluded %d pps for for too low accuracy \n',sum(LA))

% imagery check
ima_fail = contains(ima_check,'no','IgnoreCase',1)';
fprintf('Excluded %d pps for failing the imagery_check \n',sum(ima_fail&~LA))

% incorrrect response
IR = R(:,1) == 2;
fprintf('Excluded %d pps for incorrect response \n',sum(IR&~(ima_fail|LA)))

% exclude them
EX = LA | ima_fail | IR;
CT(EX,:) = []; V(EX,:) = []; P(EX,:) = []; condition(EX) = [];
A(EX) = []; CF(EX,:) = []; V2(EX,:,:) = []; R(EX,:) = []; cong(EX) = [];

fprintf('Continuing with %d participants \n',length(A));

%% Plot reality judgements and vividness 
if plotMainFigure
    
PP(1) = sum(R(cong==1,1)==1)/sum(cong==1);
PP(2) = sum(R(cong==2,1)==1)/sum(cong==2);

f = figure; 

subplot(1,4,3:4); b = bar(1,PP(1)); b.FaceColor = 'b'; hold on;
b = bar(2,PP(2)); b.FaceColor = 'r'; 
set(gca,'XTick',1:2);
set(gca,'XTickLabels',{'Congruent','Incongruent'});
ylabel('p(real)'); ylim([0 0.5]); title('Reality monitoring')

% plot avg vividness per condition
VC = zeros(2,2,2);
for c = 1:2
    for r = 1:2
        idx = cong==c & R(:,1)==(r-1);
        VC(c,r,1) = mean(V(idx,3));
        VC(c,r,2) = std(V(idx,3))/sqrt(sum(idx));
    end
end

subplot(1,4,1); 
b = barwitherr(squeeze(VC(1,:,2)),squeeze(VC(1,:,1))); 
b.FaceColor = 'b'; ylim([0 4]); set(gca,'XTickLabels',{'Imagined','Real'});
ylabel('Vividness'); title('Congruent') 

subplot(1,4,2); 
b = barwitherr(squeeze(VC(2,:,2)),squeeze(VC(2,:,1))); 
b.FaceColor = 'r'; ylim([0 4]); set(gca,'XTickLabels',{'Imagined','Real'});
ylabel('Vividness'); title('Incongruent') 

f.Position = [228 258 964 424];

set(gcf,'NumberTitle','off') 
set(gcf,'Name','Figure 2B')

end 

% multinomial regression on vividness
fprintf('Running multinomial regression on vividness \n')
Ord_V = ordinal(V(:,3),{'1','2','3','4','5'},[5,4,3,2,1]);
resp = R(:,1);resp(resp==0)= -1;

[~,~,statsCong] = mnrfit([resp(cong==1)],Ord_V(cong==1),'model','ordinal','link','probit');
[~,~,statsInco] = mnrfit([resp(cong==2)],Ord_V(cong==2),'model','ordinal','link','probit');

fprintf('\t Congruent: t(1) %.3f, p %.3f \n',statsCong.t(5),statsCong.p(5))
fprintf('\t Incongruent: t(1) %.3f,p %.3f \n',statsInco.t(5),statsInco.p(5))


%% Plot vividness over time 
if plotSupplement
    VT = zeros(2,2,10,2);
    for c = 1:2
        for r = 1:2
            for v = 1:10
                idx = cong==c & R(:,1)==(r-1);
                VT(c,r,:,1) = mean(V2(idx,:,1));
                VT(c,r,:,2) = std(V2(idx,:,1))/sqrt(sum(idx));
            end
        end
    end
    
    f = figure;
    subplot(2,1,1);
    errorbar(squeeze(VT(1,1,:,1))',squeeze(VT(1,1,:,2))','b--'); hold on;
    errorbar(squeeze(VT(1,2,:,1))',squeeze(VT(1,2,:,2))','b','LineWidth', 2); ylim([1.5 3.5]); 
    xlim([0.5 10.5]); title('Congruent'); legend('Imagined','Real'); 
    ylabel('Vividness'); xlabel('Trial')
    
    subplot(2,1,2);
    errorbar(squeeze(VT(2,1,:,1))',squeeze(VT(2,1,:,2))','r--'); hold on;
    errorbar(squeeze(VT(2,2,:,1))',squeeze(VT(2,2,:,2))','r','LineWidth', 2); ylim([1.5 3.5]); 
    xlim([0.5 10.5]); title('Incongruent'); legend('Imagined','Real'); 
    ylabel('Vividness'); xlabel('Trial')
    
    f.Position = [340 264 766 452];
    
    set(gcf,'NumberTitle','off')
    set(gcf,'Name','Supplementary Figure 2A-B')
end