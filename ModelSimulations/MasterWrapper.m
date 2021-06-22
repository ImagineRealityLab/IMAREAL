%% Reality monitoring model 
% We model congruency by having two channels, one for each stimulus. We
% hypothesize 3 different possibilities for source monitoring: (H0) complete
% seperation between imagery and perception, so that they don't influence
% each other at all; (H1) Suppressive influence of imagery on perception,
% modelling "I am imagining, so what I see must be imagined"; (H2)
% Facilitating influence of imagery on perception, where external and
% internal input is combined to create one sensory experience of a certain
% vividness

% Perception and imagery are modelled by 2D Gaussians that interact in
% different ways according to the different models. The reality monitoring
% judgement, whether a stimulus is judged to be real, is determined by
% evaluating whether the sample taken from the perception Gaussian crosses 
% the reality monitoring threshold. Vividness is determined by the strength
% of the sample taken fromt the imagery Gaussian. 

clear;

N = 1000; % number of participants
group_vividness = [2.5 0]; % only imagine one stim 
group_vividness_SD = [1 0; 0 1];
threshold = 2.5; % reality threshold 

rng(1); 

Nstim = 2;
trialSD = [0.5 0; 0 0.5];

alpha = 2;

f = figure;

%% Individual differences in imagery vividness
mean_vividness = zeros(N,Nstim);

for i = 1:N
    
    mean_vividness(i,:) = mvnrnd(group_vividness, group_vividness_SD); % maybe should be a different distr where all vals > 0 
    
end

%% H0: No interaction 

% ---- Congruent condition ---- %
V_cong  = zeros(N,1,2);
SJ_cong = nan(N,1);
im_stim = 1; pres_stim = 1; % present same stim
dprime = zeros(2,1); dprime(pres_stim) = 1; 

for i = 1:N
    % Sample stimulus
    X(i,:) = alpha.*mvnrnd(dprime,[1 0; 0 1]);         
    
    % Sample imagery
    V_cong(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD);
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_cong(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_cong(i) = 2; % other
    elseif ~any(X(i,:)>threshold)
        SJ_cong(i) = 0; % none
    end
end

% get vividness imagined stim
V_cong = squeeze(V_cong(:,:,im_stim));

% remove 'other' responses
ridx = SJ_cong == 2;
V_cong(ridx,:) = [];
SJ_cong(ridx)  = [];

% ---- Incongruent condition ----- %
V_inco  = zeros(N,1,2);
SJ_inco = nan(N,1);
im_stim = 1; pres_stim = 2; % present different stim 
dprime = zeros(2,1); dprime(pres_stim) = 1; 

for i = 1:N
    % Sample stimulus
    X(i,:) = alpha.*mvnrnd(dprime,[1 0; 0 1]);            
    
    % Sample imagery
    V_inco(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD);
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_inco(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_inco(i) = 2; % other
    elseif ~any(X(i,:)>threshold)
        SJ_inco(i) = 0; % none
    end
end

% get vividness imagined stim
V_inco = squeeze(V_inco(:,:,im_stim));

% remove 'other' responses
ridx = SJ_inco == 2;
V_inco(ridx,:) = [];
SJ_inco(ridx)  = [];

% --- No stimulus ---- %
V_nost  = zeros(N,1,2);
SJ_nost = nan(N,1);
im_stim = 1; pres_stim = 1; % evaluate same stim
dprime = zeros(2,1); 

for i = 1:N
    % Sample stimulus
    X(i,:) = alpha.*mvnrnd(dprime,[1 0; 0 1]);         
    
    % Sample imagery
    V_nost(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD);
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_nost(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_nost(i) = 2; % other
    elseif ~any(X(i,:)>threshold)
        SJ_nost(i) = 0; % none
    end
end

% get vividness imagined stim
V_nost = squeeze(V_nost(:,:,im_stim));

% remove 'other' responses
ridx = SJ_nost == 2;
V_nost(ridx,:) = [];
SJ_nost(ridx)  = [];

% --- percentage real --- %
pCong = sum(SJ_cong==1)/length(SJ_cong); 
pInco = sum(SJ_inco==1)/length(SJ_inco);
pNost = sum(SJ_nost==1)/length(SJ_nost);
subplot(2,3,1);
b = bar(1,pCong.*100); b(1).FaceColor = 'b'; hold on
b = bar(2,pInco.*100); b(1).FaceColor = 'r'; hold on;
b = bar(3,pNost.*100); b(1).FaceColor = 'y';
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'Congruent','Incongruent','No stimulus'})
ylim([0 100]); ylabel('Percentage real')

title('H0: No interaction');
% --- vividness --- % 
V(1,1,1) = mean(V_cong(SJ_cong==0)); V(1,1,2) = std(V_cong(SJ_cong==0))./sqrt(sum(SJ_cong==0));
V(1,2,1) = mean(V_cong(SJ_cong==1)); V(1,2,2) = std(V_cong(SJ_cong==1))./sqrt(sum(SJ_cong==1));

V(2,1,1) = mean(V_inco(SJ_inco==0)); V(2,1,2) = std(V_inco(SJ_inco==0))./sqrt(sum(SJ_inco==0));
V(2,2,1) = mean(V_inco(SJ_inco==1)); V(2,2,2) = std(V_inco(SJ_inco==1))./sqrt(sum(SJ_inco==1));

V(3,1,1) = mean(V_nost(SJ_nost==0)); V(3,1,2) = std(V_nost(SJ_nost==0))./sqrt(sum(SJ_nost==0));
V(3,2,1) = mean(V_nost(SJ_nost==1)); V(3,2,2) = std(V_nost(SJ_nost==1))./sqrt(sum(SJ_nost==1));

subplot(2,3,4);
barwitherr(squeeze(V(:,:,2)),squeeze(V(:,:,1)));
set(gca,'Xtick',1:3); ylim([0 5])
set(gca,'XTickLabels',{'congruent','incongruent','no stimulus'}); legend('imagined','real');
ylabel('Vividness')
%% H1: Unidirectional suppressive

% ---- Congruent condition ---- %
V_cong  = zeros(N,1,2);
SJ_cong = nan(N,1);
im_stim = 1; pres_stim = 1; % present same stim
dprime = zeros(2,1); dprime(pres_stim) = 1; 

for i = 1:N          
    
    % Sample imagery
    V_cong(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD);
    
    % Sample stimulus - imagery
    X(i,:) = (alpha.*mvnrnd(dprime,[1 0; 0 1])) - squeeze(V_cong(i,1,:))';   
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_cong(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_cong(i) = 2; % other
    elseif ~any(X(i,:)> threshold)
        SJ_cong(i) = 0; % none
    end
end

% get vividness imagined stim
V_cong = squeeze(V_cong(:,:,im_stim));

% remove 'other' responses
ridx = SJ_cong == 2;
V_cong(ridx,:) = [];
SJ_cong(ridx)  = [];

% ---- Incongruent condition ----- %
V_inco  = zeros(N,1,2);
SJ_inco = nan(N,1);
im_stim = 1; pres_stim = 2; % present different stim 
dprime = zeros(2,1); dprime(pres_stim) = 1; 

for i = 1:N
    
    % Sample imagery
    V_inco(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD);
    
    % Sample stimulus - imagery
    X(i,:) = (alpha.*mvnrnd(dprime,[1 0; 0 1])) - squeeze(V_inco(i,1,:))';  
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_inco(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_inco(i) = 2; % other
    elseif ~any(X(i,:)>threshold)
        SJ_inco(i) = 0; % none
    end
end

% get vividness imagined stim
V_inco = squeeze(V_inco(:,:,im_stim));

% remove 'other' responses
ridx = SJ_inco == 2;
V_inco(ridx,:) = [];
SJ_inco(ridx)  = [];


% ---- No stimulus ----- %
V_nost  = zeros(N,1,2);
SJ_nost = nan(N,1);
im_stim = 1; pres_stim = 1; % evaluate same stim 
dprime = zeros(2,1); 

for i = 1:N
    
    % Sample imagery
    V_nost(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD);
    
    % Sample stimulus - imagery
    X(i,:) = (alpha.*mvnrnd(dprime,[1 0; 0 1])) - squeeze(V_nost(i,1,:))';  
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_nost(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_nost(i) = 2; % other
    elseif ~any(X(i,:)>threshold)
        SJ_nost(i) = 0; % none
    end
end

% get vividness imagined stim
V_nost = squeeze(V_nost(:,:,im_stim));

% remove 'other' responses
ridx = SJ_nost == 2;
V_nost(ridx,:) = [];
SJ_nost(ridx)  = [];

% --- percentage imagined --- %
pCong = sum(SJ_cong==1)/length(SJ_cong); 
pInco = sum(SJ_inco==1)/length(SJ_inco);
pNost = sum(SJ_nost==1)/length(SJ_nost);
subplot(2,3,2);
b = bar(1,pCong.*100); b(1).FaceColor = 'b'; hold on
b = bar(2,pInco.*100); b(1).FaceColor = 'r'; hold on;
b = bar(3,pNost.*100); b(1).FaceColor = 'y';
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'Congruent','Incongruent','No stimulus'})
ylim([0 100]); ylabel('Percentage real')
title('H1: Unidirectional suppression');

% --- vividness --- % 
V(1,1,1) = mean(V_cong(SJ_cong==0)); V(1,1,2) = std(V_cong(SJ_cong==0))./sqrt(sum(SJ_cong==0));
V(1,2,1) = mean(V_cong(SJ_cong==1)); V(1,2,2) = std(V_cong(SJ_cong==1))./sqrt(sum(SJ_cong==1));

V(2,1,1) = mean(V_inco(SJ_inco==0)); V(2,1,2) = std(V_inco(SJ_inco==0))./sqrt(sum(SJ_inco==0));
V(2,2,1) = mean(V_inco(SJ_inco==1)); V(2,2,2) = std(V_inco(SJ_inco==1))./sqrt(sum(SJ_inco==1));

V(3,1,1) = mean(V_nost(SJ_nost==0)); V(3,1,2) = std(V_nost(SJ_nost==0))./sqrt(sum(SJ_nost==0));
V(3,2,1) = mean(V_nost(SJ_nost==1)); V(3,2,2) = std(V_nost(SJ_nost==1))./sqrt(sum(SJ_nost==1));

subplot(2,3,5);
barwitherr(squeeze(V(:,:,2)),squeeze(V(:,:,1)));
set(gca,'Xtick',1:3); ylim([0 5])
set(gca,'XTickLabels',{'congruent','incongruent','no stimulus'}); legend('imagined','real');
ylabel('Vividness')
%% H2: Complete mixing bidirectional 

% ---- Congruent condition ---- %
V_cong  = zeros(N,1,2);
SJ_cong = nan(N,1);
im_stim = 1; pres_stim = 1; % present same stim
dprime = zeros(2,1); dprime(pres_stim) = 1; 

for i = 1:N          
    
    % Sample imagery + stimulus
    V_cong(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD) + (alpha.*mvnrnd(dprime,[1 0; 0 1]));
    
    % Decide based on total sensory strength
    X(i,:) = squeeze(V_cong(i,1,:))';   
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_cong(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_cong(i) = 2; % other
    elseif ~any(X(i,:)> threshold)
        SJ_cong(i) = 0; % none
    end
end

% get vividness imagined stim
V_cong = squeeze(V_cong(:,:,im_stim));

% remove 'other' responses
ridx = SJ_cong == 2;
V_cong(ridx,:) = [];
SJ_cong(ridx)  = [];


% ---- Incongruent condition ----- %
V_inco  = zeros(N,1,2);
SJ_inco = nan(N,1);
im_stim = 1; pres_stim = 2; % present different stim 
dprime = zeros(2,1); dprime(pres_stim) = 1; 

for i = 1:N
    
    % Sample imagery + stimulus
    V_inco(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD) + (alpha.*mvnrnd(dprime,[1 0; 0 1]));
    
    % Decide based on total sensory strength
    X(i,:) = squeeze(V_inco(i,1,:))';  
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_inco(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_inco(i) = 2; % other
    elseif ~any(X(i,:)> threshold)
        SJ_inco(i) = 0; % none is above threshold - so lower viv in general 
    end
end

% get vividness imagined stim
V_inco = squeeze(V_inco(:,:,im_stim));

% remove 'other' responses
ridx = SJ_inco == 2;
V_inco(ridx,:) = [];
SJ_inco(ridx)  = [];


% ---- No stimulus ----- %
V_nost  = zeros(N,1,2);
SJ_nost = nan(N,1);
im_stim = 1; pres_stim = 1; % evaluate imagined stim 
dprime = zeros(2,1); 

for i = 1:N
    
    % Sample imagery + stimulus
    V_nost(i,1,:) = mvnrnd(mean_vividness(i,:), trialSD) + (alpha.*mvnrnd(dprime,[1 0; 0 1]));
    
    % Decide based on total sensory strength
    X(i,:) = squeeze(V_nost(i,1,:));%(alpha.*mvnrnd(dprime,[1 0; 0 1])) + squeeze(V_nost(i,1,:))';  
    
    % source judgement 
    if X(i,pres_stim) > threshold
        SJ_nost(i) = 1; % presented
    elseif X(i,3-pres_stim) > threshold
        SJ_nost(i) = 2; % other
    elseif ~any(X(i,:)>threshold)
        SJ_nost(i) = 0; % none is above threshold - so lower viv in general 
    end
end

% get vividness imagined stim
V_nost = squeeze(V_nost(:,:,im_stim));

% remove 'other' responses
ridx = SJ_nost == 2;
V_nost(ridx,:) = [];
SJ_nost(ridx)  = [];


% --- percentage real --- %
pCong = sum(SJ_cong==1)/length(SJ_cong); 
pInco = sum(SJ_inco==1)/length(SJ_inco);
pNost = sum(SJ_nost==1)/length(SJ_nost);
subplot(2,3,3);
b = bar(1,pCong.*100); b(1).FaceColor = 'b'; hold on
b = bar(2,pInco.*100); b(1).FaceColor = 'r'; hold on;
b = bar(3,pNost.*100); b(1).FaceColor = 'y';
set(gca,'XTick',1:3); set(gca,'XTickLabels',{'Congruent','Incongruent','No stimulus'})
ylim([0 100]); ylabel('Percentage real')
title('H2: Complete mixing');


% --- vividness --- % 
V(1,1,1) = mean(V_cong(SJ_cong==0)); V(1,1,2) = std(V_cong(SJ_cong==0))./sqrt(sum(SJ_cong==0));
V(1,2,1) = mean(V_cong(SJ_cong==1)); V(1,2,2) = std(V_cong(SJ_cong==1))./sqrt(sum(SJ_cong==1));

V(2,1,1) = mean(V_inco(SJ_inco==0)); V(2,1,2) = std(V_inco(SJ_inco==0))./sqrt(sum(SJ_inco==0));
V(2,2,1) = mean(V_inco(SJ_inco==1)); V(2,2,2) = std(V_inco(SJ_inco==1))./sqrt(sum(SJ_inco==1));

V(3,1,1) = mean(V_nost(SJ_nost==0)); V(3,1,2) = std(V_nost(SJ_nost==0))./sqrt(sum(SJ_nost==0));
V(3,2,1) = mean(V_nost(SJ_nost==1)); V(3,2,2) = std(V_nost(SJ_nost==1))./sqrt(sum(SJ_nost==1));

subplot(2,3,6);
barwitherr(squeeze(V(:,:,2)),squeeze(V(:,:,1)));
set(gca,'Xtick',1:3); ylim([0 5])
set(gca,'XTickLabels',{'congruent','incongruent','no stimulus'}); legend('imagined','real');
ylabel('Vividness')

f.Position = [1 41 1440 783];