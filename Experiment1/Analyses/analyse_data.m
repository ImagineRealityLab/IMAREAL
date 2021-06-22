function [CT,CF,V,P,IC,A,C,V2] = analyse_data(data)

% - imagery vividness - %
V(1) = mean(data.VVIQ);
V(2) = mean(data.imagery(1:end-1,1));
V(3) = data.imagery(end,1); % vividness CT

V2   = data.imagery; % all 10 trls 

% - critical trial - %
CT(1) = data.critical_trial(1,1);
CT(2) = data.critical_trial(2,1);

CF(1) = data.critical_trial(1,3);
CF(2) = data.critical_trial(2,3);

% - age and ima check - %
A = data.age;
IC = contains(data.ima_check,'yes','IgnoreCase',1);
C  = data.comments;

% - discimination performance - % 
P(1) = mean(data.discrimination(:,4));

R    = zeros(length(data.discrimination),1);
R(data.discrimination(:,3)==2 & data.discrimination(:,4)==1) = 2;
R(data.discrimination(:,3)==2 & data.discrimination(:,4)==0) = 1;
R(data.discrimination(:,3)==1 & data.discrimination(:,4)==1) = 1;
R(data.discrimination(:,3)==1 & data.discrimination(:,4)==0) = 2;

% Hits
hits = sum(data.discrimination(:,3)==1 & R==1);
if hits == 0; hits = 0.5; % floor
elseif hits/sum(data.discrimination(:,3)==1)==1 % ceiling
    hits = hits-0.5; 
end
H = hits/sum(data.discrimination(:,3)==1);


% FA's
false_alarms =  sum(data.discrimination(:,3)==2 & R==1);
if false_alarms == 0 % add count if 0
    false_alarms = 0.5;
elseif false_alarms/sum(data.discrimination(:,3)==2)==1
    false_alarms = false_alarms-0.5;
end
FA = false_alarms/sum(data.discrimination(:,3)==2);

P(2) = dprime(H,FA);
P(3) = P(2)/sqrt(2);
