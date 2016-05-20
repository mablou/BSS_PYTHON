clear all
close all

%%
Inline_log = 641.1;
Xline_log = 400.1;

% On charge les donnees de porosite de puits 
Davis_phi_500 = importdata('../data/Davis_NPOR_time_500us');
Phi_log_500 = [Inline_log*ones(length(Davis_phi_500.data),1)...
                     Xline_log*ones(length(Davis_phi_500.data),1)...
                     Davis_phi_500.data(:,1) Davis_phi_500.data(:,2)];


load pdf2d

marg=mean(pdf2d);
marg=marg./sum(marg);

X=Phi_log_500(Phi_log_500(:,3)>=652 & ...
                                        Phi_log_500(:,3)<=686,4);
stddev=sqrt(var(X));

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   ATTENTION!!!!!!!!! FUNCTION GMDISTRIBUTION DOES %%%%%%%%%%%%%%%%%%%
%%%%%   NOT ALWAYS RETURNS SOMETHING GOOD- CHECK INVERSION SUCCESS%%%%%%%%%




GM_model = gmdistribution.fit(X,2);
prop1 = GM_model.PComponents(1);
prop2 = GM_model.PComponents(2);
m1=GM_model.mu(1);
m2=GM_model.mu(2);
std1=prop1*stddev;
std2=prop2*stddev;



% 
% std1=GM_model.Sigma(1);
% std2=GM_model.Sigma(2);

x2=[m1 m2 std1 std2];
%%
figure(2)
grid on
fam1= normpdf(vec_phi,x2(1),x2(3));
fam2= normpdf(vec_phi,x2(2),x2(4));

fam11=fam1/max(fam1)*prop1;
fam22=fam2/max(fam2)*prop2;

% plot(vec_v2,f,'r')
% hold on

hold on

plot(vec_phi,fam11/max(fam11),'g')
hold on
plot(vec_phi,fam22/max(fam11),'k')
hold on
plot(vec_phi,(fam11+fam22)/max(fam11),'r')

% Mean and std deviation of each familiy are saved in stat_poro
save('stat_poro.mat','m1','m2','std1','std2');