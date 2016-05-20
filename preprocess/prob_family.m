clear all
close all
    

%%
% Loading well data
Davis_phi_low_def= importdata('../data/Davis_NPOR_low_def');
Davis_IP_low_def = importdata('../data/Davis_IP_low_def');

%Restrict well data to the MF formation                 
PHI_TOT = Davis_phi_low_def.data(Davis_phi_low_def.data(:,1)>=652 & ...
                                        Davis_phi_low_def.data(:,1)<=686,2);
IA_TOT = Davis_IP_low_def.data(Davis_IP_low_def.data(:,1)>=652 & ...
                                        Davis_IP_low_def.data(:,1)<=686,2);
          
%%
                 
load pdf2d % Kernel 2D savec previously
figure
imagesc(vec_phi,vec_ia,pdf2d)
colorbar
ylabel('IP')
xlabel('Porosity')
% caxis([0 9e-6])
set(gca,'YDir','normal','Fontname','Helvetica','fontsize',20)

hold on
plot(PHI_TOT,IA_TOT,'w+')
% printFigureToPdf([cd,'/img','pdf2d.pdf'], [12,8],'in');
% clf
clear vec_ia vec_phi;

%% Classify data in family
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOU NEED TO START THIS WHOLE CELL IN ONE 
%   EXECUTION AND FOLLOW THE INSTRUCTIONS
%   IN THE MATLAB COMMAND WINDOW
%   PRESS ENTER WHEN YOU SEE FAMILY ONE 
%   TO START DEFINING THE POINTS POLYGON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


matrix_tot = [PHI_TOT,IA_TOT];
sy = ['y*';'k.';'mo'];% add symbols for more families
IDX=zeros(length(matrix_tot),2);
for i = 1:2 % add a greater number than 2 for more families

    % pick values in each class
    input(['FAMILY ',int2str(i)])
    figure(1)
    x = ginput(20);
    x(end+1,:) = x(1,:);
    in = inpolygon(matrix_tot(:,1),matrix_tot(:,2),x(:,1),x(:,2));
    IDX(:,i) = in;
    plot(matrix_tot(in,1),matrix_tot(in,2),sy(i,:),'markersize',15)

    IA_TOT = matrix_tot(in,2);       % IP data in family
    PHI_family = matrix_tot(in,1);     % Phi data in family

    %% Definition des largeurs de bande par 3 differents auteurs
    % 
    % % Largeur de bande proposee par Silverman (1986)
    l_ia1=0.9 * (min(std(IA_TOT),iqr(IA_TOT)/1.349)) * length(IA_TOT)^(-1/5);
    l_phi1=0.9 * (min(std(PHI_family),iqr(PHI_family)/1.349)) * length(PHI_family)^(-1/5);
    %
    % % par Bowman et Foster
    l_ia2=std(IA_TOT)*length(IA_TOT)^(-1/6);
    l_phi2=std(PHI_family)*length(PHI_family)^(-1/6);
    % 
    % % par regle de pouce (Deheurels (1977))
    l_ia3=1.06*std(IA_TOT)* length(IA_TOT)^(-1/5);
    l_phi3=1.06*std(PHI_family)* length(PHI_family)^(-1/5);
    % 
    disp('l_ia')
    disp(l_ia1)
    disp(l_ia2)
    disp(l_ia3)
    l_ia=input('Entrez la valeur de l_ia : ');

    disp('l_phi')
    disp(l_phi1)
    disp(l_phi2)
    disp(l_phi3)
    l_phi=input('Entrez la valeur de l_phi : ');


  %% Preparation de la grille pour le kernel
% 
 min_ia=(round(min(IA_TOT)/100))*100;
 max_ia=(round(max(IA_TOT)/100))*100;
 
% % min_ia=(round(min(cube_ai))); %donnees provenant du cube AI et non des forages
% % max_ia=(round(max(cube_ai)));
% 
 min_phi=round(min(PHI_TOT)*100)/100;
 max_phi=round(max(PHI_TOT)*100)/100;
 

gril2d=grille2(min_phi,max_phi,0.001,7000,20000,10);


 l=length(gril2d);

    %% CALCUL DU KERNEL (sur donnes par facies)
    %
    r=zeros(l,1);
    %
    for q=1:length(IA_TOT)
            g=(1/ l_ia)*(1/ l_phi)*(1/sqrt(2*pi)).*exp(-0.5.*...
                (((IA_TOT(q)-gril2d(:,2))./ l_ia).^2)).*(1/sqrt(2*pi)).*...
                        exp(-0.5.*(((PHI_family(q)-gril2d(:,1))./l_phi).^2));
            r=r+g;
    end
    pdf22=(1/length(IA_TOT)).*r;

    %% RESHAPE DU PDF 2D EN PLAN

    % % Definition des vecteurs pour trouver les indices dans simul_bayes
     vec_ia=(7000:10:20000)';
     vec_phi=(min_phi:0.001:max_phi)';
      
     
     if i ==1
         vec_ia1=vec_ia;
         vec_phi1=vec_phi;
     elseif i==2
         vec_ia2=vec_ia;
         vec_phi2=vec_phi;
         
     end
 pdf2d=reshape(pdf22,length(vec_phi),length(vec_ia));
 pdf2d=pdf2d';
    PDF2D(:,:,i) = pdf2d;

    %% AFFICHAGE 3D DU KERNEL AVEC DES SLICES 
    %
    figure(i+1)
    imagesc(vec_phi,vec_ia,pdf2d)
    colorbar
    ylabel('IP')
    xlabel('PHI')
    caxis([0 max(max(pdf2d))])
    set(gca,'YDir','normal','Fontname','Helvetica','fontsize',14)
%    hold on
%    plot(matrix_tot(in,4),matrix_tot(in,5),'wo')
%      printFigureToPdf([cd,'/fig','kernel_slice',int2str(i),'.pdf'], [12,8],'in');

end
s = squeeze(sum(PDF2D,2));
sk = sum(s,2);
sk = sk*ones(1,2);
prob = s./sk;
sumprob = sum(prob,2);

figure
hold on
plot(vec_ia,prob(:,1),'r-')
plot(vec_ia,prob(:,2),'b-')
ylabel('Probability','FontSize',14)
xlabel('IP','FontSize',14)
axis([8e03,2e04,0,1])
legend({'Famille 1','Famille 2'},'FontSize',14)

% printFigureToPdf([cd,'/fig','prob.pdf'], [6,4],'in');

%dfdsa
save('prob_2facies_phi','prob','IDX','PDF2D');
prob1=prob(:,1);
prob2=prob(:,2);
save('prob','prob1','prob2');
