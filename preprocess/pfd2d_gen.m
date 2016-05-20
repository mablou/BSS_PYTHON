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

%% PLOT 2D
%scatter plot of both variables to see correlation
   
plot(PHI_TOT ,IA_TOT,'o');

%%
% %% Bandwidth definition for kernel
% 
% Bandwitdh calculation Silverman (1986)
 l_ia1=0.9 * (min(std(IA_TOT),iqr_mb(IA_TOT)/1.349)) * length(IA_TOT)^(-1/5);
 l_phi1=0.9 * (min(std(PHI_TOT),iqr_mb(PHI_TOT)/1.349)) * length(PHI_TOT)^(-1/5);
% 
% % Bowman et Foster
 l_ia2=std(IA_TOT)*length(IA_TOT)^(-1/6);
 l_phi2=std(PHI_TOT)*length(PHI_TOT)^(-1/6);
% 
% % Regle de pouce (Deheurels (1977))
 l_ia3=1.06*std(IA_TOT)* length(IA_TOT)^(-1/5);
 l_phi3=1.06*std(PHI_TOT)* length(PHI_TOT)^(-1/5);
% 
   disp('l_ia')
   disp(l_ia1)
   disp(l_ia2)
   disp(l_ia3)
%    
   l_ia=input('Enter value for IP bandwidth : ');
   
   disp('l_phi')
   disp(l_phi1)
   disp(l_phi2)
   disp(l_phi3)
%    
   l_phi=input('Enter value for PHI bandwidth : ');
% 
 %% Preparation de la grille pour le kernel
% 
 
% % min_ia=(round(min(cube_ai))); %donnees provenant du cube AI et non des forages
% % max_ia=(round(max(cube_ai)));

% Choose boundaries and grid spacing for kernel grid

min_phi = 0.12;
max_phi = 0.3;
min_IP = 9800;
max_IP = 18500;
N_grid = 200;

% The function grille2 just prepare a grid to calculate the kernel on
gril2d=grille2(min_phi,max_phi,(max_phi-min_phi)/(N_grid-1),...
                    min_IP,max_IP,(max_IP-min_IP)/(N_grid-1));

 l=length(gril2d);
%% CALCUL DU KERNEL
% 
 r=zeros(l,1);
% 
for q=1:length(PHI_TOT)
    
        g=(1/ l_ia)*(1/ l_phi)*(1/sqrt(2*pi)).*exp(-0.5.*...
                    (((IA_TOT(q)-gril2d(:,2))./ l_ia).^2)).*...
                                     (1/sqrt(2*pi)).*exp(-0.5.*...
                           (((PHI_TOT(q)-gril2d(:,1))./l_phi).^2));
 
        r=r+g;
end

pdf22=(1/length(IA_TOT)).*r;
 
%% RESHAPE DU PDF 2D EN PLAN
 
% % Definition of vector for each variable with regular steps identical to
% the kernel
    vec_phi=(min_phi:(max_phi-min_phi)/(N_grid-1):max_phi)';
    vec_ia=(min_IP:(max_IP-min_IP)/(N_grid-1):max_IP)';
    
 pdf2d=reshape(pdf22,length(vec_phi),length(vec_ia));
 pdf2d=pdf2d';
 
 %%
 
% Plot of the kernel == joint prob function
% 
figure(3)
imagesc(vec_phi,vec_ia,pdf2d)
colorbar
ylabel('IP','fontsize',18)
xlabel('Porosity','fontsize',18)
% caxis([0 1e-3])
set(gca,'YDir','normal','Fontname','Helvetica','fontsize',14)

%outputs to be used by the main algorithm routine
save('pdf2d','pdf2d','vec_ia','vec_phi');