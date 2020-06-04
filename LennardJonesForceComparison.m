% Lennard-JOnes potential which is at least roughly realistic for cells

   sigma_max = 15;%2^(-1/6); 
   sigma_max2 = 40; %leave 40, this is for scaling
   dx =1:0.1:200;
   f0 =1;
   n = 2;
   m = 2 * n;
   
   % matrix to store values of potential for different eps
   n_eps = 3; % number of different eps values
    Z = zeros(n_eps,length(dx));
   %eps_ij = 10;
   eps_ij(1)= 3*50*sigma_max/sigma_max2;
   eps_ij(2) =3*150*sigma_max/sigma_max2; 
   eps_ij(3) = 3*200*sigma_max/sigma_max2;


   
   
   z_ij = zeros(1,length(dx));
   z_ijder = zeros(1,length(dx));
   force_ij = zeros(1,length(dx));
   Force = zeros(n_eps,length(dx));
   for j = 1:n_eps
       for i=1:length(dx)
    
           %CiL and CoA
     
        z_ij(i) = eps_ij(j)*((sigma_max/dx(i))^m- (sigma_max/dx(i))^n); 
        z_ijder(i) = n*eps_ij(j)/dx(i)*(2*(sigma_max/dx(i))^m - (sigma_max/dx(i))^n); 
        force_ij(i) = -f0 * z_ijder(i);
       end
       Z(j,:) = z_ij;
       Force(j,:) = force_ij;
   end

%    figure
%    plot(dx,Z(1,:),'LineWidth',3)
%  ylim([-10 10]);
% xlim([0 100])
%  
%  change= [7,9];
%  distance1 = norm(change);
% 
%    title(['m = ' num2str(m)]) % only repulsion
% 
%  % title(['n = ' num2str(n) ', m = ' num2str(m)])
%  xlabel(['Distance between the cells, r (',char(181),'m)'])
%  ylabel('Potential, Z')
%  
%  set(gca,'FontSize',36)
%  ax=gca;
%  box on
% set(gca,'linewidth',4) 
% grid on
%  hold on 
%  
% %  plot(dx, Z(2,:),'--', 'LineWidth',3)
% %  
% %  plot(dx, Z(3,:),':', 'LineWidth',3)
%   
% h4 = plot([sigma_max,sigma_max],[-100,100],'k', 'LineWidth',3)
%  
% %  legend ('\epsilon = 1', '\epsilon = 5', '\epsilon = 10')

 
%  % plot force
  figure 
 
 plot(dx,Force(1,:),'LineWidth',3)
  xlabel('Distance between cells')
 ylabel('Force strength')
 
 hold on
 
 plot(dx,Force(2,:),'--','LineWidth',3)
 
 plot(dx,Force(3,:),':','LineWidth',3)
 
 plot([sigma_max,sigma_max],[-10,10],'k','LineWidth',3)
 
   ylim([-4 4]);
   xlim([0 100])

    set(gca,'FontSize',36)
 ax=gca; 
 set(gca,'linewidth',4) 
 box on
 grid on

 legend ('\epsilon = 19', '\epsilon = 75', '\epsilon = 94')
