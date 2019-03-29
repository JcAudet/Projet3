

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation de psy initial
% Psy0=[];
% for j = 1 : Jmax
%     for i = 1 : Imax 
%         %La boucle if sert a placer les Conditions Frontieres
%         if i==1
%             Psy0=[Psy0 (Psy_haut(j)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
%         elseif i==Imax
%             Psy0=[Psy0 (Psy_bas(j)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
%         elseif j==1;
%             Psy0=[Psy0 (Psy_gauche(i)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
%         elseif j== Jmax
%             Psy0=[Psy0 (Psy_droite(i)+CI(i*dx,j*dy,x_0,y_0,sig,k))/2];
%         else
%     Psy0=[Psy0 CI(i*dx,j*dy,x_0,y_0,sig,k)];
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% b et f
for j = 1 : Jmax
    for i = 1 : Imax 
       var=(1+1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m)+1i*dt*V_potentiel(i*dx,j*dy)/(2*hbar));
       b_vect=[b_vect var];
       var2=(1-1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m)-1i*dt*V_potentiel(i*dx,j*dy)/(2*hbar));
       f_vect=[f_vect var2];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Re-ecrire les vecteur comme des matrices pour ensuite les ploter
% cell=[];
% s=0;
% for k= 1:Kmax
%     for j=1 : Jmax
%         for i=1 : Imax
%             s=s+1;
%             cell(i,j,k)= Psy(s,k);
%         end
%     end
%     s=0;
% end
 
% for k=1 : Kmax
%     surf(x,y,real(cell(:,:,k)));
%     pause(0.02);
% end





