%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lam

Lam=[lamda lamda*2 lamda*4 lamda*8];

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;

% Init matrices
norm_l=zeros(length(Lam),length(t));
norm_th_l=zeros(length(Lam),length(t));
err_l=zeros(length(Lam),length(t));

% M et M2
v_mat=zeros(length(y),length(x));
V=v_mat(:);

b = 1 + 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) + 1i*dt*V/(2*hbar);
f = 1 - 1i*hbar*dt*(1/dx^2 + 1/dy^2)/(2*m) - 1i*dt*V/(2*hbar);

c=-1i*hbar*dt*(1/dx^2)/(4*m);
d=-1i*hbar*dt*(1/dy^2)/(4*m);

g=1i*hbar*dt*(1/dx^2)/(4*m);
k=1i*hbar*dt*(1/dy^2)/(4*m);

diag=ones(1,length(V)-1)*d;
diag1=ones(1,length(V)-1)*d;
diag(ly:ly:length(V)-1)=0;
diag1(ly:ly:length(V)-1)=0;

diag2=ones(1,length(V)-1)*k;
diag3=ones(1,length(V)-1)*k;
diag2(ly:ly:length(V)-1)=0;
diag3(ly:ly:length(V)-1)=0;

M=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[b' ones(1,length(V)-ly)*c ones(1,length(V)-ly)*c diag diag1]);
M2=sparse([1:length(V) ly+1:length(V) 1:length(V)-ly (2:length(V))-1 (1:length(V)-1)+1],...
    [1:length(V) (ly+1:length(V))-ly (1:length(V)-ly)+ly 2:length(V) 1:length(V)-1]...
    ,[f' ones(1,length(V)-ly)*g ones(1,length(V)-ly)*g diag2 diag3]);

for w=1:length(Lam)

    kp=2*(pi/Lam(w))*[1,0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONDITIONS INITIALES
    
    [psy,norm_l(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);

    Psy_mat=zeros(length(y),length(x),length(t));
    Psy_mat(:,:,1)=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
    % v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') );
    % v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
    V=v_mat(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Init theorie
    
    Psy_th=zeros(length(t),length(y),length(x));
    
    [Psy_th(1,:,:),norm_th_l(w,1)]=analy(x,y,t(1));
    
    norm_th_l(w,1)=trapeze_2D(abs(squeeze(Psy_th(1,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    err_l(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

%     figure()
    tic
    for j = 1 : length(t)-1
        
        b=M2*Psy(:,j);
        Psy(:,j+1) = mldivide(M,b);
        
        psy=vec2mat(Psy(:,j+1),length(y))';
        Psy_mat(:,:,j+1)=psy;
        
        [Psy_th(j+1,:,:),norm_th_l(w,j+1)]=analy(x,y,t(j+1));
        
        norm_l(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat(:,:,j+1))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        err_l(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-squeeze(abs(Psy_mat(:,:,j+1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        
    end
    toc
    w    
end

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_l(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la resolution en fonction de la longueur d''onde')

figure()
hold on
for i=1:length(Dx)
    plot(t,norm_th_l(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Dx(1)),sprintf('dx = dy = %d',Dx(2)),sprintf('dx = %d',Dx(3)))
title('Evolution de la norme de la solution analyique en fonction de la longueur d''onde')

figure()
hold on
for i=1:length(Lam)
    plot(t,err_l(i,:))
end
hold off
legend(sprintf('dx = dy = %d',Lam(1)),sprintf('dx = dy = %d',Lam(2)),sprintf('dx = %d',Lam(3)),sprintf('dx = %d',Lam(4)))
title('Evolution de l''erreur avec differentes valeurs de longueur d''onde')

P_l=polyfit(log10(Lam),log10(err_t(:,end)'),1)

figure()
loglog(Lam,err_l(:,end),'ro')
hold on
loglog(Lam,10.^(P_l(1)*log10(Dt) +P_l(2)),'b','LineWidth',1.5)
xlabel('Log(dt)')
ylabel('log(Err)')
legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dt+%.3f',P_l(1),P_l(2)),'Location','Best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sigma

Sig=[sig/2 sig sig*2 sig*4];

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;

% Init matrices
norm_s=zeros(length(Sig),length(t));
norm_th_s=zeros(length(Sig),length(t));
err_s=zeros(length(Sig),length(t));


for w=1:length(Sig)

    sig=Sig(w);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONDITIONS INITIALES
    
    [psy,norm_s(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);

    Psy_mat=zeros(length(y),length(x),length(t));
    Psy_mat(:,:,1)=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Init theorie
    
    Psy_th=zeros(length(t),length(y),length(x));
    
    [Psy_th(1,:,:),norm_th_s(w,1)]=analy(x,y,t(1));
    
    norm_th_s(w,1)=trapeze_2D(abs(squeeze(Psy_th(1,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    err_s(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

%     figure()
    tic
    for j = 1 : length(t)-1
        
        b=M2*Psy(:,j);
        Psy(:,j+1) = mldivide(M,b);
        
        psy=vec2mat(Psy(:,j+1),length(y))';
        Psy_mat(:,:,j+1)=psy;
        
        [Psy_th(j+1,:,:),norm_th_s(w,j+1)]=analy(x,y,t(j+1));
        
        norm_s(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat(:,:,j+1))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        err_s(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-squeeze(abs(Psy_mat(:,:,j+1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        
%         subplot(1,2,1)
%         surf(x,y,abs(squeeze(Psy_mat(j,:,:))).^2,'edgecolor','none');
%         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_x(w,j)))
%         view(0,90);
%         subplot(1,2,2)
%         surf(x,y,abs(squeeze(Psy_th(j,:,:))).^2,'edgecolor','none');
%         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_th_x(w,j)))
%         view(0,90);
%         hold off
%         
%         pause(0.01)
        
    end
    toc
    w    
end

figure()
hold on
for i=1:length(Sig)
    plot(t,norm_s(i,:))
end
hold off
legend(sprintf('Sigma = %d',Sig(1)),sprintf('Sigma = %d',Sig(2)),sprintf('Sigma = %d',Sig(3)))
title('Evolution de la norme de la resolution en fonction de l''ecart type')

figure()
hold on
for i=1:length(Sig)
    plot(t,norm_th_s(i,:))
end
hold off
legend(sprintf('Sigma = %d',Sig(1)),sprintf('Sigma = %d',Sig(2)),sprintf('Sigma = %d',Sig(3)))
title('Evolution de la norme de la solution analyique en fonction de l''ecart type')

figure()
hold on
for i=1:length(Sig)
    plot(t,err_s(i,:))
end
hold off
legend(sprintf('Sigma = %d',Lam(1)),sprintf('Sigma = %d',Lam(2)),sprintf('Sigma = %d',Lam(3)),sprintf('dx = %d',Lam(4)))
title('Evolution de l''erreur avec differentes valeurs de longueur d''onde')

P_s=polyfit(log10(Lam),log10(err_t(:,end)'),1)

figure()
loglog(Lam,err_l(:,end),'ro')
hold on
loglog(Lam,10.^(P_s(1)*log10(Dt) +P_s(2)),'b','LineWidth',1.5)
xlabel('Log(sig)')
ylabel('log(Err)')
legend('Erreur après un temps fixe de propagation',sprintf('Régression y=%.3f*dt+%.3f',P_s(1),P_s(2)),'Location','Best')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DY

Dy=[dy dy*2 dy*4];
y1=0:Dy(end):yf;

x=0:dx:xf;
y=0:dy:yf; ly=length(y);
t=0:dt:tf;

% Init matrices
norm_y=zeros(length(Dy),length(t));
norm_th_y=zeros(length(Dy),length(t));
err_y=zeros(length(Dy),length(t));

for w=1:length(Dy)

    y=0:Dy(w):yf;
    
    ly=length(y);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CONDITIONS INITIALES
    
    [psy,norm_y(w,1)] = wp_ini_2D(x,y,sig,sig,kp,x0,y0);

    Psy=[];
    Psy(:,1)=psy(:);

    Psy_mat=zeros(length(y),length(x),length(t));
    Psy_mat(:,:,1)=psy;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul des facteurs

    v_mat=zeros(length(y),length(x));
    % v_mat=( barr(x,y,1000,x(floor(length(x)/2)),4e-11,y(floor(length(y)/2)),15e-11,10e-11,'Carre') );
    % v_mat=sparse( ( barr_simple(x,y,1000,x(length(x)/2),1e-14,y(length(y)/2),8e-14,'Carre') )' );
    V=v_mat(:);

    b = 1 + 1i*hbar*dt*(1/dx^2 + 1/Dy(w)^2)/(2*m) + 1i*dt*V/(2*hbar);
    f = 1 - 1i*hbar*dt*(1/dx^2 + 1/Dy(w)^2)/(2*m) - 1i*dt*V/(2*hbar);

    c=-1i*hbar*dt*(1/dx^2)/(4*m);
    d=-1i*hbar*dt*(1/Dy(w)^2)/(4*m);

    g=1i*hbar*dt*(1/dx^2)/(4*m);
    k=1i*hbar*dt*(1/Dy(w)^2)/(4*m);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Creation de M et M2, v et v2

    %M=sparse(ones(length(b_vect),1),ones(length(b_vect),1),b_vect);
    M=diag(b);

    %M2=sparse(ones(length(b_vect),1),ones(length(b_vect),1),f_vect);
    M2=diag(f);

    tic
    for j = 1 : length(y)
        for i = 1 : length(x)

            ind=indexeur(j,i);

            if i > 1
                M(ind,indexeur(j,i-1)) = c;
                M2(ind,indexeur(j,i-1)) = g;
            end

            if i < length(x)
                M(ind,indexeur(j,i+1)) = c;
                M2(ind,indexeur(j,i+1)) = g;
            end

            if j < length(y)
                M(ind,indexeur(j+1,i)) = d;
                M2(ind,indexeur(j+1,i)) = k;
            end

            if j > 1
                M(ind,indexeur(j-1,i)) = d;
                M2(ind,indexeur(j-1,i)) = k;
            end
        end
    end
    toc

    tic
    M=sparse(M);
    M2=sparse(M2);
    toc

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Init theorie
    
    Psy_x=zeros(1,length(x));
    Psy_y=zeros(1,length(y));
    Psy_th=zeros(length(t),length(y),length(x));
    
    [Psy_th(1,:,:),norm_th_y(w,1)]=analy(x,y,t(1));
    
    norm_th_y(w,1)=trapeze_2D(abs(squeeze(Psy_th(1,:,:))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
    err_y(w,1)=trapeze_2D(abs(squeeze(abs(Psy_th(1,:,:)))-squeeze(abs(Psy_mat(:,:,1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calcul

%     figure()
    tic
    for j = 1 : length(t)-1
        
        b=M2*Psy(:,j);
        Psy(:,j+1) = mldivide(M,b);
        
        psy=vec2mat(Psy(:,j+1),length(y))';
        Psy_mat(:,:,j+1)=psy;
        
        [Psy_th(j+1,:,:),norm_th_y(w,j+1)]=analy(x,y,t(j+1));
        
        norm_y(w,j+1)=trapeze_2D(abs(squeeze(Psy_mat(:,:,j+1))).^2,x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        err_y(w,j+1)=trapeze_2D(abs(squeeze(abs(Psy_th(j+1,:,:)))-squeeze(abs(Psy_mat(:,:,j+1)))),x(1),x(end),y(1),y(end),length(x)-1,length(y)-1);
        
%         subplot(1,2,1)
%         surf(x,y,abs(squeeze(Psy_mat(j,:,:))).^2,'edgecolor','none');
%         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_x(w,j)))
%         view(0,90);
%         subplot(1,2,2)
%         surf(x,y,abs(squeeze(Psy_th(j,:,:))).^2,'edgecolor','none');
%         title(sprintf('Temps = %e  Norme: %.10f',t(j),norm_th_x(w,j)))
%         view(0,90);
%         hold off
%         
%         pause(0.01)
        
    end
    toc  
    w    
end

figure()
hold on
for i=1:length(Dy)
    plot(t,norm_x(i,:))
end
hold off
legend(sprintf('dy = %d',Dy(1)),sprintf('dy = %d',Dy(2)),sprintf('dy = %d',Dy(3)))
title('Evolution de la norme de la resolution en fonction de DY')

figure()
hold on
for i=1:length(Dy)
    plot(t,norm_th_y(i,:))
end
hold off
legend(sprintf('dy = %d',Dy(1)),sprintf('dy = %d',Dy(2)),sprintf('dy = %d',Dy(3)))
title('Evolution de la norme de la solution analyique en fonction de DY')

figure()
hold on
for i=1:length(Dy)
    loglog(t,err_y(i,:))
end
hold off
legend(sprintf('dy = %d',Dy(1)),sprintf('dy = %d',Dy(2)),sprintf('dy = %d',Dy(3)))
title('Evolution de l''erreur en fonction de DY')

P_y=polyfit(log(Dy),log(err_y(:,end)'),1)
w
figure()
plot(log(Dy),log(err_y(:,end)))
hold on
plot(log(Dy),P_y(1)*log(Dy)+P_y(2))

