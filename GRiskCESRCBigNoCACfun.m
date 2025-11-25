function [ktildez0,cz0,vz0,kmax,kbar,ksd,kssM,kssno,cbar,csd,ibar,isd,ybar,ysd,v1,v2,v3,v4,iter,trashit,trashmax,v,ktilde,c,z,kmin,csdcompdata,caccompdata,cac,sigma,sdz] = GRiskCESRCBigfunNoCaC(sdz,sigma,A,alpha,rho,beta,mu,delta,gp,gpz,T,Ti,itermax,k,v,UsekssTh,kssTh,Reps)
% Same as standard, but stores more stuff.
% This solves the Markow growth problem with CES-utility by iterating on a
% discrete grid.
% Resource constraint
% Uncertain productivity z, n=gpz different states
% Two options for doing N approx, rouwen and markovappr

% %Parameters, CHECK THAT USUAL RATIOS WORK OUT

iter=0;
trashit=0;
trashmax=0;

% This part is for regular
%[P, z]=rouwen(rho, 0, sdz, gpz); % P(n,m) = Prob(z(t+1)=z_n|z(t)=z_m)
if gpz>31
    [P,z,~]=markovapprcj(rho,sdz*sqrt(1-rho^2),4,gpz);
else
    [P,z,~]=markovapprcj(rho,sdz*sqrt(1-rho^2),3,gpz);
end
z=z'; % Required with markovapprcj
P=P'; % Required with markovapprcj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is to test out CE, just two sdz values, >0 one must be changed
% manually in second line.
% if sdz>0
%     [P, z]=rouwen(rho, 0, sdz, gpz); % P(n,m) = Prob(z(t+1)=z_n|z(t)=z_m)
% else
%     [P, z]=rouwen(rho, 0, .3, gpz); % P(n,m) = Prob(z(t+1)=z_n|z(t)=z_m)
%     P=eye(gpz);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% v=zeros(gp,gpz); % v(i,m) value function before iteration, initially all zeros
nextv=ones(gp,gpz); %value function after iteration, initially all ones (must differ from initial v)
U0=nan(gp,gp,gpz); %U(i,j,m) stores present period utility of going from k(i),z(m) to k(j)

disp('Setting up U0...');
%Compute u0 first because it will not change as we iterate, so only do once
for i=1:gp
    for j=1:gp
        for m=1:gpz
%            if exp(z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-k(j) < 0 % or <= Non-pos consumption
            if (1+z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-k(j) < 0 % or <= Non-pos consumption
                U0(i,j,m)=-inf;
            else
                % U0(i,j,m)=((exp(z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-k(j))^(1-sigma)-1)/(1-sigma);
                % U0(i,j,m)=((exp(z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-k(j))^(1-sigma))/(1-sigma); % So U0 doesn't change sign around c=1
                % U0(i,j,m)=(((1+z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-k(j))^(1-sigma)-1)/(1-sigma);
                U0(i,j,m)=(((1+z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-k(j))^(1-sigma))/(1-sigma); % So U0 doesn't change sign around c=1
            end
        end
    end
end
disp('Done setting up U0!');
disp('Starting iterations...');


% keep iterating until the two values are close enough
% while max(abs(nextv-v),[],'all')/min(abs(nextv),[],'all')>tol
% while max(max(abs(nextv-v)))/min(min(abs(nextv)))>tol
while max(max(abs(nextv-v)))>0
    
    if iter~=0 %update v, except in first iteration
        v=nextv;
    end
    
    iter=iter+1; %counter
    
    if iter>itermax
        break; %stop if max iterations reached
    end

    V=nan(gp,gp,gpz);

    for m=1:gpz 
        M=zeros(gp,gpz);
        M(:,m)=ones(gp,1);
        V(:,:,m)=U0(:,:,m)+beta*M*P'*v';
    end
    clear M;

    %gp=100, gpz=11: Parfor2 .053, Parfor4 .089, for .003
    %gp=1000, gpz=21: Parfor2 1.58, Parfor4 1.44, for .55
    %gp=2000, gpz=21: Parfor2 6.7, Parfor4 5.8, for 2.2
    %gp=3000, gpz=21: Parfor2 17.2, Parfor4 15.2, for 4.5
    
%     for i=1:gp;
%         for j=1:gp;
%             for m=1:gpz;
%                 %V(i,j,m)= U0(i,j,m)+beta*(P(:,m)'*(v(j,:))');
%                 V(i,j,m)= U0(i,j,m)+beta*(P(1,m)*v(j,1)+P(2,m)*v(j,2)+P(3,m)*v(j,3));
%             end;
%         end;
%     end;
    
    V(isnan(V))=-inf;
    V=permute(V,[1 3 2]);
    [nextv, ktildeposition]=max(V,[],3);
    clear V; 
end

clear U0;
iter

while iter >= itermax
    disp('XXXXXXXXXXXX DID NOT CONVERGE XXXXXXXXXXXX');
    trashit=1;
    break;
end

while max(max(ktildeposition)) >=  gp
    disp('XXXXXXXXXXXX Kmax too low, increase kmaxgp XXXXXXXXXXXX');
    trashmax=1;
    break;
end


% Things below this line are for plots and displaying results

ktilde=k(ktildeposition);
c=nan(gp,gpz);

for i=1:gp
    for m=1:gpz
%        c(i,m)=exp(z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-ktilde(i,m);
        c(i,m)=(1+z(m))*A*(k(i))^alpha+(1-delta)*k(i)+mu-ktilde(i,m);
    end
end


% Finding the steady-state from computed rules, and feeding it into the
% simulation, or using theoretical ss
kssMpoints=find((1:gp)'==ktildeposition(:,round(gpz/2)));
kssM=sum(k(kssMpoints))/(size(kssMpoints,1)-1);
kssno=size(kssMpoints,1)-1;

if UsekssTh==1
    Dist1=abs(k-kssTh);
    idk1=round(find(Dist1==min(Dist1))); % theoretical kss
else
    Dist1=abs(k-kssM);
    idk1=find(Dist1==min(Dist1)); % computed kss, median of all ss points exluding 0
end


kmax_s=nan(Reps,1);
kmin_s=nan(Reps,1);
kbar_s=nan(Reps,1);
ksd_s=nan(Reps,1);
cbar_s=nan(Reps,1);
csd_s=nan(Reps,1);
ibar_s=nan(Reps,1);
isd_s=nan(Reps,1);
ybar_s=nan(Reps,1);
ysd_s=nan(Reps,1);
cac_s=nan(Reps,1);
csdcompdata_s=nan(Reps,1);
caccompdata_s=nan(Reps,1);
v3_s=nan(Reps,1);

disp('Starting simulations...');

for r = 1:Reps
    
    sim=rand(Ti+T,1);
    Zsim=hitm_zn(z,P,sim); %simulated values of z, initial always z=E(z)
    Zpositionsim=hitm_sn(P,sim); %simulated positions in z-vector, consistent with sims from previous line
    clear sim;
    
    Kpositionsim=nan(Ti+T+1,1); %will hold simulated values for k
    Kpositionsim(1,1)=idk1(1,1); %initial capital stock position
    
    % Applying optimal decision rule
    for t=1:Ti+T
        Kpositionsim(t+1,1)=ktildeposition(Kpositionsim(t,1),Zpositionsim(t));
    end
    
    Kpositionsim=Kpositionsim(Ti+1:Ti+T);
    Ksim=k(Kpositionsim(1:T));
    clear Kpositionsim Zpositionsim;
    Zsim=Zsim(Ti+1:Ti+T);
    
    
    % for t=1:T-1
    %     Csim(t,1)=exp(Zsim(t,1))*A*(Ksim(t,1))^alpha+(1-delta)*Ksim(t,1)+mu-Ksim(t+1,1);
    %     Isim(t,1)=Ksim(t+1,1)-(1-delta)*Ksim(t,1);
    %     Ysim(t,1)=exp(Zsim(t,1))*A*(Ksim(t,1))^alpha;
    % end
    
    %Ysim=A*exp(Zsim).*(Ksim.^alpha);
    Ysim=A*(1+Zsim).*(Ksim.^alpha);
    Csim=Ysim+(1-delta)*Ksim+mu-lagmatrix(Ksim,-1);
    Csim=Csim(1:T-1);
    Isim=Ysim(1:T-1)-Csim;
    
%     kmax=max(Ksim);
%     kmin=min(Ksim);
%     kbar=mean(Ksim);
%     ksd=std(Ksim);
%     cbar=mean(Csim);
%     csd=std(Csim);
%     ibar=mean(Isim);
%     isd=std(Isim);
%     ybar=mean(Ysim);
%     ysd=std(Ysim);
    
    kmax_s(r)=max(Ksim);
    kmin_s(r)=min(Ksim);
    kbar_s(r)=mean(Ksim);
    ksd_s(r)=std(Ksim);
    cbar_s(r)=mean(Csim);
    csd_s(r)=std(Csim);
    ibar_s(r)=mean(Isim);
    isd_s(r)=std(Isim);
    ybar_s(r)=mean(Ysim);
    ysd_s(r)=std(Ysim);
    
    clear Ksim Isim Ysim;
    
    % acCsim = autocorr(Csim,'NumLags',1);
    % %cac = acCsim(2);
    % cac_s(r) = acCsim(2);
    cac_s(r) = nan;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % HP-filtering and % deviaiton from trend before computing stats, as in
    % % data
    csdcompdata = nan;
    caccompdata = nan;
    
    % % [Ksimtrend, Ksimcycle]=hpfilter(Ksim,1600);
    % % Ksimhat=Ksimcycle./Ksimtrend*100;
    % [Csimtrend, Csimcycle]=hpfilter(Csim,1600);
    % Csimhat=Csimcycle./Csimtrend*100;
    % % [Isimtrend, Isimcycle]=hpfilter(Isim,1600);
    % % Isimhat=Isimcycle./Isimtrend*100;
    % % [Ysimtrend, Ysimcycle]=hpfilter(Ysim,1600);
    % % Ysimhat=Ysimcycle./Ysimtrend*100;
    % clear Csimtrend Csimcycle Csim
    
    % % ksd=std(Ksimhat);
    % % csd=std(Csimhat);
    % % isd=std(Isimhat);
    % % ysd=std(Ysimhat);
    
    % % csdcompdata = std(Csimhat);
    % csdcompdata_s(r) = std(Csimhat);
    % acCsimhat = autocorr(Csimhat,'NumLags',1);
    % % caccompdata = acCsimhat(2);
    % caccompdata_s(r) = acCsimhat(2);

    clear Csimhat;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Simulated value function to test computed by iterations
    % Sim=1e5;
    % Tsim=2e4+1;
    % Kpositionsim=nan(Tsim,1); %will hold simulated values for k
    % Kpositionsim(1,1)=idk1(1,1); %initial capital stock position
    % Beta=beta.^(0:Tsim-2);
    % Vsim=nan(Sim,1);
    % Kbar=nan(Sim,1);
    % Ksd=nan(Sim,1);
    % Cbar=nan(Sim,1);
    % Csd=nan(Sim,1);
    % for s=1:Sim
    %     sim=rand(Tsim,1);
    %     Zsim=hitm_zn(z,P,sim); %simulated values of z, initial always z=E(z)
    %     Zpositionsim=hitm_sn(P,sim); %simulated positions in z-vector, consistent with sims from previous line
    %     for t=1:Tsim-1
    %         Kpositionsim(t+1,1)=ktildeposition(Kpositionsim(t,1),Zpositionsim(t));
    %     end
    % Ksim=k(Kpositionsim(1:Tsim)); % Applying optimal decision rule
    % Csim=A*exp(Zsim).*(Ksim.^alpha)+(1-delta)*Ksim+mu-lagmatrix(Ksim,-1);
    % Usim=((Csim(1:Tsim-1)).^(1-sigma)-1)/(1-sigma);
    % Vsim(s,1)=Beta*Usim;
    % Kbar(s,1)=mean(Ksim);
    % Ksd(s,1)=std(Ksim);
    % Cbar(s,1)=mean(Csim(1:Tsim-1));
    % Csd(s,1)=std(Csim(1:Tsim-1));
    % end
    % kbarsim=mean(Kbar)
    % ksdsim=mean(Ksd)
    % cbarsim=mean(Cbar)
    % csdsim=mean(Csd)
    % meanv=mean(Vsim)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Dist2=abs(k-kbar);
    Dist2=abs(k-kbar_s(r));
    idk2=find(Dist2==min(Dist2)); % kbar mean simulated k
    % v3=v(idk2(1,1),round(gpz/2));  % v(kbar,E(z))
    v3_s(r)=v(idk2(1,1),round(gpz/2));  % v(kbar,E(z))
    % v4=v(idk2(1,1),:)*diag(P10000); % E(v(kbar,z))
    % v4_s(r)=v(idk2(1,1),:)*diag(P10000); % E(v(kbar,z))
    
end 

disp('Done with simulations!');
clear ktildeposition;

kmax=max(kmax_s);
kmin=min(kmin_s);
kbar=mean(kbar_s);
ksd=mean(ksd_s);
cbar=mean(cbar_s);
csd=mean(csd_s);
ibar=mean(ibar_s);
isd=mean(isd_s);
ybar=mean(ybar_s);
ysd=mean(ysd_s);
cac=mean(cac_s);
csdcompdata=mean(csdcompdata_s);
caccompdata=mean(caccompdata_s);
v3=mean(v3_s);


v1=v(idk1(1,1),round(gpz/2)); % v(kss,E(z))
v2=nan; % E(v(kss,z))
v4=nan; % E(v(kbar,z))
%P10000=P^10000;
%v2=v(idk1(1,1),:)*diag(P10000); % E(v(kss,z))

ktildez0=ktilde(:,round(gpz/2));
cz0=c(:,round(gpz/2));
vz0=v(:,round(gpz/2));