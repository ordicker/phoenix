% 1) According to the notation and coordinate system of "An_Introduction_to_Graphene_Plasmonics" page 42
% 2) Another good source is "Transfer matrix method for optics in graphene layers", Journal of Physics: Condensed Matter, Volume 25, Number 21

clear;
close all;

%% Constants
constants;

%% Params
n=500;
% lambda = 7.437e-6;
lambda = linspace(3e-6,5.5e-6,n); %in [m]
% lam_sp = linspace(10e-9,1000e-9,n);
% q=2*pi./lam_sp;
%f = linspace(20,70,n)*1e12; %in [Hz]
%omega = 1*linspace(0.*Ef/hbar,1.5*Ef/hbar,n); %in [m]
%f=omega/2/pi;
f=c./lambda; % frequencies in Hz
lambda = c./f; %in [m]
omega=2*pi*f;
k0=2*pi./lambda;
theta = 0*linspace(1,89,n); % angle in degrees 
S=0;
%% Layers properties
% Define epsilons of layers

% eps2xx=1*eps_Si(f);
% eps2zz=1*eps_Si(f);
eps3xx=1*ones(1,n);
eps3zz=1*ones(1,n);
eps2xx=(-10.645+17.191i)*ones(1,n);
eps2zz=(-10.645+17.191i)*ones(1,n);
eps1xx=1*ones(1,n);
eps1zz=1*ones(1,n);
% eps2xx=eps_BN_x;
% eps2zz=eps_BN_z;
% eps3xx=eps_BN_x;
% eps3zz=eps_BN_z;
% eps2xx=1.*ones(1,n);
% eps2zz=1.*ones(1,n);
% eps4xx=-10000;
% eps4zz=-10000;
% eps3xx=eps_Ge(f);
% eps3zz=eps_Ge(f);
% eps4xx=ones(1,n);
% eps4zz=ones(1,n);
% eps4xx=eps_BN_x;
% eps4zz=eps_BN_z;
% eps3xx=ones(1,n);% omega in 1/cm
% eps3zz=ones(1,n);% omega in 1/cm
% eps5xx=ones(1,n);
% eps5zz=ones(1,n);
% eps6xx=ones(1,n);
% eps6zz=ones(1,n);

q=k0.*eps1xx.^0.5.*sind(theta); % for regular TMM
[Q,O] = meshgrid(q,omega); % make matrices for the TM

% Define number of layers
N=3;

% Transfrom eps to matrixes
for i=1:N
    k = [ ' EPS' num2str(i) 'XX=repmat(eps' num2str(i) 'xx,n,1).'';' ];
    eval(k);
    k = [ ' EPS' num2str(i) 'ZZ=repmat(eps' num2str(i) 'zz,n,1).'';' ];
    eval(k);
end

% Define thicknesses, size N-2
d = [  20.0  ]*1e-9; % in [m]

% Define at which interface is there a graphene layer, size N-1
Ng = [  ];

%% Transsmision matrix method for TM single layer graphene
for m=N:-1:1
    
    % put sigma at interfaces with graphene
    if isempty(Ng(Ng == m-1)) == 1
        flag1 = 0;
    else
        flag1 = 1;
    end
    if m > 1
        k = [ '[TM' num2str(m-1) num2str(m) ',k' num2str(m-1) 'z, k' num2str(m) 'z]=TM2L(Q,O,S*flag1,EPS' num2str(m-1) 'XX,EPS' num2str(m-1) 'ZZ,EPS' num2str(m) 'XX,EPS' num2str(m) 'ZZ);' ];
        eval(k);
    end
    if m > 2
        k = [ ' PM' num2str(m-1) '=PM(d(' num2str(m-2) '),k' num2str(m-1) 'z);' ];
        eval(k);
    end
    
end

%% perform multiplicatino of all the matrices.
k = [ 'LAST=TM' num2str(N-1) num2str(N) ';' ];
eval(k);
if N == 2
    PM=PM(0,zeros(size(k1z)));
    FIRST = PM;
    d=0;
else
    k = [ 'FIRST=PM' num2str(N-1) ';' ];
    eval(k);
end
TT = Block_multi(FIRST, LAST, n);
LAST=TT;

for l=N-1:-1:2
    
    k = [ 'FIRST=TM' num2str(l-1) num2str(l) ';' ];
    eval(k);
    TT = Block_multi(FIRST, LAST, n);
    LAST=TT;
    if l > 2
        k = [ 'FIRST=PM' num2str(l-1) ';' ];
        eval(k);
        TT = Block_multi(FIRST, LAST, n);
        LAST=TT;
    end

end

% Sort element blocks of the final matrix
T11=TT(1:n,1:n);
T12=TT(1:n,n+1:2*n);
T21=TT(n+1:2*n,1:n);
T22=TT(n+1:2*n,n+1:2*n);

k = ['kOUTz=k' num2str(N) 'z;'];
eval(k);
k = ['EPS_OUT=repmat(eps' num2str(N) 'xx,n,1).'';'];
eval(k);
EPS_IN=repmat(eps1xx,n,1).';
t=(1./T11); % t
T=EPS_IN.*kOUTz./(EPS_OUT.*k1z).*abs(1./T11).^2; % T
r=(T21./T11); % reflection ceoficient
R=abs(T21./T11).^2; % R
A=1-T-R; % absrobtion


%% plot
set(0,'DefaultFigureWindowStyle','normal');
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineMarkerSize',8)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','white');
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')

figure, plot(lambda, R(:,1),lambda, T(:,1),lambda, A(:,1));  legend('R','T','A')

% figure, imagesc(q/qF,omega*hbar_eV/Ef_ev, log10(10+imag(r))/log10(100)); set(gca,'YDir','normal'); colormap pink(512) ;colorbar;xlabel('q/qF'); ylabel('hbar*omega/Ef')
% figure, imagesc(q/qF,omega*hbar_eV/Ef_ev, imag(r)); set(gca,'YDir','normal'); colormap pink(512);colorbar;xlabel('q/qF'); ylabel('hbar*omega/Ef')
% hold on; line([0 2],[0 2],'LineStyle','-','Color','r'); line([0 1],[2 1],'LineStyle','-','Color','r'); %caxis([0.5 0.501])
% figure, imagesc(q,f*1e-12, imag(r)); colormap pink(512);colorbar;xlabel('q[1/m]'); ylabel('f [THz]'); %set(gca,'YDir','normal');
%  hold on; line([0 6e8],[26.068 26.068],'LineStyle','-','Color','b');

%% analytical solution

% %% for singel layer graphene, from page 73, only for Drude conductivity!
% % alpha=e^2/(4*pi*eps0*hbar*c);
% epsAvg=(eps1xx+eps2xx)/2;
% % OMEGA=sqrt(2*alpha*Ef*hbar*c*q./epsAvg)/hbar;
% qq=1i*2*omega.*epsAvg.*eps0./sigma;
%
% hold on;
% plot(qq/qF, omega*hbar_eV/Ef_ev);
% % plot(q, 1e4./(2*pi*c./OMEGA)*1e-6); xlim([0 1e8]); ylim([600 1800]);%for Basov


% %% solve numerically for double layer graphene , 3 materials, page 89
% % does not work for BN, maybe this eq does not take into account phonons
%
% maxNumOfModes=15;
% for i=1:length(omega)
%     k1z=(-omega(i).^2/c^2.*eps1xx + q.^2.*eps1xx./eps1zz).^0.5;
%     k2z=(-omega(i).^2/c^2.*eps2xx + q.^2.*eps2xx./eps2zz).^0.5;
%     k3z=(-omega(i).^2/c^2.*eps3xx + q.^2.*eps3xx./eps3zz).^0.5;
%     DR = exp(k2z*d).*(eps3xx./k3z+eps2xx./k2z+1i*sigma(i)./omega(i)/eps0).*(eps1xx./k1z+eps2xx./k2z+1i*sigma(i)./omega(i)/eps0)...
%         -exp(-k2z*d).*(eps3xx./k3z-eps2xx./k2z+1i*sigma(i)./omega(i)/eps0).*(eps1xx./k1z-eps2xx./k2z+1i*sigma(i)./omega(i)/eps0);
%     criteria=find(diff(sign(real(DR)))~=0);
%     %     figure, plot(q,real(DR));
%     %         figure, plot(q,real(exp(k2z*d).*(eps3xx./k3z+eps2xx./k2z+1i*sigma1(i)./omega(i)/eps0).*(eps1xx./k1z+eps2xx./k2z+1i*sigma1(i)./omega(i)/eps0))); hold on;
%     %         plot(q,real(exp(-k2z*d).*(eps3xx./k3z-eps2xx./k2z+1i*sigma1(i)./omega(i)/eps0).*(eps1xx./k1z-eps2xx./k2z+1i*sigma1(i)./omega(i)/eps0)));
%
%     if length(criteria) < maxNumOfModes
%         criteria(length(criteria)+1:maxNumOfModes)=1;
%     end
%     for j=1:maxNumOfModes
%         k = [ ' q' num2str(j) '(:,i)=q(criteria(' num2str(j) '));' ];
%         eval(k);
%     end
% end
%
% hold on;
% plot(q1/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q2/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q3/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q4/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q5/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q6/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q7/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q8/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q9/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q10/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q11/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q12/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q13/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q14/qF, omega*hbar_eV/Ef_ev,'b');
% plot(q15/qF, omega*hbar_eV/Ef_ev,'b');
