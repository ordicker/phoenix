% 1) According to the notation and coordinate system of "An_Introduction_to_Graphene_Plasmonics" page 42
% 2) Another good source is "Transfer matrix method for optics in graphene layers", Journal of Physics: Condensed Matter, Volume 25, Number 21

clear;
close all;

%% Constants
constants;

%% Params
n=1; 
%lambda = linspace(3e-6,5.5e-6,n); %in [m]
lambda = 4e-6;
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
eps3xx=1.5209^2*ones(1,n);
eps3zz=1.5209^2*ones(1,n);
eps2xx=(4.90+11.70i)^2*ones(1,n);
eps2zz=(4.90+11.70i)^2*ones(1,n);
eps1xx=3.4699^2*ones(1,n);
eps1zz=3.4699^2*ones(1,n);


q=k0.*eps1xx.^0.5.*sind(theta); % for regular TMM
[Q,O] = meshgrid(q,omega); % make matrices for the TM

% Define number of layers
N=3;

% Transfrom eps to matrices
for i=1:N
    k = [ ' EPS' num2str(i) 'XX=repmat(eps' num2str(i) 'xx,n,1).'';' ];
    eval(k);
    k = [ ' EPS' num2str(i) 'ZZ=repmat(eps' num2str(i) 'zz,n,1).'';' ];
    eval(k);
end

NN = 500
AAA=zeros(NN,1) 
RRR=zeros(NN,1)
TTT=zeros(NN,1)
DDD=linspace(1e-4,1e-1,NN)
counter = 1
for dd=linspace(1e-4,1e-1,NN)
% Define thicknesses, size N-2
d = [  dd  ]*1e-6; % in [m]
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


TTT(counter)=T;
RRR(counter)=R;
AAA(counter)=A;
counter=counter+1;
end


%% plot
set(0,'DefaultFigureWindowStyle','normal');
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultLineMarkerSize',8)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureColor','white');
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
 
%figure, plot(lambda, R(:,1),lambda, T(:,1),lambda, A(:,1));  legend('R','T','A')

figure, plot(DDD, RRR, DDD, TTT, DDD, AAA);  legend('R','T','A')
