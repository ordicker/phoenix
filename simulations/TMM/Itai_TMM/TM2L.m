% Calulates the transfer matrix between two material with a graphene layer in the middle

function [TM, k1z, k2z] = TM2L(q, omega, s, eps1xx, eps1zz, eps2xx, eps2zz)
% eps_in is the first layer by definition 

c=299792458; %m/s
eps0=8.854e-12;

% k1z=sqrt(omega.^2/c^2.*eps1xx - q.^2.*eps1xx./eps1zz).*sqrt(1-((eps_in./eps1xx).^0.5.*sind(theta)).^2);
% k2z=sqrt(omega.^2/c^2.*eps2xx - q.^2.*eps2xx./eps2zz).*sqrt(1-((eps_in./eps2xx).^0.5.*sind(theta)).^2);
k1z=sqrt(omega.^2/c^2.*eps1xx - q.^2.*eps1xx./eps1zz);
k2z=sqrt(omega.^2/c^2.*eps2xx - q.^2.*eps2xx./eps2zz);
eta=eps1xx.*k2z./(eps2xx.*k1z);
xsi=s.*k2z./(omega.*eps0.*eps2xx);

T11=0.5*(1+eta+xsi);
T12=0.5*(1-eta-xsi);
T21=0.5*(1-eta+xsi);
T22=0.5*(1+eta-xsi);

TM = [T11 T12; T21 T22];

end 
