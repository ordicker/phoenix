% performs block multiplication of two square matrices  of size nXn
function [TT] = Block_multi(A, B, n)

A1=A(1:n,1:n);
A2=A(1:n,n+1:2*n);
A3=A(n+1:2*n,1:n);
A4=A(n+1:2*n,n+1:2*n);

B1=B(1:n,1:n);
B2=B(1:n,n+1:2*n);
B3=B(n+1:2*n,1:n);
B4=B(n+1:2*n,n+1:2*n);

T11=A1.*B1+A2.*B3;
T12=A1.*B2+A2.*B4;
T21=A3.*B1+A4.*B3;
T22=A3.*B2+A4.*B4;

TT=[T11 T12; T21 T22];
end