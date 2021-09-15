clc;
clear all;
%Matlab_code_3node_triangle_heat_transfer/Thermal.m
k = 30;  h = 60; Tinf = 25; Q2 = 10^6; Q1=0; qn=2*10^5;
% Gauss Points and weights for one point gauss quadrature
% for triangular area integration
% -------------------------------------------------------
xi1 = 1/3;  eta1 = 1/3;  w1 = 1;

% Gausspoints for line integration
% --------------------------------
beta1 = -0.774597;   ew1 = 5/9;
beta2 = 0.774597;    ew2 = 5/9;
beta3 = 0;           ew3 = 8/9;

% Element 1
coord = [0.0, 0.0;
         0.5, 0.0;
         0.0, 0.6];
[ke1d,fe1d] = domain(xi1, eta1, coord, k, Q1);

[ke1hg1,fe1hg1] = gamah(beta1, coord, h, Tinf, 1);
[ke1hg2,fe1hg2] = gamah(beta2, coord, h, Tinf, 1);
[ke1hg3,fe1hg3] = gamah(beta3, coord, h, Tinf, 1);
ke1h = ke1hg1*ew1 + ke1hg2*ew2 + ke1hg3*ew3;
fe1h = fe1hg1*ew1 + fe1hg2*ew2 + fe1hg3*ew3;

ke1 = ke1d + ke1h;
fe1 = fe1d + fe1h;


% Element 2
coord = [0.5, 0.0;
         0.8, 0.0;
         0.8, 0.3];
[ke2d,fe2d] = domain(xi1, eta1, coord, k, Q2);

[ke2hg1,fe2hg1] = gamah(beta1, coord, h, Tinf, 1);
[ke2hg2,fe2hg2] = gamah(beta2, coord, h, Tinf, 1);
[ke2hg3,fe2hg3] = gamah(beta3, coord, h, Tinf, 1);
ke2h = ke2hg1*ew1 + ke2hg2*ew2 + ke2hg3*ew3;
fe2h = fe2hg1*ew1 + fe2hg2*ew2 + fe2hg3*ew3;

[fe2q1] = gamaq(beta1, coord, qn, 2);
[fe2q2] = gamaq(beta2, coord, qn, 2);
[fe2q3] = gamaq(beta3, coord, qn, 2);
fe2q = fe2q1*ew1 + fe2q2*ew2 + fe2q3*ew3;

ke2 = ke2d + ke2h;
fe2 = fe2d + fe2h +fe2q;

% Element 3
coord = [0.8, 0.3;
         0.5, 0.3;
         0.5, 0.0];
[ke3d,fe3d] = domain(xi1, eta1, coord, k, Q2);

[fe3q1] = gamaq(beta1, coord, qn, 1);
[fe3q2] = gamaq(beta2, coord, qn, 1);
[fe3q3] = gamaq(beta3, coord, qn, 1);
fe3q = fe3q1*ew1 + fe3q2*ew2 + fe3q3*ew3;

ke3 = ke3d ;
fe3 = fe3d + fe3q;


% Element 4
coord = [0.5, 0.0;
         0.5, 0.6;
         0.0, 0.6];
[ke4d,fe4d] = domain(xi1, eta1, coord, k, Q1);

[ke4hg1,fe4hg1] = gamah(beta1, coord, h, Tinf, 2);
[ke4hg2,fe4hg2] = gamah(beta2, coord, h, Tinf, 2);
[ke4hg3,fe4hg3] = gamah(beta3, coord, h, Tinf, 2);
ke4h = ke4hg1*ew1 + ke4hg2*ew2 + ke4hg3*ew3;
fe4h = fe4hg1*ew1 + fe4hg2*ew2 + fe4hg3*ew3;

ke4 = ke4d + ke4h;
fe4 = fe4d + fe4h;


% Assembly
K = zeros(7,7);
F = zeros(7,1);

K([1,2,7],[1,2,7]) = ke1(1:3,1:3);
K([2,3,4],[2,3,4]) = K([2,3,4],[2,3,4]) + ke2(1:3,1:3);
K([4,5,2],[4,5,2]) = K([4,5,2],[4,5,2]) + ke3(1:3,1:3);
K([2,6,7],[2,6,7]) = K([2,6,7],[2,6,7]) + ke4(1:3,1:3);

F([1,2,7]) = fe1(1:3);
F([2,3,4]) = F([2,3,4]) + fe2(1:3);
F([4,5,2]) = F([4,5,2]) + fe3(1:3);
F([2,6,7]) = F([2,6,7]) + fe4(1:3);


% Imposition of B.C.
Kreduce = K(2:6,2:6)
Freduce = F(2:6) - (K(1,2:6)*200+ K(7,2:6)*200)'

% Finding Solution
ureduce = inv(Kreduce)*Freduce;
un = [200; ureduce ; 200];
un


% % ==============================================
% %% Printing Intermediate Result to The Output File
% % ------------------------------------------------
fid=fopen('Steps_Q1_(a).txt','w');
fprintf(fid,'The Element Stiffness matrices are\n');
fprintf(fid,'===================================\n');
fprintf(fid,'k = %12.4e, h = %12.4e, Tinf = %12.4e, Q = %12.4e\n\n',k,h,Tinf,Q2);
fprintf(fid,'Ke1 \n');
fprintf(fid,'----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',ke1(i,1:3));
end
fprintf(fid,'Ke2 \n');
fprintf(fid,'----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',ke2(i,1:3));
end
fprintf(fid,'Ke3 \n');
fprintf(fid,'----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',ke3(i,1:3));
end
fprintf(fid,'Ke4 \n');
fprintf(fid,'----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\n\n',ke4(i,1:3));
end


fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');

fprintf(fid,'fe1 \n');
fprintf(fid,'-----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',fe1(i));
end
fprintf(fid,'fe2 \n');
fprintf(fid,'-----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',fe2(i));
end
fprintf(fid,'fe3 \n');
fprintf(fid,'-----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',fe3(i));
end
fprintf(fid,'fe4 \n');
fprintf(fid,'-----\n\n');
for i = 1:3
   fprintf(fid,'%14.4e\n\n',fe4(i));
end


fprintf(fid,'\n\nThe Global Stiffness matrix is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K\n');
fprintf(fid,'--\n');
for i = 1:7
   fprintf(fid,'%12.4e\t8%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\n\n',K(i,1:7));
end

fprintf(fid,'\n\nThe Global Load Vector is\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'F\n');
fprintf(fid,'--\n');
for i = 1:7
   fprintf(fid,'%12.4e\n\n',F(i));
end

fprintf(fid,'\n\nImposition of Boundary Condition\n');
fprintf(fid,'==================================\n\n');
fprintf(fid,'K u = F\n');
fprintf(fid,'--------\n');
for i = 1:7
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',K(i,1:7),F(i));
end

fprintf(fid,'\n\nReduced Equations\n');
fprintf(fid,'=======================\n\n');
fprintf(fid,'K_reduce u = F_reduce\n');
fprintf(fid,'--------\n');
for i = 1:5
   fprintf(fid,'%12.4e\t%12.4e\t%12.4e\t%12.4e\t%12.4e\t\t\t\t%12.4e\n\n',Kreduce(i,1:5),Freduce(i));
end

fprintf(fid,'\n\nThe Final Solution\n');
fprintf(fid,'=========================\n\n');
fprintf(fid,'Tn\n');
fprintf(fid,'--\n');
for i = 1:7
   fprintf(fid,'%12.4e\n\n',un(i));
end
function [ke,fe] = domain(xi, eta, coord, k, Q)

% This function calculate element Stiffness Matrix

% INPUT:
% ======
% k = Thermal_conductivity
% Q = Heat generation per unit volume
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% ke = Element stiffness matrix
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);


J = [x(1)-x(3), y(1)-y(3); x(2)-x(3), y(2)-y(3)];
N = [xi, eta, 1-xi-eta];
B = inv(J)*[1, 0, -1; 0, 1, -1];


ke = 0.5*k*det(J)*B'*B;

fe = N'*0.5*Q*det(J);
end
function [ke,fe] = gamah(beta, coord, h, Tinf, edgeno)

% This function calculate loadvector on gamaq part

% INPUT:
% ======
% h = Convective heat transfer coefficient
% Tinf = Ambient temperature
% beta = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge heat flux is prescribed
%
% OUTPUT:
% =======
% ke = Element stiffness matrix
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);
xi = (1+beta)/2;

if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N = [1-xi, xi, 0];
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N = [0, 1-xi, xi];
elseif edgeno == 3
    ledge = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
    N = [xi, 0, 1-xi];
end

fe = N'*h*ledge*Tinf*0.5;
ke = h*N'*N*ledge*0.5;
  
end
function [fe] = gamaq(beta, coord, qn, edgeno)

% This function calculate loadvector on gamaq part

% INPUT:
% ======
% qn = Heat flux 
% beta = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge heat flux is prescribed
%
% OUTPUT:
% =======
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);
xi = (1+beta)/2;

if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N = [1-xi, xi, 0];
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N = [0, 1-xi, xi];
elseif edgeno == 3
    ledge = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
    N = [xi, 0, 1-xi];
end

fe = N'*qn*ledge*0.5;
end


