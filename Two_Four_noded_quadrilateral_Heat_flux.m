clc;
clear all;

k = 30;  h = 60; Tinf = 25; Q1 = 0; Q2 = 10^6 ; qn= 2*10^5;


% Gausspoints for 3x3 line integration
% --------------------------------
xi = [-0.774597, 0, 0.774597];
xiw = [5/9, 8/9, 5/9];

eta = [-0.774597, 0, 0.774597];
etaw = [5/9, 8/9, 5/9];

% Element 1
coord = [0.0, 0.0;
         0.5, 0.0;
         0.5, 0.6;
         0.0, 0.6];
     
 % Sum over gauss points:
 ke1 = zeros(4,4);
 fe1 =zeros(4,1);
 for i = 1:3
     for j = 1:3
         [ke1d , fe1d] = domain(xi(i), eta(j), coord, k, Q1);
         ke1 = ke1 + ke1d*xiw(i)*etaw(j);
         fe1 = fe1 + fe1d*xiw(i)*etaw(j);
     end
 end
 
  for i = 1:3
     for j = 1:3
         [ke1h , fe1h] = gamah(xi(i), coord, h, Tinf, 1);
         ke1 = ke1 + ke1h*xiw(i)*etaw(j);
         
     end
     fe1 = fe1 + fe1h*xiw(i);
 end
 
   for i = 1:3
     for j = 1:3
         [ke1h , fe1h] = gamah(xi(i), coord, h, Tinf, 3);
         ke1 = ke1 + ke1h*xiw(i)*etaw(j);
         
     end
     fe1 = fe1 + fe1h*xiw(i);
 end
 

 
 % Element 2
coord = [0.5, 0.0;
         0.8, 0.0;
         0.8, 0.3;
         0.5, 0.3];
     
 % Sum over gauss points:
 ke2 = zeros(4,4);
 fe2 =zeros(4,1);
 for i = 1:3
     for j = 1:3
         [ke2d , fe2d] = domain(xi(i), eta(j), coord, k, Q2);
         ke2 = ke2 + ke2d*xiw(i)*etaw(j);
         fe2 = fe2 + fe2d*xiw(i)*etaw(j);
     end
 end
 
   for i = 1:3
     for j = 1:3
         [ke2h , fe2h] = gamah(xi(i), coord, h, Tinf, 1);
         ke2 = ke2 + ke2h*xiw(i)*etaw(j);
        
     end
      fe2 = fe2 + fe2h*xiw(i);
 end
    
 
     for i = 1:3
         [ fe2q] = gamaq(xi(i), coord, qn, 2);
           fe2 = fe2 + fe2q*xiw(i);
     end
     


  
     for i = 1:3
         [ fe2q] = gamaq(xi(i), coord, qn, 3);
        
         fe2 = fe2 + fe2q*xiw(i);
     end




% Assembly
K = zeros(7,7);
F = zeros(7,1);

K([1,2,6,7],[1,2,6,7]) = ke1(1:4,1:4);
K([2,3,4,5],[2,3,4,5]) = K([2,3,4,5],[2,3,4,5]) + ke2(1:4,1:4);


F([1,2,6,7]) = fe1(1:4);
F([2,3,4,5]) = F([2,3,4,5]) + fe2(1:4);



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
fid=fopen('Steps_Q1_(c).txt','w');
fprintf(fid,'The Element Stiffness matrices are\n');
fprintf(fid,'===================================\n');
fprintf(fid,'k = %12.4e, h = %12.4e, Tinf = %12.4e, Q = %12.4e\n\n',k,h,Tinf,Q2);
fprintf(fid,'Ke1 \n');
fprintf(fid,'----\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',ke1(i,1:4));
end
fprintf(fid,'Ke2 \n');
fprintf(fid,'----\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\t%14.4e\t%14.4e\t%14.4e\n\n',ke2(i,1:4));
end



fprintf(fid,'The Element Load Vector are\n');
fprintf(fid,'===================================\n\n');

fprintf(fid,'fe1 \n');
fprintf(fid,'-----\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\n\n',fe1(i));
end
fprintf(fid,'fe2 \n');
fprintf(fid,'-----\n\n');
for i = 1:4
   fprintf(fid,'%14.4e\n\n',fe2(i));
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


function [ke, fe] = domain(xi, eta, coord, k, Q)

% This function calculate element Stiffness Matrix

% INPUT:
% ======
% C = Constitutive Matrix
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% ke = Element stiffness matrix

x = coord(:,1);
y = coord(:,2);

ddmat = 0.25* [-(1-eta), (1-eta), (1+eta), -(1+eta);
               -(1-xi), -(1+xi),  (1+xi),   (1-xi)];

N1=0.25*(1-xi)*(1-eta);
N2=0.25*(1+xi)*(1-eta);
N3=0.25*(1+xi)*(1+eta);
N4=0.25*(1-xi)*(1+eta);
N=[N1, N2, N3, N4];

J = ddmat*coord;

B = inv(J)*ddmat;

ke = k*det(J)*(B'*B);
fe= N'*Q*det(J);


end


function [ke, fe] = gamah(xi, coord, h, Tinf, edgeno)

% This function calculate loadvector on gamaq part

% INPUT:
% ======
% tvec = traction components
% xi = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge traction is prescribed
%
% OUTPUT:
% =======
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);
xi=(1+xi)/2;
if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N1 = (1-xi)/2; N2 = (1+xi)/2; N3 = 0; N4 = 0;
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N1 = 0; N2 = (1-xi)/2; N3 = (1+xi)/2; N4 = 0;
elseif edgeno == 3
    ledge = sqrt((x(4)-x(3))^2 + (y(4)-y(3))^2);
    N1 =0; N2=0; N3 = (1+xi)/2; N4 = (1-xi)/2;
elseif edgeno == 4
    ledge = sqrt((x(4)-x(1))^2 + (y(4)-y(1))^2);
    N1 = (1-xi)/2; N2 = 0; N3 = 0; N4 = (1+xi)/2;
end

N = [N1, N2, N3, N4];


fe = 0.5*ledge*N'*h*Tinf;
ke= 0.5*h*N'*N*ledge;
  
end

function [ fe] = gamaq(xi, coord, qn, edgeno)

% This function calculate loadvector on gamaq part

% INPUT:
% ======
% tvec = traction components
% xi = Gausspoint for the limit -1 to 1.
% coord = Nodal coordinates of the element
% edgeno = On which edge traction is prescribed
%
% OUTPUT:
% =======
% fe = Element load vector

x = coord(:,1);
y = coord(:,2);
xi=(1+xi)/2;
if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N1 = (1-xi)/2; N2 = (1+xi)/2; N3 = 0; N4 = 0;
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N1 = 0; N2 = (1-xi)/2; N3 = (1+xi)/2; N4 = 0;
elseif edgeno == 3
    ledge = sqrt((x(4)-x(3))^2 + (y(4)-y(3))^2);
    N1 =0; N2=0; N3 = (1+xi)/2; N4 = (1-xi)/2;
elseif edgeno == 4
    ledge = sqrt((x(4)-x(1))^2 + (y(4)-y(1))^2);
    N1 = (1-xi)/2; N2 = 0; N3 = 0; N4 = (1+xi)/2;
end

N = [N1, N2, N3, N4];


fe = 0.5*ledge*N'*qn;

  
end
