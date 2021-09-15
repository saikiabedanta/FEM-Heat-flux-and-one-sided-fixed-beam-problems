% Parameters
E = 200e9; nu = 0.3; P0 = 2e3; P1 = 4e3; theta = 20;
nele = 2;         % no of elements
nnode = 9;        % no of nodes

% Gauss Points and weights for one point gauss quadrature
% for triangular area integration
% -------------------------------------------------------
ngaussd = 4;
xivec = [1/3, 0.6, 0.2, 0.2];
etavec = [1/3, 0.2, 0.6, 0.2];
wvec = [-27/48, 25/48, 25/48, 25/48];

% Gausspoints for line integration
% --------------------------------
ngaussb = 3;
betavec = [-0.774597, 0, 0.774597];
ewvec = [5/9, 8/9, 5/9];

% co-ordinates for elements
% -------------------------
coord = [1, 0.0, 0.0;    % 1st column is element no., 2nd column is x-coordinate, 3rd column is y-coordinate
         2, 0.6, 0.0; 
         3, 0.6, 0.4; 
         4, 0.0, 0.4;
         5, 0.3, 0.0;
         6, 0.6, 0.2;
         7, 0.3, 0.4;
         8, 0.0, 0.2;
         9, 0.3, 0.2];

connect = [1, 1, 2, 4, 5, 9, 8;     % 1st column is element no., other columns are global node no.s of that particular element
           2, 2, 3, 4, 6, 7, 9];
       
       planeflg=1;
[C] = constitutive(E, nu, planeflg);

fp = fopen('Steps_Q2_(b).txt','w');

% Global stiffness matrix and global load vector
% ----------------------------------------------
K = zeros(2*nnode,2*nnode);
F = zeros(2*nnode,1);


% Calculation of element stiffness matrix and element load vector
% ---------------------------------------------------------------
fprintf(fp,'Elemental stiffness matrices and Elemental load vectors:\n');
fprintf(fp,'=======================================================');

for el = 1:nele     
    
    nd1 = connect(el,2);
    nd2 = connect(el,3);
    nd3 = connect(el,4);
    nd4 = connect(el,5);
    nd5 = connect(el,6);
    nd6 = connect(el,7);

    xy = [coord(nd1,2), coord(nd1,3);
          coord(nd2,2), coord(nd2,3);
          coord(nd3,2), coord(nd3,3);
          coord(nd4,2), coord(nd4,3);
          coord(nd5,2), coord(nd5,3);
          coord(nd6,2), coord(nd6,3)];

    vec = [2*nd1-1, 2*nd1, 2*nd2-1, 2*nd2, 2*nd3-1, 2*nd3, 2*nd4-1, 2*nd4, 2*nd5-1, 2*nd5, 2*nd6-1, 2*nd6];       % Global DOF                  
    
    kele = zeros(12,12);
    fele = zeros(12,1);

    for gpd = 1:ngaussd     % loop on gauss points for domain
        
        xi = xivec(gpd);    eta = etavec(gpd);     w = wvec(gpd);
            kele(1:12,1:12) = kele(1:12,1:12) + elestiff(xi, eta, xy, C)*w;
    end

    for gpb = 1:ngaussb     % loop on gauss points for boundary terms
        
        beta = betavec(gpb);       ew = ewvec(gpb);
        
        if el==2
            fele(1:12) = fele(1:12) + gamat(beta, xy, P0, P1, theta, 1)*ew;
        end

    end
    
    fprintf(fp,'\n\nElement %d:\n',el);
    fprintf(fp,'-----------\n');
    fprintf(fp,'\nStiffness matrix:\n');

    for i = 1:12
        for j = 1:12
            fprintf(fp,'%14.4e\t',kele(i,j));
        end
        fprintf(fp,'\n');
    end

    fprintf(fp,'Load vector:\n');
    for i = 1:12
        fprintf(fp,'%14.4e\n',fele(i));
    end
    
    % Assembly using connectivity data
    for ii = 1:12
        for jj = 1:12
            K(vec(ii),vec(jj)) = K(vec(ii),vec(jj)) + kele(ii,jj);
        end
        F(vec(ii)) = F(vec(ii)) + fele(ii);
    end
    
end

fprintf(fp,'\n\nGlobal stiffness matrix:\n');
fprintf(fp,'========================\n\n');
for i = 1:2*nnode
    for j = 1:2*nnode
        fprintf(fp,'%14.4e\t',K(i,j));
    end
    fprintf(fp,'\n');
end

fprintf(fp,'\n\nGlobal load vector:\n');
fprintf(fp,'==================\n\n');
for i = 1:2*nnode
    fprintf(fp,'%14.4e\n',F(i));
end


% Imposition of B.C.
Kreduce = K([3:6,9:14,17:18],[3:6,9:14,17:18]);
Freduce = F([3:6,9:14,17:18]);

% Finding Solution
ureduce = Kreduce\Freduce;
un = [0;0;ureduce(1:4);0;0;ureduce(5:10);0;0;ureduce(11:12)];

% Stresses at Gauss Points
Stress_Gauss = zeros(4,1,3);     % No. of elements x No. of gauss points x 3 stress components
for ele = 1:nele
    vec = connect(ele,:);
    
    nd1 = connect(el,2);
    nd2 = connect(el,3);
    nd3 = connect(el,4);
    nd4 = connect(el,5);
    nd5 = connect(el,6);
    nd6 = connect(el,7);

    xy = [coord(nd1,2), coord(nd1,3);
          coord(nd2,2), coord(nd2,3);
          coord(nd3,2), coord(nd3,3);
          coord(nd4,2), coord(nd4,3);
          coord(nd5,2), coord(nd5,3);
          coord(nd6,2), coord(nd6,3)];
    
    uvec = zeros(12,1);
    for i = 1:6
       uvec(2*i-1) = un(2*vec(i)-1);
       uvec(2*i) = un(2*vec(i));
    end
    ig = 0;
    for i = 1:4
        for j = 1:4
            ig = ig + 1;
            [strs] = stress(xivec(i), etavec(j), xy, C, uvec);
            Stress_Gauss(ele,ig,1:3) = strs;
        end
    end
end


fprintf(fp,'\n\nFinal displacements:\n');
fprintf(fp,'========================\n\n');
for i = 1:2*nnode
    fprintf(fp,'%14.4e\n',un(i));
end

fprintf(fp,'\n\nStresses at Gauss Points\n');
fprintf(fp,'=========================\n\n');
for ele = 1:nele
    fprintf(fp,'Element No.: %d\n',ele);
    
    ig = 0;
    for i = 1:4
        for j = 1:4
            ig = ig + 1;
            fprintf(fp,'xi - %10.3e, eta - %10.3e, stresses - %10.3e   %10.3e   %10.3e\n',xivec(i),etavec(j), Stress_Gauss(ele,ig,1:3));
        end
    end
    fprintf(fp,"\n");
end

fclose(fp);

function [C] = constitutive(E, nu, planeflg)

% This function calculate Constitutive Matrix

% INPUT:
% ======
% E = Young's Modulus
% nu = Poisson's ratio
% planeflg = 1  for plane stress
%            2  for plane strain
%
% OUTPUT:
% =======
% C = Constitutive Matrix

if (planeflg == 1) 
   
    C = [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];
    C = C*E/(1-nu^2);
    
elseif (planeflg == 2)
    
    C = [1-nu, nu, 0; nu, 1-nu, 0; 0, 0, (1-2*nu)/2];
    C = C*E/(1+nu)/(1-2*nu);
    
end
end

function [ke] = elestiff(xi, eta, coord, C)

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

alpha = 1 - xi - eta ;
ddmat = [4*xi - 1, 0, 1 - 4*alpha, 4*eta, -4*eta, 4-8*xi-4*eta;
         0,   4*eta - 1,  1-4*alpha, 4*xi, 4-8*eta-4*xi, -4*xi];

J = ddmat*coord;

R1 = [1, 0, 0, 0;
      0, 0, 0, 1;
      0, 1, 1, 0];
  
R2 = zeros(4,4);  R2(1:2,1:2) = inv(J);  R2(3:4,3:4) = inv(J);

R3 = zeros(4,12);
R3(1,1:2:12) = ddmat(1,1:6);
R3(2,1:2:12) = ddmat(2,1:6);
R3(3,2:2:12) = ddmat(1,1:6);
R3(4,2:2:12) = ddmat(2,1:6);

B = R1*R2*R3;

ke = 0.5*det(J)*(B'*C*B);

end

function [fe] = gamat(beta, coord, P0, P1, theta, edgeno)

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
xi = (1+beta)/2;

if edgeno == 1
    ledge = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N1 = (1-xi)*(1-2*xi); N2 = xi*(2*xi-1); N3 = 0; N4 = 4*xi*(1-xi); N5 = 0; N6 = 0;
elseif edgeno == 2
    ledge = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N1 = 0; N2 = (1-xi)*(1-2*xi); N3 = xi*(2*xi-1); N4 = 0; N5 = 4*xi*(1-xi); N6 = 0;
elseif edgeno == 3
    ledge = sqrt((x(4)-x(3))^2 + (y(4)-y(3))^2);
    N1 = xi*(2*xi-1); N2=0; N3 = (1-xi)*(1-2*xi); N4 = 0; N5 = 0; N6 = 4*xi*(1-xi);
end

N = zeros(2,12);
N(1,1:2:12) = [N1, N2, N3, N4, N5, N6];
N(2,2:2:12) = [N1, N2, N3, N4, N5, N6];



y = N1*y(1) + N2*y(2) + N3*y(3) + N4*y(4) + N5*y(5)+ N6*y(6);

theta = theta*(pi/180);
tvec= [P1*cos(theta) + (P0 - P1)*cos(theta)*y/0.4;
       P1*sin(theta) + (P0 - P1)*sin(theta)*y/0.4];

fe = ledge*N'*tvec*0.5;
  
end

function [strs] = stress(xi, eta, coord, C, uvec)

% This function calculate element Stiffness Matrix

% INPUT:
% ======
% C = Constitutive Matrix
% xi, eta = Gausspoint
% coord = Nodal coordinates of the element
%
% OUTPUT:
% =======
% strs = Stress

x = coord(:,1);
y = coord(:,2);

alpha = 1 - xi - eta ;

ddmat = [4*xi - 1, 0, 1 - 4*alpha, 4*eta, -4*eta, 4-8*xi-4*eta;
         0,   4*eta - 1,  1-4*alpha, 4*xi, 4-8*eta-4*xi, -4*xi];

J = ddmat*coord;

R1 = [1, 0, 0, 0;
      0, 0, 0, 1;
      0, 1, 1, 0];
  
R2 = zeros(4,4);  R2(1:2,1:2) = inv(J);  R2(3:4,3:4) = inv(J);

R3 = zeros(4,12);
R3(1,1:2:12) = ddmat(1,1:6);
R3(2,1:2:12) = ddmat(2,1:6);
R3(3,2:2:12) = ddmat(1,1:6);
R3(4,2:2:12) = ddmat(2,1:6);


B = R1*R2*R3;

strs = C*B*uvec;


end

