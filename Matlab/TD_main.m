%**************************************************************************
%*    THE FOLLOWING PROGRAM USES FOURTH-ORDER RUNGE-KUTTA METHOD          *
%*    METHOD TO COMPUTE WIND-INDUCED EFFECTS IN THE TIME-DOMAIN           *
%*                                                                        *
%*                             -*-                                        *
%*                      SEPTEMEBER 24, 2011                               * 
%*                       WORKAMAW WARSIDO                                 *
%**************************************************************************
clear all; 
close all; 
clc;
%--------------------------------------------------------------------------
% Input analysis parameters
%--------------------------------------------------------------------------
display('Choose the numerical scheme to be used');
sc =input('Press "i" for implicit or "e" for explicit scheme --->','s');
while(sc~='i'&& sc~='e')
    sc =input('Wrong character entered press "i" for implicit or "e" for explicit scheme --->','s');
end 
n = input('Enter the number of floors of the building(n) --->');
k = input('Enter the number of modes to be considered for analysis(k) --->');
T = input('Enter the test duration in sec(T) --->');
sf = input('Enter the test sampling frequency in Hz(sf) -->');
N = T*sf;                               %Length of time history data
dt = 1/(sf/100);                        %time interval to be used for numerical analysis obtained using time scale of 1:100
%--------------------------------------------------------------------------
%  DATA READING
%--------------------------------------------------------------------------	
mv=xlsread('mass.xls');                
m=mat_diag(mv,3*n);                    %Diagonal mass matrix (3n x 3n);                  
h=xlsread('height.xls');               %Height vextor (n x 1);
ra=xlsread('radius.xls');              %Torsional radius vector (n x 1);
e=xlsread('centroid.xls');             %Eccentricity of C.M. from origin (n x 2);
Wv=xlsread('omega.xls');
W=mat_diag(Wv,k);                      %Diagonal matrix of omega(W) (k x k);
phi=xlsread('mode.xls');               %Mode shape matrix (3n x k);
SIv=xlsread('damping.xls');            
SI=mat_diag(SIv,k);                    %Diagonal damping matrix (k x k);
W2=W*W;	                               %Diagonal matrix of omega2(W2) (k x k);
%--------------------------------------------------------------------------
% Read load file names
%--------------------------------------------------------------------------
dirData = dir('G:\RWDI\Dissertation work\Chapter 3\ANALYSIS\FL load for diff reduced velocity\70ms');      % Get directory contents 
dirData = dirData(~[dirData.isdir]);  % Use only the file data
fileNames = {dirData.name};           % Get file names 
file_name = char(fileNames);
out_q(:,1:14)=strcat(file_name(:,1:8),'_','q','.xls');
out_u(:,1:14)=strcat(file_name(:,1:8),'_','u','.xls');
out_a(:,1:14)=strcat(file_name(:,1:8),'_','a','.xls');
out_v(:,1:14)=strcat(file_name(:,1:8),'_','v','.xls');
out_m(:,1:14)=strcat(file_name(:,1:8),'_','m','.xls');
out_t(:,1:14)=strcat(file_name(:,1:8),'_','t','.xls');

for I = 1:length(fileNames)                 %computation for each angle done by this loop
%change the directory to where load files are located
cd('G:\RWDI\Dissertation work\Chapter 3\ANALYSIS\FL load for diff reduced velocity\70ms');     
pv=xlsread(file_name(I,:),'F');   
p=pv';                                 %Load matrix  (3n x N);
%change the directory to where the analysis files are loated
cd('G:\RWDI\Dissertation work\Chapter 3\ANALYSIS\Analysis for diff velocity\TD analysis');  
%--------------------------------------------------------------------------
% Mean and fluctuating components of the applied force
%--------------------------------------------------------------------------
avg_p = mean(p,2);                     % mean component
for i=1:3*n
    flu_p(i,:)=p(i,:) - avg_p(i);      % fluctuating component
end
%--------------------------------------------------------------------------
% Mean and fluctuating components of the generalizeg force
%--------------------------------------------------------------------------
avg_P = phi'*avg_p;                    % mean component 
P = phi'*p;                            % fluctuating component
%--------------------------------------------------------------------------
% Modified mode shapes to consider eccentricity 
%--------------------------------------------------------------------------
for i=1:n
    for j=1:k
        phi(i,j)=phi(i,j)-e(i,2)*phi(2*n+i,j);            %x-direction
        phi(n+i,j)=phi(n+i,j)+ e(i,1)*phi(2*n+i,j);       %y-direction
        phi(2*n+i,j)=phi(2*n+i,j);                        %Rotation
    end
end
phit=phi';                     % Transpose of the modified modeshape matrix
%--------------------------------------------------------------------------
% Modal properties
%--------------------------------------------------------------------------
M=phit*m*phi;                  % Modal mass matrix (k x k)
K=W2*M;                        % Modal stiffness matrix (k x k)
C=2*SI*M*W;                    % Modal damping matrix (k x k)
%--------------------------------------------------------------------------
% Modal participation coefficients
%--------------------------------------------------------------------------
[ga_d,ga_a,ga_v,ga_m,ga_t] = gamma(phi,W2,m,h,e,n,k);
%--------------------------------------------------------------------------
% Compute q, qd and qdd
%--------------------------------------------------------------------------
if(sc=='i')
    [q,qd,qdd] =q_implicit(M,P,C,K,k,N,dt);
else
    [q,qd,qdd] =q_explicit(M,P,C,K,k,N,dt);
end
%--------------------------------------------------------------------------
% Statistical values of the generalized coordinate
%--------------------------------------------------------------------------
avg_q = mean(q,2);                 %Average of the generalized coordinate 
std_q = std(q,0,2);                %Stadndard deviation of the generalized coordinate
%--------------------------------------------------------------------------
% Modal participation coefficients
%--------------------------------------------------------------------------
[ga_d,ga_a,ga_v,ga_m,ga_t]=gamma(phi,W2,m,h,e,n,k);
%--------------------------------------------------------------------------
% Displacement response  
%--------------------------------------------------------------------------
[um,u]=td_response(ga_d,q,3*n,k,N);
%Statistics of modal displacements 
avg_um = mean(um,3);
std_um = std(um,0,3);
%Statistics of total dispalcement 
avg_u = mean(u,2);
std_u = std(u,0,2);
max_u = max(u,[],2);
min_u = min(u,[],2);
[max_u1,min_u1] = peak(u,1);
%--------------------------------------------------------------------------
% Acceleration response  
%--------------------------------------------------------------------------
[am,a]=td_response(ga_a,q,3*n,k,N);
%Statistics of modal accelerations 
avg_am = mean(am,3);
std_am = std(am,0,3);
%statistics of total acceleration  
avg_a = mean(a,2);
std_a = std(a,0,2);
max_a = max(a,[],2);
min_a = min(a,[],2);
[max_a1,min_a1] = peak(a,1);
%--------------------------------------------------------------------------
% Storey shear response  
%--------------------------------------------------------------------------
[vm,v]=td_response(ga_v,q,2*n,k,N);
%Statistics of modal shears
avg_vm = mean(vm,3);
std_vm = std(vm,0,3);
%Statistics of total shear
avg_v = mean(v,2);
std_v = std(v,0,2);
max_v = max(v,[],2);
min_v = min(v,[],2);
[max_v1,min_v1] = peak(v,1);
%--------------------------------------------------------------------------
% Storey bending moment response  
%--------------------------------------------------------------------------
[mbm,bm]=td_response(ga_m,q,2*n,k,N);
%Statistics of modal bending moments 
avg_mm = mean(mbm,3);
std_mm = std(mbm,0,3);
%Statistics of total bending moment 
avg_m = mean(bm,2);
std_m = std(bm,0,2);
max_m = max(bm,[],2);
min_m = min(bm,[],2);
[max_m1,min_m1] = peak(bm,1);
%--------------------------------------------------------------------------
% Storey torsion response  
%--------------------------------------------------------------------------
[tm,t]=td_response(ga_t,q,n,k,N);
%Statistics of modal torsional moments 
avg_tm = mean(tm,3);
std_tm = std(tm,0,3);
%Statistics of total torsional moments
avg_t = mean(t,2);
std_t = std(t,0,2);
max_t = max(t,[],2);
min_t = min(t,[],2);
[max_t1,min_t1] = peak(t,1);
%--------------------------------------------------------------------------
% Write the results on Excel files
%--------------------------------------------------------------------------
%Generalized coordinate (q)
XLSWRITE(out_q(I,:),avg_q,'avg_q');
XLSWRITE(out_q(I,:),std_q,'std_q');
%Modal floor displacement
XLSWRITE(out_u(I,:),avg_um,'avg_um');
XLSWRITE(out_u(I,:),std_um,'std_um');
%Total floor displacement
XLSWRITE(out_u(I,:),avg_u,'avg_u');
XLSWRITE(out_u(I,:),std_u,'std_u');
XLSWRITE(out_u(I,:),max_u,'max_u');
XLSWRITE(out_u(I,:),max_u1,'max_u1');
XLSWRITE(out_u(I,:),min_u,'min_u');
XLSWRITE(out_u(I,:),min_u1,'min_u1');
%Modal floor acceleration
XLSWRITE(out_a(I,:),avg_am,'avg_am');
XLSWRITE(out_a(I,:),std_am,'std_am');
%Total floor acceleration
XLSWRITE(out_a(I,:),avg_a,'avg_a');
XLSWRITE(out_a(I,:),std_a,'std_a');
XLSWRITE(out_a(I,:),max_a,'max_a');
XLSWRITE(out_a(I,:),max_a1,'max_a1');
XLSWRITE(out_a(I,:),min_a,'min_a');
XLSWRITE(out_a(I,:),min_a1,'min_a1');
%Modal storey shear
XLSWRITE(out_v(I,:),avg_vm,'avg_vm');
XLSWRITE(out_v(I,:),std_vm,'std_vm');
%Total storey shear
XLSWRITE(out_v(I,:),avg_v,'avg_v');
XLSWRITE(out_v(I,:),std_v,'std_v');
XLSWRITE(out_v(I,:),max_v,'max_v');
XLSWRITE(out_v(I,:),max_v1,'max_v1');
XLSWRITE(out_v(I,:),min_v,'min_v');
XLSWRITE(out_v(I,:),min_v1,'min_v1');
%Modal storey bending moment
XLSWRITE(out_m(I,:),avg_mm,'avg_mm');
XLSWRITE(out_m(I,:),std_mm,'std_mm');
%Total storey bending moment
XLSWRITE(out_m(I,:),avg_m,'avg_m');
XLSWRITE(out_m(I,:),std_m,'std_m');
XLSWRITE(out_m(I,:),max_m,'max_m');
XLSWRITE(out_m(I,:),max_m1,'max_m1');
XLSWRITE(out_m(I,:),min_m,'min_m');
XLSWRITE(out_m(I,:),min_m1,'min_m1');
%Modal storey torsion
XLSWRITE(out_t(I,:),avg_tm,'avg_tm');
XLSWRITE(out_t(I,:),std_tm,'std_tm');
%Total storey torsion
XLSWRITE(out_t(I,:),avg_t,'avg_t');
XLSWRITE(out_t(I,:),std_t,'std_t');
XLSWRITE(out_t(I,:),max_t,'max_t');
XLSWRITE(out_t(I,:),max_t1,'max_t1');
XLSWRITE(out_t(I,:),min_t,'min_t');
XLSWRITE(out_t(I,:),min_t1,'min_t1');

end





