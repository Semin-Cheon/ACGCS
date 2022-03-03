%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                    %%%
%%%      ACGCS (Atom Coordinates Generator for Crystal Structure       %%%
%%%                                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2021/10/04 (Y/M/D)

% This script generates atom positions based on primitive vectors of the
% unit cell of your crystal material.

% This script is a modifided version of this script:'create_FePt_model_sized.m'
% which is Written by Yongsoo Yang, the Principal invesgator of MDAIL.

% Author: Semin Cheon. MD/Ph.D integrated course. MDAIL(KAIST). chsm0338@kaist.ac.kr

close all
clear
clc
check(1);

%% User defined parameters
% This is console sector. you can manipulated Basic informations of your
% bulk material, growth number etc...

% ParticleRadius = 1.2; % Bulk material- Particle's Radius in Nanometers
alpha = 1.79;
ParticleRadius = 1.2 * alpha/2; % Particle Radius in Nanometers

lattice_para = 3.87; % FCC lattice parameter in Anstrom

ShiftPara = [0; 0; 0];% shift amount in Angstrom (to break symmetry)

Structure_type = 1;% Select a Structure type: 1 - FCC
%                                             2 - BCC
%                                             3 - Hexagonal 
%                                             4 - Graphene

N = 20; % # of growth sequence

PLOT_YN = 1; % turn on / off the final report
%% Dummy space
P = []; % Dummy space for atom positions ex) P = [P(1), P(2), P(3), ...]
PP = [];

%% Set seed points and primitive vectors

if  Structure_type == 1 % FCC lattice structure (4 atoms per unit cell)
    % FCC lattice structure (4 atoms per unit cell)
    % Primitive vectors (x,y,z)
    a1 = [1,1,0]*lattice_para/2; %(xy plane)
    a2 = [0,1,1]*lattice_para/2; % (yz plane)
    a3 = [1,0,1]*lattice_para/2; % (zx plane)
    zerovec = [0,0,0]*lattice_para/2; % Origin
    Structure_type = sprintf('FCC');
    
elseif Structure_type == 2 % BCC lattice structure (2 atoms per unit cell)
    a1 = [1,1,-1]*lattice_para/2; %
    a2 = [-1,1,1]*lattice_para/2; % 
    a3 = [1,-1,1]*lattice_para/2; % 
    zerovec = [0,0,0]*lattice_para; % Origin
    Structure_type = sprintf('BCC');
    
elseif Structure_type == 3 % Hexagonal system
    
    a1 = [cosd(30),sind(30),0]*lattice_para/2;
    a2 = [-1*sind(60),cosd(60),0]*lattice_para/2;
    a3 = [0 0 1]*lattice_para;
    zerovec = [0,0,0]*lattice_para;
    Structure_type = sprintf('Hexagonal system');
    
elseif Structure_type == 4 % Graphene
    % Primitive vectors (x,y,z)
    a1 = [-1*sind(60),cosd(60), 0]*sqrt(3)*lattice_para/2;
    a2 = [0, 1, 0]*sqrt(3)*lattice_para/2;
    a3 = [0 0 0]*lattice_para/2;
    zerovec = [0,0,0]*lattice_para/2;
    b1 = [cosd(30), sind(30) , 0]*sqrt(3)*lattice_para/2;
    b2 = [0,1 ,0 ]*sqrt(3)*lattice_para/2;
    b3 = [0,0,0]*sqrt(3)*lattice_para/2;
    % the other Seed point
    bb1 = [-1/sqrt(3)+cosd(30), sind(30) , 0]*sqrt(3)*lattice_para/2;
    bb2 = [-1/sqrt(3),1 ,0 ]*sqrt(3)*lattice_para/2;
    bb3 = [-1/sqrt(3),0,0]*sqrt(3)*lattice_para/2;
    
    PP(:,1) = bb1';
    PP(:,2) = bb2';
    PP(:,3) = bb3';
    Structure_type = sprintf('Graphene');
end
check(1);
pause(1)
%% Starts growth sequence

P(:,1) = zerovec'; %Seed positions
LengthofP = 0; % Initial length of the old Ps
LengthofPP = 0; % Initial length of the old PPs This is only used for b's

if Structure_type == 4
for i= 1 : N
    
    L = size(P,2);
    LL = size(PP,2);
    
    Bowl1(:,:) = a1' + P(:,(LengthofP+1:L));
    Bowl2(:,:) = a2' + P(:,(LengthofP+1:L));
    Bowl3(:,:) = a3' + P(:,(LengthofP+1:L));
    Bowl4(:,:) = -a1' + P(:,(LengthofP+1:L));
    Bowl5(:,:) = -a2' + P(:,(LengthofP+1:L));
    Bowl6(:,:) = -a3' + P(:,(LengthofP+1:L));
    
    %%%% For Graphene
    Bowl7(:,:) = b1' + PP(:,(LengthofPP+1:LL));
    Bowl8(:,:) = -b1' + PP(:,(LengthofPP+1:LL));
    Bowl9(:,:) = b2' + PP(:,(LengthofPP+1:LL));
    Bowl10(:,:) = -b2' + PP(:,(LengthofPP+1:LL));
    Bowl11(:,:) = b3' + PP(:,(LengthofPP+1:LL));
    Bowl12(:,:) = -b3' + PP(:,(LengthofPP+1:LL));
   
    KingBowl = [Bowl1,Bowl2,Bowl3,Bowl4,Bowl5,Bowl6];
    Dummy = unique(KingBowl','rows','stable')'; % deleting the overlaped postions
    
    small = [ Bowl7, Bowl8, Bowl9, Bowl10, Bowl11,Bowl12];
    smalldummy = unique(small','rows','stable')'; % deleting the overlaped postions
    
    JJ = size(smalldummy(:,:),2); % the # of the new positions
    LengthofPP = size(PP,2); % length of the old PPs
    
    J = length(Dummy(:,:)); % the # of the new positions
    LengthofP = size(P,2); % length of the old Ps
    
    for j = 1:J
        P(:,L+j) = Dummy(:,j);
    end
    
    for j = 1:JJ
        PP(:,LL+j) = smalldummy(:,j);
    end
    
    clear Bowl7 Bowl8 Bowl9 Bowl10 Bowl11 Bowl12
    
    clear Bowl1 Bowl2 Bowl3
    clear Bowl4 Bowl5 Bowl6
end

else
    for i= 1 : N
    
    L = size(P,2);
    
    Bowl1(:,:) = a1' + P(:,(LengthofP+1:L));
    Bowl2(:,:) = a2' + P(:,(LengthofP+1:L));
    Bowl3(:,:) = a3' + P(:,(LengthofP+1:L));
    Bowl4(:,:) = -a1' + P(:,(LengthofP+1:L));
    Bowl5(:,:) = -a2' + P(:,(LengthofP+1:L));
    Bowl6(:,:) = -a3' + P(:,(LengthofP+1:L));
    
    KingBowl = [Bowl1,Bowl2,Bowl3,Bowl4,Bowl5,Bowl6];
    Dummy = unique(KingBowl','rows','stable')'; % deleting the overlaped postions
 
    J = length(Dummy(:,:)); % the # of the new positions
    LengthofP = size(P,2); % length of the old Ps
    
    for j = 1:J
        P(:,L+j) = Dummy(:,j);
    end
    
    clear Bowl1 Bowl2 Bowl3
    clear Bowl4 Bowl5 Bowl6
    end
end
clear Dummy KingBowl i j J L N 
check(1);

%%
P = round(P*10^(6))/10^(6); % to make elements ~10^(-12) to zero which the error is originated from '부동소수점 계산'
PP = round(PP*10^(6))/10^(6);
P = unique(P','rows','stable')' ; % deleting the overlaped postions
PP = unique(PP','rows','stable')' ; % deleting the overlaped postions

%% Apply Rotation matrix to Atom coordinates

vec1 = [0 0 1]; % X
vec2 = [0 1 0]; % Y
vec3 = [1 0 0]; % Z

if Control_On ==0
    th = 0 ; 
    ph = 0 ; 
    ps = 0 ;
else
    ph = 0 ; 
    ps = 0 ;
end

R = [cosd(ps)*cosd(th),cosd(ps)*sind(th)*sind(ph)-cosd(ph)*sind(ps), sind(ps)*sind(ph)+cosd(ps)*cosd(ph)*sind(th)
        cosd(th)*sind(ps), cosd(ps)*cosd(ph)+sind(ps)*sind(th)*sind(ph), cosd(ph)*sind(ps)*sind(th)-cosd(ps)*sind(ph)
-sind(th),cosd(th)*sind(ph), cosd(th)*cosd(ph) ];

P1 = R*[P,PP]; % Final atoms coordinates

%% User defined area
% make hollow space
% Norms = sqrt( sum(P.^2,1));
% Inside = find(Norms < 4*c);
% P(:,Inside) = [];
%
Norms = sqrt( sum(P1.^2,1)); 
ParticleInd  = find(Norms < ParticleRadius*10);
atomcoordinates = P1(:,ParticleInd);

% % Rantangular shape
% ParticleIndX  = find(P1(1,:).^2 < 100);
% ParticleIndY  = find(P1(2,:).^2 < 100);
% ParticleIndZ  = find(P1(3,:).^2 < 300);
% ParticleInd = intersect(ParticleIndX, ParticleIndY);
% ParticleInd = intersect(ParticleInd,ParticleIndZ);
% 
% atomcoordinates = P1(:,ParticleInd);
% floor = atomcoordinates(3,:);
% floor = unique(floor);
% floor = length(floor);

check(1);
%% set atom types

NumberofPositions = length(P1);
Atomtypes = ones(1,NumberofPositions); %(Default is designated to Fe atom)
atomtype = Atomtypes(ParticleInd);

D = randi([0 1],1,length(atomtype));
atomtype = atomtype + round(D);

% save('AtomicStructure_Simulated_FePt.mat','atomcoordinates','atomtype');
% save(sprintf('AtomicStructure_Simulated_FePt %1.f .mat',floor),'atomcoordinates','atomtype');

if PLOT_YN ==1
check(1);
fprintf('\n Atoms coordinates are generated succesfully. \n %d atom positons are generated.',size(ParticleInd,2))
fprintf('\n\n Crystal structure type is %s',Structure_type)
fprintf('\n Lattice parameter = %.3f A \n Particle radius = %.3f nm \n',lattice_para,ParticleRadius)
fprintf('\n please continue to the next step - "generate_XYZ_forPrismatic.m" for Multislice simulation.\n')
fprintf('                      Or the step - "obtain_3D_volume_testScript.m" to make a 3D volume.\n\n')
check(1);
figure(1);
scatter3(P1(1,:),P1(2,:),P1(3,:),'filled','MarkerFaceColor','g','MarkerEdgeColor','k')
xlim([-50 50])
ylim([-50 50])
zlim([-50 50])
title('Crystal you growth :)')
check(1);
figure(2);
scatter3(atomcoordinates(1,:),atomcoordinates(2,:),atomcoordinates(3,:),'filled','MarkerEdgeColor',[0 0 0])
xlim([-1*ParticleRadius*10 ParticleRadius*10])
ylim([-1*ParticleRadius*10 ParticleRadius*10])
zlim([-1*ParticleRadius*10 ParticleRadius*10])
title(sprintf('Coordinates of Atoms: %d atoms',size(ParticleInd,2)))
end