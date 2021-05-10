%%
% Thermal conductivity is W/(m*C)
% Mass density is 3 kg/m^3
% Specific heat is 0.3 J/(kg*C)
clear
clc

% dX=10; %deltaX
% dY=10; %deltaY
% dZ=2; %deltaZ
% rSource = 15; %jari titik source

dX=25; %deltaX
dY=25; %deltaY
dZ=5; %deltaZ
rSource = 20; %jari titik source


sourceX = 745650;
sourceY = 153059;
sourceZ = 100;

lengthX=1000;
lengthY=1000;
lengthZ=200;

boundaryX = [sourceX-lengthX/2 sourceX+lengthX/2]; clear lengthX;
boundaryY = [sourceY-lengthY/2 sourceY+lengthY/2]; clear lengthY;
boundaryZ = [sourceZ-lengthZ/2 sourceZ+lengthZ/2]; clear lengthZ;

% 1. X
% 2. Y
% 3. Depth
% 4. PHIE
% 5. Permeability
% 6. RHOB
% 7. Young Modulus
% 8. Bulk
% 9. Shear
% 10. Poison Ratio
% 11. Overburden
% 12. Porepressure
% 13. SHmin
% 14. SHmin100_Arah
% 15. SHmax
% 16. Thermal Conductivity
data=textread('/Drive/D/Work/Partial Differential Equation/DataExample/DataChevron#Area1');
xDat = data(:,1);
yDat = data(:,2);
zDat = data(:,3);
poisonRatio = data(:,10);
density = data(:,6).*1000;
youngModulus = data(:,7);

% Create a 3-D rectangular mesh grid
[xg, yg, zg] = meshgrid(boundaryX(1):dX:boundaryX(2),boundaryY(1):dY:boundaryY(2),boundaryZ(1):dZ:boundaryZ(2));
Pcube = [xg(:), yg(:), zg(:)]; clear xg yg zg;

% Extract the grid points located outside of the unit spherical region.
% ukur jadi2 dari titik pusat dengan euclidean
dist = [Pcube(:,1)-sourceX Pcube(:,2)-sourceY Pcube(:,3)-sourceZ];
Pcavitycube = Pcube(vecnorm(dist') > rSource,:); clear Pcube;

% Create points on the unit sphere
[x1,y1,z1] = sphere(24);
Psphere = [x1(:) y1(:) z1(:)];
Psphere = unique(Psphere,'rows');

Psphere(:) = Psphere(:).*rSource;
Psphere(:,1) = Psphere(:,1)+sourceX;
Psphere(:,2) = Psphere(:,2)+sourceY;
Psphere(:,3) = Psphere(:,3)+sourceZ;

% Combine the coordinates of the rectangular grid (without the points inside the sphere) and the surface coordinates of the unit sphere
Pcombined = [Pcavitycube;Psphere];

% Create an alphaShape object representing the cube with the spherical cavity
shpCubeWithSphericalCavity = alphaShape(Pcombined(:,1), ...
                                        Pcombined(:,2), ...
                                        Pcombined(:,3));

% figure
% plot(shpCubeWithSphericalCavity,'FaceAlpha',0.4)
% title('alphaShape: Cube with Spherical Cavity')

% Recover the triangulation that defines the domain of the alphaShape object
[tri,loc] = alphaTriangulation(shpCubeWithSphericalCavity);

% Create a PDE model
modelCube = createpde('thermal','transient');

% Create a geometry from the mesh and import the geometry and the mesh into the model
[gCube,mshCube] = geometryFromMesh(modelCube,loc',tri');

% Plot the resulting geometry.
figure
pdegplot(modelCube,'FaceAlpha',0.5,'FaceLabels','on')
title('PDEModel: Cube with Spherical Cavity')

% generateMesh(modelCube,'Hmax',50,'Hgrad',2,'GeometricOrder','quadratic');
% figure();pdemesh(modelCube,'FaceAlpha',0.5);

%Define thermal material properties.
F1=scatteredInterpolant(xDat,yDat,zDat,thermalConductivity);
k=@(location,state) F1(location.x, location.y, location.z);
F2=scatteredInterpolant(xDat,yDat,zDat,density);
rho=@(location,state) F2(location.x, location.y, location.z);


% cp = 920;
% rho = 2700;
% k = 210;
thermalProperties(modelCube,'ThermalConductivity',k,'SpecificHeat',heatCapacity,'MassDensity',rho)

%Define boundary conditions
thermalBC(modelCube,'Face',7,'Temperature',1000)
% Set IC as 25
h = @(location,state) 25+25.*location.z./50;
thermalIC(modelCube,h)


tlist = 0:2592000:155520000;
R = solve(modelCube,tlist);
figure
pdeplot3D(modelCube,'ColorMapData',R.Temperature(:,end))

