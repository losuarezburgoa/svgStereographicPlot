function [ stereoPlaneStrCell, svgCompleteText ] =svgstereographicplot ...
    ( dipdirDipArray, svgFileNameString, plotStereoTemplateTrue, colorStringCell )
%
% Description:
% Creates a Scalable Vector Graphic (SVG) file of the graphic where the
% orientation of three dimensional planes are represented in an equalangle
% spherical projection. Files with SVG extension can be shown in any Web
% Browser (e.g. Firefox), its edition can be performed in Inkscape program.
% 
% Default values:
% If svgFileNameString is not defined, the program puts the name of
% 'stereographicRepPlanes.svg' to the output file.
% If plotStereoTemplateTrue is not defined, it is by default false.
% If colorStringCell is not defined, the program generates random colors.
%
% External sub-function(s):
% randomsvgcolors.
%
% Input(s):
% Dip direction and dip of the plane is wanted to plot, in an 1x2 array
% (dipdirDipArray).
%
% Name --without extension-- of the SVG file to be created
% (svgFileNameString).
%
% A boolean value of true if the stereographic template is wanted to show
% (plotStereoTemplateTrue). This template is compossed by: the great
% circle, its center, the north, east, south and western ticks, the north
% symbol, and the comments texts.
%
% Cell of texts or hexagesimal numbers that specifies the color of each
% plane. For more color names see for example:
%   http://www.december.com/html/spec/colorsvg.html
% The hexagesimal numbers are given in the following format: "#000000"
% (colorStringCell).
%
% Output(s):
% Cell that stores in each element the structures of the elements and
% variables of each plane, that can be used for other matlab calculations
% (stereoPlaneStrCell).
%  stereoPlaneSTR.planeTrace: stores the arcSTR, where:
%      arcSTR.center: coordiantes of the arc center;
%      arcSTR.firstPt: coordinates of the arc start-point;
%      arcSTR.endPt: coordinates of the arc end-point.
%  stereoPlaneSTR.pole: stores the poleSTR, where:
%      poleSTR.pole: vector of 2x1 of the pole location in ENL coordinate
%      system.
%  stereoPlaneSTR.specialPts: stores the specialPtsSTR structure, where:
%      specialPtsSTR.dipdirPt: coordiantes of the middle point of the arc.
%  stereoPlaneSTR.equation: stores the equationSTR structure, where:
%      equationSTR.polynomCoeffs: coefficients of the polynomial equation
%      of the arc.
%
% Text chain with the structure of graphic in SVG (svgCompleteText).
%
% Example1:
% dipdirDipArray =[ 236, 25; 124, 16; 56, 85; 180, 10; 305, 45; ...
%   25, 90; 190, 0 ];
% svgFileNameString ='example1';
% plotStereoTemplateTrue =true;
% colorStringCell ={ 'blue', 'red', 'green', 'purple', 'sandybrown', ...
%   'darkolivegreen', 'gold' };
%
% %%%%%%%%%%%%%%%%
% svgstereographicplot( dipdirDipArray );
% %%%%%%%%%%%%%%%%

%% % Number of data %
numPlanes =size( dipdirDipArray, 1 );

%% % Input management %
if nargin <2
    svgFileNameString ='stereographicRepPlanes';
    plotStereoTemplateTrue =false;
    colorStringCell =randomsvgcolors( 1, numPlanes );
elseif nargin < 3
    plotStereoTemplateTrue =false;
    colorStringCell =randomsvgcolors( 1, numPlanes );
elseif nargin <4
    colorStringCell =randomsvgcolors( 1, numPlanes );
end

%% % Doing the calculations for each plane %
stereoPlaneStrCell =cell( 1, numPlanes );
svgFileTextCell =cell( 1, numPlanes );

for i=1 :numPlanes
    [ stereoPlaneSTR, svgFileText ] =createplanescripts ...
        ( dipdirDipArray(i,:), i, colorStringCell{i} );
    stereoPlaneStrCell{i} =stereoPlaneSTR;
    svgFileTextCell(i) ={svgFileText};
end

%% % Script for plotting in SVG %
% Create the svg code head %
svgHeadText =[ ...
    sprintf( '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' ), ...
    sprintf( '<svg\n' ), ...
    sprintf( 'xmlns="http://www.w3.org/2000/svg"\n' ), ...
    sprintf( 'viewBox="0 0 2500 2500"\n' ), ...
    sprintf( 'version="1.1"\n' ), ...
    sprintf( 'id="basicStereograph"\n' ), ...
    sprintf( 'width="100%%"\n' ), ...
    sprintf( 'height="100%%">\n\n' ), ...
    sprintf( '<title id="title3015">Stereographic representation</title>\n\n' ), ...
    sprintf( '/* Begins drawing /\n' ), ...
    sprintf( '/* Put here all the plane-groups /\n' ) ];

% Create the svg code for each plane %
svgPlanesText =cell2mat( svgFileTextCell );

% Create the svg code for the stereographic template %
if plotStereoTemplateTrue
    comment1String ='Project name';
    comment2String ='Date';
    % script to draw the template %
    svgStereoTemplateText =[ ...
      sprintf( '/* Here the stereographic background is put /\n' ), ...
      sprintf( '<g id="stereographicPlate" fill="none" stroke="black" >\n' ), ...
      sprintf( '  <path id="greatCircle" d="M1250,1250 m-1000,0 a 1000,1000 0 1,0 2000,0 a 1000,1000 0 1,0 -2000,0" stroke-width="20"/>\n' ), ...
      sprintf( '  <g id="centerCross" stroke-width="8" >\n' ), ...
      sprintf( '    <path id="centerHorz" d="M1250,1250 m-50,0 l100,0" />\n' ), ...
      sprintf( '    <path id="centerVert" d="M1250,1250 m0,-50 l0,100" />\n' ), ...
      sprintf( '  </g>\n' ), ...
      sprintf( '  <path id="northLine" d="M1250,1250 m0,-1080 l0,80" stroke-width="10" />\n' ), ...
      sprintf( '  <path id="eastLine" d="M1250,1250 m1000,0 l80,0" stroke-width="10" />\n' ), ...
      sprintf( '  <path id="soithLine" d="M1250,1250 m0,1000 l0,80" stroke-width="10" />\n' ), ...
      sprintf( '  <path id="westLine" d="M1250,1250 m-1080,0 l80,00" stroke-width="10" />\n' ), ...
      sprintf( '  <text x="1250" y="120" font-family="Sans" font-size="120" fill="black" stroke-width="none" style="dominant-baseline: central; text-anchor: middle;"> N </text>\n' ), ...
      sprintf( '  <text x="0" y="2400" font-family="monospace" font-size="48" fill="black"> Lower hemisphere </text>\n' ), ...
      sprintf( '  <text x="0" y="2450" font-family="monospace" font-size="48" fill="black"> equalangle spherical projection </text>\n' ), ...
      sprintf( '  <text x="2500" y="2400" font-family="monospace" font-size="48" fill="black" text-anchor="end"> %s </text>\n', comment1String ), ...
      sprintf( '  <text x="2500" y="2450" font-family="monospace" font-size="48" fill="black" text-anchor="end"> %s </text>\n', comment2String ), ...
      sprintf( '</g>\n' ), ...
     ];
else
    svgStereoTemplateText ='';
end

% Create tje svg code for the foot %
svgFoodText =[ ...
    sprintf( '/* Ends drawing /\n' ), ...
    sprintf( '</svg>\n' ), ...
   ];

% Concatenate the complete dvg code %
svgCompleteText =[ svgHeadText, svgStereoTemplateText, svgPlanesText, ...
    svgFoodText ];

% Display on screen %
% display( svgCompleteText );

%% % Saving into a file %
fileID =fopen( [svgFileNameString,'.svg'], 'w', 'n', 'UTF-8' );
fprintf( fileID, '%s', svgCompleteText );
status =fclose( fileID );
if ~status
    display('SVG file created, and successfully closed!');
else
    display('File could not be closed!');
end

end
%% FUNCTION
function [ stereoPlaneSTR, svgFileText ] =createplanescripts( dipdirDip, ...
    planeCounter, colorString )
% Input(s):
% Dip direction and dip of the plane is wanted to plot, in an 1x2 array
% (dipdirDipArray).
%
% Plane name id counter, an integer number (planeCounter).
%
% Text or hexagesimal number that specifies the color of the plane
% (colorString). For more color names see for example:
% http://www.december.com/html/spec/colorsvg.html
% The hexagesimal number is given in the following format: "#000000".
%
% Output(s):
% Structure of the elements and variables of the plane (stereoPlaneSTR).
%
% Text chain with the structure of graphic in SVG (svgFileText).
%
% Example1:
% dipdirDipArray =[ 236, 25 ];
% createplanescripts( [ 236, 25 ], 1, 'blue' );
%
% %%%%%%%%%%%%%%%%
% createplanescripts( dipdirDip, planeCounter, colorString )
% %%%%%%%%%%%%%%%%

%% Making difference between an arc and a line in the projection

if and(dipdirDip(2) >=0, dipdirDip(2)<90 )
    [ stereoPlaneSTR, svgFileText ] =circlearcasplane( dipdirDip, ...
        planeCounter, colorString );
elseif dipdirDip(2) ==90
    [ stereoPlaneSTR, svgFileText ] =linesegmentasplane( dipdirDip, ...
        planeCounter, colorString );
else
    error('Dip angle can not be negative in this semiespherical projection!');
end

end

function [ stereoPlaneSTR, svgFileText ] =circlearcasplane( dipdirDip, ...
    planeCounter, colorString )
%% Calculations
% dipdirDeg =dipdirDip(1);
% dipDeg =dipdirDip(2);

% Plane pole calculation %
trendDeg =mod( (dipdirDip(1) +180), 360 );
plungeDeg =90 -dipdirDip(2);
trendPlungeArray =[ trendDeg, plungeDeg ];

% all into radians %
dipdirArrayRad =dipdirDip *pi/180;
trendPlungeArrayRad =trendPlungeArray *pi/180;

% Center and radius of the plane trace %
Rg =1 *sec( dipdirArrayRad(2) );
rg =1 *tan( dipdirArrayRad(2) );

% Radial longitude from net center to the complement angle of the plunge %
rp =1 *tan( (pi/2 -trendPlungeArrayRad(2))/2 );

%% % Special points given in ENL coordinate system %
% ENL coordiante system: x=East, y=North, z=eLevation. %
transMat1 =[ sin(dipdirArrayRad(1)), -cos(dipdirArrayRad(1)); ...
    cos(dipdirArrayRad(1)), sin(dipdirArrayRad(1)) ];

% Point at initial arc in ENL coordiante system %
p1Vec =transMat1 *[0, -1]';
% Point at end arc in ENL coordiante system %
p2Vec =transMat1 *[0, 1]';
% Point at middle arc in ENL coordiante system %
pmedVec =transMat1 *[(Rg-rg), 0]';
% Point of the pole %
poleVec =transMat1 *[-rp, 0]';

%% % Circle equation in ENL coordinate system %
% canonical equation parameters %
xo =0;
yo =0;
r =Rg;

% second order polynomial equation %
% A*1 +B*x +C*y +D*x^2 +E*xy +F*y^2 =0 %
A =xo^2 +yo^2 -r^2;
B =-2*xo;
C =-2*yo;
D =1;
E =0;
F =1;
coeffVec =[ A, B, C, D, E, F ];

%% % Creating an object structure %
arcSTR =struct( 'center', [0, 0]', 'firstPt', p1Vec, 'endPt', p2Vec );
poleSTR =struct( 'pole', poleVec );
specialPtsSTR =struct( 'dipdirPt', pmedVec );
equationSTR =struct( 'polynomCoeffs', coeffVec );

stereoPlaneSTR =struct( 'planeTrace', arcSTR, 'pole', poleSTR, ...
    'specialPts', specialPtsSTR, 'equation', equationSTR );

%% % Script for plotting in SVG %
% drawing the arc and the pole %
rotAngleDeg =90 -dipdirDip(1);
squareSize =30;
squareFillString ='white';

svgFileText =[ ...
    sprintf('<g id="plane%d" fill="none" stroke="%s" transform="rotate(%7.4f 1250 1250)">\n', ...
      planeCounter, colorString, -rotAngleDeg ), ...
    sprintf('  <path id="trace" d="M1250,1250 m0,-1000 a%4.1f,%4.1f 0 0,1 0,2000" stroke-width="10" />\n', ...
      Rg*1000, Rg*1000 ), ...
    sprintf('  <path id="midTick" d="M1250,1250 m%4.1f,0 l80,0" stroke-width="5" />\n', ...
      (Rg-rg)*1000 ), ...
    sprintf('  <rect id="rectPole" transform="rotate(%7.4f %4.1f 1250)" x="%4.1f" y="%d" width="%d" height="%d" fill="%s" stroke-width="5" />\n', ...
      rotAngleDeg, 1250-rp*1000, 1250-rp*1000-squareSize/2, 1250-squareSize/2, squareSize, squareSize, squareFillString ), ...
    sprintf('  <text id="poleNum" transform="rotate(%7.4f %4.1f 1250)" x="%4.1f" y="%d" font-family="monospace" font-size="48" fill="%s" stroke-width="none" style="dominant-baseline: central; text-anchor: end;"> %s </text>\n', ...
      rotAngleDeg, 1250-rp*1000, 1250-rp*1000-1.1*squareSize, 1250-squareSize/2, colorString, num2str(planeCounter) ), ...
    sprintf('</g>\n\n') ];
end
%% FUNCTION
function [ stereoPlaneSTR, svgFileText ] =linesegmentasplane( dipdirDip, ...
    planeCounter, colorString )
%% Calculations
% dipdirDeg =dipdirDip(1);
% dipDeg =dipdirDip(2);

% all into radians %
dipdirArrayRad =dipdirDip *pi/180;

%% % Special points given in ENL coordinate system %
% ENL coordiante system: x=East, y=North, z=eLevation. %
transMat1 =[ sin(dipdirArrayRad(1)), -cos(dipdirArrayRad(1)); ...
    cos(dipdirArrayRad(1)), sin(dipdirArrayRad(1)) ];

% Point at initial arc in ENL coordiante system %
p1Vec =transMat1 *[0, -1]';
% Point at end arc in ENL coordiante system %
p2Vec =transMat1 *[0, 1]';
% Point at middle arc in ENL coordiante system %
pmedVec =transMat1 *[0, 0]';
% Point of the pole %
poleVec =transMat1 *[-1, 0]';
% Center %
centerVec =[-Inf, 0]';

%% % Line equation in ENL coordinate system %
% canonical equation parameters %
xo =0;
yo =0;
m =1 /tan( dipdirDip(1)*180/pi );

% first order polynomial equation %
% A*1 +B*x +C*y =0 %
A =yo -m*xo;
B =-m;
C =1;
coeffVec =[ A, B, C ];

%% % Creating an object structure %
arcSTR =struct( 'center', centerVec, 'firstPt', p1Vec, 'endPt', p2Vec );
poleSTR =struct( 'pole', poleVec );
specialPtsSTR =struct( 'dipdirPt', pmedVec );
equationSTR =struct( 'polynomCoeffs', coeffVec );

stereoPlaneSTR =struct( 'planeTrace', arcSTR, 'pole', poleSTR, ...
    'specialPts', specialPtsSTR, 'equation', equationSTR );

%% % Script for plotting in SVG %
% drawing the line segment and the pole %
rotAngleDeg =90 -dipdirDip(1);
squareSize =30;
squareFillString ='white';

svgFileText =[ ...
    sprintf('<g id="plane%d" fill="none" stroke="%s" transform="rotate(%7.4f 1250 1250)">\n', ...
      planeCounter, colorString, -rotAngleDeg ), ...
    sprintf('  <path id="trace" d="M1250,1250 m0,-1000 l0,2000" stroke-width="10" />\n' ), ...
    sprintf('  <path id="midTick" d="M1250,1250 m0,0 l80,0" stroke-width="5" />\n' ), ...
    sprintf('  <rect id="rectPole" transform="rotate(%7.4f 250 1250)" x="%d" y="%d" width="%d" height="%d" fill="%s" stroke-width="5" />\n', ...
      rotAngleDeg, 250-squareSize/2, 1250-squareSize/2, squareSize, squareSize, squareFillString ), ...
    sprintf('  <text id="poleNum" transform="rotate(%7.4f %4.1f 1250)" x="%4.1f" y="%d" font-family="monospace" font-size="48" fill="%s" stroke-width="none" style="dominant-baseline: central; text-anchor: end;"> %s </text>\n', ...
      rotAngleDeg, 1250-1000, 1250-1000-1.1*squareSize, 1250-squareSize/2, colorString, num2str(planeCounter) ), ...
    sprintf('</g>\n\n') ];
end
% Written by: Ludger O. Suarez-Burgoa, Assistant Professor.
% Universidad Nacional de Colombia, Medellin (See license.txt file).
%
% Copyright (c) 2014 Universidad Nacional de Colombia
% Copyright (c) 2014 Ludger Suarez-Burgoa
% All rights reserved.
%
% The svgstereographicplot function is free software:
% you can redistribute it and/or modify it under the terms of 
% the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details (file: gpl-3.0.txt).
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
