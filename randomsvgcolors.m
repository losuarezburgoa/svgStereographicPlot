function [ colorStringCell ] =randomsvgcolors( varargin )
%
% Description:
% Generate random colors expressed in hexagesimal format used in SVG
% drawings. 
%
% Input(s):
% Size of the cell to generate (varargin). In case of only one random
% number put it as (1,1).
%
% Output(s):
% Cell of random colors in strings (colorStringCell).
%
% Example1:
% Create a cell of random colors of size (1,2,3).
% colorStringCell =randomsvgcolors( 1, 2, 3 );
%
%%%%%%%%%%%%%%%%%
% [ colorStringCell ] =randomsvgcolors( varargin )
%%%%%%%%%%%%%%%%%

numDims =size(varargin,2);
numValsPerDim =zeros( 1, numDims );

for i=1 :numDims
    numValsPerDim(i) =varargin{i};
end
numElements =prod( numValsPerDim );

colorStringCell =cell( 1, numElements );
for i=1 :numElements
    colorString =['#', randsample('0123456789abcdef',6,true)];
    colorStringCell(i) ={colorString};
end

xx ='';
for i=1 :numDims-1
    xx =[xx, [ num2str(numValsPerDim(i),'%d'), ', ']];
end
xx =[xx, num2str(numValsPerDim(end),'%d')];

evalString =['colorStringCell =reshape( colorStringCell, ', xx, ' );'];
eval( evalString );

end
% Written by: Ludger O. Suarez-Burgoa, Assistant Professor.
% Universidad Nacional de Colombia, Medellin (See license.txt file).
%
% Copyright (c) 2014 Universidad Nacional de Colombia
% Copyright (c) 2014 Ludger Suarez-Burgoa
% All rights reserved.
%
% The randomsvgcolors function is free software:
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
