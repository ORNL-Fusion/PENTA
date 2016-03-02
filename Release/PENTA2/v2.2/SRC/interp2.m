function zi = interp2(varargin)
%INTERP2 2-D interpolation (table lookup).
%   Copyright 1984-2008 The MathWorks, Inc.
%   $Revision: 5.33.4.17 $
%
%   Removed lines for speed: 7/2009 JL

ExtrapVal = varargin{end}; % user specified ExtrapVal

 
x=varargin{1};
y=varargin{2};
z=varargin{3};
xi=varargin{4};
yi=varargin{5};


xx = x(1,:); yy = y(:,1);

%
% Check for non-equally spaced data.  If so, map (x,y) and
% (xi,yi) to matrix (row,col) coordinate system.
%
xx = xx.'; % Make sure it's a column.
dx = diff(xx); dy = diff(yy);
xdiff = max(abs(diff(dx))); 
ydiff = max(abs(diff(dy))); 


if (xdiff > eps(class(xx))*max(abs(xx))) || (ydiff > eps(class(yy))*max(abs(yy)))

    % Determine the nearest location of xi in x
    [xxi,j] = sort(xi(:));
    [ignore,i] = sort([xx;xxi]);
    ui(i) = 1:length(i);
    ui = (ui(length(xx)+1:end)-(1:length(xxi)))';
    ui(j) = ui;
    
    % Map values in xi to index offset (ui) via linear interpolation
    ui(ui<1) = 1;
    ui(ui>length(xx)-1) = length(xx)-1;
    ui = ui + (xi(:)-xx(ui))./(xx(ui+1)-xx(ui));
    
    % Determine the nearest location of yi in y
    [yyi,j] = sort(yi(:));
    [ignore,i] = sort([yy;yyi(:)]);
    vi(i) = 1:length(i);
    vi = (vi(length(yy)+1:end)-(1:length(yyi)))';
    vi(j) = vi;
    
    % Map values in yi to index offset (vi) via linear interpolation
    vi(vi<1) = 1;
    vi(vi>length(yy)-1) = length(yy)-1;
    vi = vi + (yi(:)-yy(vi))./(yy(vi+1)-yy(vi));
    
    [x,y] = meshgrid(ones(class(x)):size(x,2),ones(class(y)):size(y,1));
    xi(:) = ui; yi(:) = vi;
    
end


% Now do the interpolation based on method.
zi = linear(ExtrapVal,x,y,z,xi,yi);



%------------------------------------------------------
function F = linear(ExtrapVal,arg1,arg2,arg3,arg4,arg5)
%LINEAR 2-D bilinear data interpolation.
%   ZI = LINEAR(EXTRAPVAL,X,Y,Z,XI,YI) uses bilinear interpolation to
%   find ZI, the values of the underlying 2-D function in Z at the points
%   in matrices XI and YI.  Matrices X and Y specify the points at which
%   the data Z is given.  X and Y can also be vectors specifying the
%   abscissae for the matrix Z as for MESHGRID. In both cases, X
%   and Y must be equally spaced and monotonic.
%
%   Values of EXTRAPVAL are returned in ZI for values of XI and YI that are
%   outside of the range of X and Y.
%
%   If XI and YI are vectors, LINEAR returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%
%   ZI = LINEAR(EXTRAPVAL,Z,XI,YI) assumes X = 1:N and Y = 1:M, where
%   [M,N] = SIZE(Z).
%
%   ZI = LINEAR(EXTRAPVAL,Z,NTIMES) returns the matrix Z expanded by
%   interleaving bilinear interpolates between every element, working
%   recursively for NTIMES. LINEAR(EXTRAPVAL,Z) is the same as
%   LINEAR(EXTRAPVAL,Z,1).
%
%   See also INTERP2, CUBIC.


% linear(extrapval,x,y,z,s,t), X and Y specified.
[nrows,ncols] = size(arg3);
s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols));
if ~isempty(sout), s(sout) = 1; end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows));
if ~isempty(tout), t(tout) = 1; end

% Matrix element indexing
ndx = floor(t)+floor(s-1)*nrows;

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Now interpolate.
onemt = 1-t;

F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
    ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;


% Now set out of range values to ExtrapVal.
if ~isempty(sout), F(sout) = ExtrapVal; end
if ~isempty(tout), F(tout) = ExtrapVal; end
