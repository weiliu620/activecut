function [ O,O_mat ] = outer( fun, vec1, vec2 )
%[O,O_MAT] = OUTER(FUN,VEC1,VEC2)
%Generalized outer product, like calculating VEC1 * VEC2' but instead of
%using multiplication to combine the elements of VEC1 and VEC2, the
%function provided by function pointer fun is called for each pair of
%elements and the results are stored in the N-by-M cell matrix O. N is the
%length of VEC1 and M is the length of VEC2.
%The function works primarily on cell arrays VEC1 and VEC2, but inputs in
%the form of numeric vectors/matrices are accepted too. When using numeric
%inputs, the function assumes the first dimension of VEC1 and VEC2 to
%contain the elements to combine.
%   FUN:    operator function pointer, the function must be a binary
%           function operating on inputs of the type contained in the cells
%           of VEC1 and VEC2, its output type is arbitrary
%   VEC1:   first cell array or numeric array/matrix
%   VEC2:   second cell array or numeric array/matrix
%   O:      resulting 2D cell array ('cell matrix')
%   O_MAT:  result cast as numeric matrix, NaN if impossible to cast
%
%Author: Philipp Keller
%Date:   2010-11-16
%

% if input is provided as numerical matrix, each element of dimension 1
% will become the content of one cell
if (iscell(vec1))
    cell1 = vec1;
elseif (isnumeric(vec1))
    cell1 = mat2cell(vec1,ones(1,size(vec1,1)));
    cell1 = cellfun(@squeeze,cell1,'uniformoutput',false);
else
    error('OUTER: vec1 needs to be a cell array or a numeric matrix');
end
if (iscell(vec2))
    cell2 = vec2;
elseif (isnumeric(vec2))
    cell2 = mat2cell(vec2,ones(1,size(vec2,1)));
    cell2 = cellfun(@squeeze,cell2,'uniformoutput',false);
else
    error('OUTER: vec2 needs to be a cell array or a numeric matrix');
end
% create horizontal cell arrays from input vectors
if (size(cell1,1) ~= 1)
    cell1 = cell1(:)';
end
if (size(cell2,1) ~= 1)
    cell2 = cell2(:)';
end

% call the provided function inside a nested call of cellfun.
% use anonymous functions for the cellfun's, to be able to directly insert
% the cell1 and cell2 arrays.
O = cellfun(@(q)(cellfun(@(p)fun(p,q),cell1,'uniformoutput',false))', ...
            cell2, 'uniformoutput',false);
% flatten the inner dimension to create a 2D cell array ('cell matrix')
O = [O{:}];
% try to create the numerical matrix
try
    O_mat = cell2mat(O);
catch
    O_mat = NaN;
end

end

