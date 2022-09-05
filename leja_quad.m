function [pts,wts] = leja_quad(level,xdims)

knots_t = @(n) knots_leja(n,0,1,'line');    %use leja points on [0,1] for the time
knots_x = @(n) knots_leja(n,-1,1,'line');   %use leja points on [-1,1] for spatial dimensions

%the grid can now be constructed completely by the Sparse Grids Matlab Kit
if xdims==1
    S = smolyak_grid(2,level,{knots_t,knots_x},@lev2knots_lin);     %grid with one spatial dimension
end

if xdims==2
    S = smolyak_grid(3,level,{knots_t,knots_x,knots_x},@lev2knots_lin); %grid with two spatial dimensions
end

if xdims==3
    S = smolyak_grid(4,level,{knots_t,knots_x,knots_x,knots_x},@lev2knots_lin); %grid with three spatial dimensions
end

Sr = reduce_sparse_grid(S);         %remove duplicates for computational efficiency

pts = Sr.knots;                     %output points of the grid
wts = Sr.weights/sum(Sr.weights);   %we scale weights to a sum of one for consistency between quadrature formulas