function [pts,wts] = leja_quad_initial(level,xdims)

knots_t = @(n) knots_leja(n,0,1,'line');    %use leja points on [0,1] for the time
knots_x = @(n) knots_leja(n,-1,1,'line');   %use leja points on [-1,1] for spatial dimensions

%we construct the index set explicitly, it needs to be ordered
%here, we want the index in time to always be 2
%because this is not accepted by the Sparse Grids Matlab Kit we instead set
%it to 1 and fix the relevant values later
if xdims==1
    num = 0;
    for j = 1:level-2
        num = num+1;
        C(num,:) = [1 j];   %all spatial indices needed
    end
    %based on this index set the Sparse Grids Matlab Kit can now construct
    %the desired quadrature formula
    S = smolyak_grid_multiidx_set(C,{knots_t,knots_x},@lev2knots_lin);
end

%with 2 and 3 spatial dimensions the procedure is essentially the same
if xdims==2
    num = 0;
    for j = 1:level-3
        for k = 1:level-2-j
            num = num+1;
            C(num,:) = [1 j k];
        end
    end
    S  = smolyak_grid_multiidx_set(C,{knots_t,knots_x,knots_x},@lev2knots_lin);
end

if xdims==3
    num = 0;
    for j = 1:level-4
        for k = 1:level-3-j
            for l = 1:level-2-j-k
                num = num+1;
                C(num,:) = [1 j k l];
            end
        end
    end
    S  = smolyak_grid_multiidx_set(C,{knots_t,knots_x,knots_x,knots_x},@lev2knots_lin);
end

Sr = reduce_sparse_grid(S);         %remove duplicates for computational efficiency

pts = Sr.knots;                     %output points of the grid
wts = Sr.weights/sum(Sr.weights);   %we scale weights to a sum of one for consistency between quadrature formulas
pts(1,:) = zeros(num,1);            %fixing the fact that the index was 2 rather than 1
                                    %by explicitly setting the values of t to 0