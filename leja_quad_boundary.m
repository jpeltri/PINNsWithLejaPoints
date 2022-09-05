function [pts,wts] = leja_quad_boundary(level,xdims)

knots_t = @(n) knots_leja(n,0,1,'line');    %use leja points on [0,1] for the time
knots_x = @(n) knots_leja(n,-1,1,'line');   %use leja points on [-1,1] for spatial dimensions

%we construct the index set explicitly, it needs to be ordered
%here, we want the index of the spatial dimension to be 1 or 2
if xdims==1
    C = zeros(2*level - 3,2);   %initialization
    for i=1:level-2             %indices of the time dimension
        C(2*i-1,:) = [i 1];     %adding both desired indices
        C(2*i,:) = [i 2];
    end
    C(2*level-3,:) = [level-1 1];   %adding the last single index
    
    %based on this index set the Sparse Grids Matlab Kit can now construct
    %the desired quadrature formula
    S = smolyak_grid_multiidx_set(C,{knots_t,knots_x},@lev2knots_lin);
end

%with two spatial dimensions, we want one of their indices to be 1 or 2
%this requires some complexity since the set shall still be ordered
if xdims==2
    num = 0;
    for i = 1:level-4
        %first the indices with 1 in the first spatial dimensions
        for k = 1:level-1-i
            num = num+1;
            C(num,:) = [i 1 k];
        end
        
        %then those with 2 in the first spatial dimensions
        for k = 1:level-2-i
            num = num+1;
            C(num,:) = [i 2 k];
        end
        
        %then those that are on the boundary in the second spatial
        %dimension
        for j = 3:level-2-i 
            num = num+1;           
            C(num,:) = [i j 1];
            C(num+1,:) = [i j 2];
            num = num+1;
        end

        num = num+1;
        C(num,:) = [i level-1-i 1];
    end
    
    %finally, some more single indices to be added
    C(num+1,:) = [level-3 1 1];
    C(num+2,:) = [level-3 1 2];
    C(num+3,:) = [level-3 2 1];
    C(num+4,:) = [level-2 1 1];


    S = smolyak_grid_multiidx_set(C,{knots_t,knots_x,knots_x},@lev2knots_lin);
end

%in three dimensions, the logic is similar to before but even more involved
if xdims==3
    num = 0;
    for i = 1:level-6
        for k = 1:level-2-i
            for l = 1:level-1-i-k
                num = num+1;
                C(num,:) = [i 1 k l];
            end
        end
        
        for k = 1:level-3-i
            for l = 1:level-2-i-k
                num = num+1;
                C(num,:) = [i 2 k l];
            end
        end

        for j = 3: level-i-4
            for l= 1:level-i-j-1
                num = num+1;
                C(num,:) = [i j 1 l];
            end
            for l= 1:level-i-j-2
                num = num+1;
                C(num,:) = [i j 2 l];
            end
            for k = 3:level-i-j-1
            end 
            for k = 3:level-i-j-2
                num = num+1;
                C(num,:) = [i j k 1];
                num = num+1;
                C(num,:) = [i j k 2];
            end  
            num = num+1;
            C(num,:) = [i j level-i-j-1 1];
        end

        num = num+1;
        C(num,:) = [i level-i-3 1 1];
        num = num+1;
        C(num,:) = [i level-i-3 1 2];
        num = num+1;
        C(num,:) = [i level-i-3 2 1];
        num = num+1;
        C(num,:) = [i level-i-2 1 1];
    end

    num = num+1;
    C(num,:) = [level-5 1 1 1];
    num = num+1;
    C(num,:) = [level-5 1 1 2];
    num = num+1;
    C(num,:) = [level-5 1 1 3];
    num = num+1;
    C(num,:) = [level-5 1 2 1];
    num = num+1;
    C(num,:) = [level-5 1 2 2];
    num = num+1;
    C(num,:) = [level-5 1 3 1];
    num = num+1;
    C(num,:) = [level-5 2 1 1];
    num = num+1;
    C(num,:) = [level-5 2 1 2];
    num = num+1;
    C(num,:) = [level-5 2 2 1];
    num = num+1;
    C(num,:) = [level-5 3 1 1];

    num = num+1;
    C(num,:) = [level-4 1 1 1];
    num = num+1;
    C(num,:) = [level-4 1 1 2];
    num = num+1;
    C(num,:) = [level-4 1 2 1];
    num = num+1;
    C(num,:) = [level-4 2 1 1];
    num = num+1;
    C(num,:) = [level-3 1 1 1];

    S = smolyak_grid_multiidx_set(C,{knots_t,knots_x,knots_x,knots_x},@lev2knots_lin);
end

Sr = reduce_sparse_grid(S);         %remove duplicates for computational efficiency

pts = Sr.knots;                     %output points of the grid
wts = Sr.weights/sum(Sr.weights);   %we scale weights to a sum of one for consistency between quadrature formulas