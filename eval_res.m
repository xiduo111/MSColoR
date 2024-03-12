function rmse = eval_res(X, Y, t, W, V, A, B)
% calculate RMSE and CC


[n, m] = size(Y{1});
c = length(Y);
 for i = 1 : c
        Y_gene{i} = Y{i} - B(i) - A(i) * t;
 end

 for j = 1:c
    for i = 1 : n
        Y_gene_pre{j}(i, :) = X(i, :) * V(:, j) + X(i, :) * W(:, j) * t ;
        rmse(i, j) = sqrt(norm((Y_gene{j}(i, :)  - Y_gene_pre{j}(i, :)).^2 ) / m); 
    end
end

rmse = mean(rmse(:));







     
     
     
     
     
     

 