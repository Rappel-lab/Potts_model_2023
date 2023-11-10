function outputvar_cell = divide_vars_single(inputvar_cell, cellid)
    outputvar_cell = cell(size(inputvar_cell));
    
    for k = 1:length(inputvar_cell)
        outputvar_cell{k} = [inputvar_cell{k}, inputvar_cell{k}(cellid)];
    end
    
end