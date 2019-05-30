function datahat = demean(data,g) 
    [~, datahat,idx] = aggregate(g, data , @(x) bsxfun(@minus, x, mean(x,1)));
    [~, isrt] = sort(cat(1, idx{:}));
    datahat = cat(1, datahat{:});
    datahat = datahat(isrt,:);    
end