function t = hist2sample(y)

t = [];
c = 0;
for k=1:length(y)
    if(y(k)>0)
        t = [t, k*ones(1,y(k))];
    end
end

end

