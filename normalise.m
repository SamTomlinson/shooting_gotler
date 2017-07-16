function [vn] = ...
          normalise(v)

maxs=[];
for j=2:length(v(1,:))-1
    if (v(1,j-1)<v(1,j) && v(1,j)>v(1,j+1))
        maxs=[j,maxs];
    end
    if (v(1,j-1)<v(1,j) && v(1,j)>v(1,j+1))
        maxs=[j,maxs];
    end
    if (v(1,j-1)<v(1,j) && v(1,j)>v(1,j+1))
        maxs=[j,maxs];
    end
    
end

% Normalise

vn=(1/v(1,maxs(1)))*v(1,:);