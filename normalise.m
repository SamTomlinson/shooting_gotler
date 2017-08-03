%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              normalise                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalises the input vector dependent on its maximum value. Outputs
% the normalise vector.

function [vn] = normalise(v)
   
    % initialise 
    maxs=[];
    % loop through finding maxima
    for j=2:length(v(1,:))-1
        if (v(1,j-1)<v(1,j) && v(1,j)>v(1,j+1))
            maxs=[j,maxs];
        end
    end
    % normalise
    vn=(1/v(1,maxs(1)))*v(1,:);
    
end