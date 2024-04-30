% pointx =[1 2 3 4 5 6];
% pointy =[2 4 6 8 10 12];
function state = lagrange_interpolation_z(pointx,pointy,z)
% z =linspace(1,10,6);
n =length(z);
state_polynomial = 1;
state =0;
    for i = 1:n
        for j = 1:n 
            if i ~= j
                state_polynomial=state_polynomial.*(z - pointx(j)) ./ (pointx(i) - pointx(j));
            end
        end
            % state_polynomial = vpa(state_polynomial);
            state = sum(pointy(i).*state_polynomial);
    end
end
