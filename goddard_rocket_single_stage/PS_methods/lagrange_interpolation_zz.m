function lagrange_polynomial = lagrange_interpolation_zz(pointx, pointy,z)
    n = length(pointx);
    % z = linspace(0, 10, length(pointx)); % Generate z values
    lagrange_polynomial = zeros(size(z));

    for i = 1:n
        basis_polynomial = ones(size(z));
        for j = 1:n 
            if i ~= j
                basis_polynomial = basis_polynomial .* (z - pointx(j)) / (pointx(i) - pointx(j));
            end
        end
        lagrange_polynomial = lagrange_polynomial + pointy(i) * basis_polynomial;
    end
end
