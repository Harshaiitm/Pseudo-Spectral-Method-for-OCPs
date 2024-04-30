function lagrange_polynomial = lagrange_interpolation(pointx, pointy)
    n = length(pointx);
    sympref('FloatingPointOutput',true);
    syms z;
    lagrange_polynomial = 0;

    for i = 1:n
        basis_polynomial = 1;
        for j = 1:n 
            if i ~= j
                basis_polynomial = basis_polynomial * (z - pointx(j)) / (pointx(i) - pointx(j));
            end
        end
        lagrange_polynomial = lagrange_polynomial + pointy(i) * basis_polynomial;
    end
    lagrange_polynomial = vpa(lagrange_polynomial);
end