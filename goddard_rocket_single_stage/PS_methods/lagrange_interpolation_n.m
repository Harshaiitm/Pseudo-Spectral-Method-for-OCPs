function interpolated_value = lagrange_interpolation_n(collocation_points, function_value, z)
    % Number of data points
    n = length(collocation_points);
    
    % Initialize the interpolated value
    interpolated_value = 0;

    for j = 1:n
        % Calculate the Lagrange basis polynomial for each data point
        Lagrange_polynomial = 1;
        for i = 1:n
            if i ~= j
                Lagrange_polynomial = Lagrange_polynomial.*(z - collocation_points(i))./(collocation_points(j)-collocation_points(i));
            end
        end
        
        % Adding current data point to the interpolation
        interpolated_value = interpolated_value + function_value(j) .* Lagrange_polynomial;
    end
end
