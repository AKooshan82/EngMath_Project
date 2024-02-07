function output = avgFilter(input)
    [m, n] = size(input);
    output = zeros(m, n);
    for i=2:m-1
        for j=2:n-1
            neighbors = input(i-1:i+1, j-1:j+1);
            output(i, j) = mean(neighbors,'all');
        end
    end
end