function conv1(x, h)
    n = length(x)
    m = length(h)
    O = zeros(1, n+m-1)
    for j = 1:m
        for i = 1:n
            O[i+j-1] += x[i]*h[j]
        end
    end
    return O
end

function conv2d(x, h)
    
end
