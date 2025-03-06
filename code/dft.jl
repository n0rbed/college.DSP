# x(n)
function dft(x)
    N = length(x)
    X = zeros(Complex{Float64}, 1, N)
    for m = 1:N
        for n = 1:N
            X[m] += x[n]*exp(-im*2*pi*(n-1)*(m-1)/N)
        end
    end
    X
end
