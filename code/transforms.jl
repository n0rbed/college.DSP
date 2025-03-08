W(m, N) = exp(-im*2pi*m/N)

function dft(x)
    N = length(x)
    X = zeros(Complex{Float64}, 1, N)
    for m = 1:N
        for n = 1:N
            X[m] += x[n]*W((n-1)*(m-1), N)
        end
    end
    X
end

function get_twiddle(r, N)
    factors = []
    for m = 0:r-1
        push!(factors, W(m, N))
    end
    return factors
end

function radix2(x)
    N = length(x)
    @assert isinteger(log2(N))
    if N == 2
        return [(x[1] + x[2]); (x[1] - x[2])]
    end

    x_even, x_odd = [], [] 

    for i = 1:length(x)
        if iseven(i-1)
            push!(x_even, x[i])
        else
            push!(x_odd, x[i])
        end
    end

    twiddle = get_twiddle(N/2, N)

    X = [radix2(x_even) .+ twiddle .* radix2(x_odd) 
        radix2(x_even) .- twiddle .* radix2(x_odd)]
    return X
end
