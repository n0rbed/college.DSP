####### TWIDDLE #######
W(m, N) = exp(-im*2pi*m/N)

function comp_twiddle(r, N)
    factors = []
    for m = 0:r-1
        push!(factors, W(m, N))
    end
    return factors
end

# for all powers of 2 in range N, comp twiddle needed only once
function get_twiddle(N)
    twiddle = []
    power_N = log2(N)
    for i = 1:power_N
        push!(twiddle, comp_twiddle((2^i)/2, 2^i))
    end
    twiddle
end


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

########## RADIX 2 FFT #############
function radix2(x)
    N = length(x)
    @assert isinteger(log2(N))
    twiddle = get_twiddle(N)
    radix2_recurse(x, N, twiddle)
end

function radix2_recurse(x, N, twiddle)
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

    A = radix2_recurse(x_even, N/2, twiddle)
    B = radix2_recurse(x_odd, N/2, twiddle) 

    X = [A .+ twiddle[Int(log2(N))] .* B
        A .- twiddle[Int(log2(N))] .* B]
    return X
end


function dft2d(x)
end
