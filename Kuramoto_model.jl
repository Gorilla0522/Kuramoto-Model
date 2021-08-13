using Statistics
using PyPlot

function Kuramoto(φ,ω,K,N)
    rx = Statistics.mean(cos.(φ))
    ry =  Statistics.mean(sin.(φ))
    f = ω + K*(ry*cos.(φ) - rx*sin.(φ))
    return f
end

function Cauthy_distribution(N,Δ,μ)
    c = Δ*tan.(π*(rand(Float64,N) .- 0.5)) .+ μ
    return c
end

function Order_parameter(φ)
    rx = Statistics.mean(cos.(φ))
    ry =  Statistics.mean(sin.(φ))
    return sqrt(rx^2 + ry^2)
end

function Runge_Kutta(dt,T,φ0,ω,K,N)
    step_number = T÷dt
    Rs = []
    φ = copy(φ0)
    for i in 1:step_number
        k1 = Kuramoto(φ,ω,K,N)
        k2 = Kuramoto(φ+0.5dt*k1,ω,K,N)
        k3 = Kuramoto(φ+0.5dt*k2,ω,K,N)
        k4 = Kuramoto(φ+dt*k3,ω,K,N)
        φ += (dt/6)*(k1 + 2k2 + 2k3 + k4)
        push!(Rs,Order_parameter(φ))
    end
    return Rs
end




#変数,定数
N = 1000
dt = 0.02
T = 100
φ0 = rand(N)*2π
ω = Cauthy_distribution(N,1,0)
Ks = collect(0:0.1:4)
Rs = []




@time　for k in Ks
    Rtemp = Runge_Kutta(dt,T,φ0,ω,k,N)
    push!(Rs,Statistics.mean(Rtemp))
end

ylim(0,1)
plot(Ks,Rs,".")