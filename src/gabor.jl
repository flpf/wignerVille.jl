include("Windows.jl")

function dgt{T<:Number}(σ::AbstractVector{T},ω::Int,τ::Int)
	σₙ=length(σ);
	σ=conv(σ,fliplr(σ));
	Δ=div(σₙ,τ);
    τₕ=div(τ,2);
	Γ=zeros(ω,σₙ);
    σ=[zeros(τₕ); σ; zeros(τₕ)];
    Θ=gaussWindow(τ,τₕ,τₕ);
	for r=1:σₙ, p=0:τ-1
		 Γ[p+1,r]=Θ[p+1].*σ[r+p];
	end
	abs(FFTW.fft(Γ,1))
end
