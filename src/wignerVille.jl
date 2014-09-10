include("hilbert.jl")

function wignerVille{T<:Number}(σ::AbstractVector{T},ω::Int,τ::Int)
#The Wigner-Ville transform without windowing.
#
Φ = hilbert(σ)	
σn=length(σ)
σp=div(σn,τ)
τₕ=div(τ,2)
Ψ=zeros(Complex,ω,σn-τₕ)
count=1;
for ι=1+τₕ:σn-τₕ-1
	Ψ[:,count]=[Φ[ι-τₕ:ι+τₕ].*conj(Φ[ι+τₕ:-1:ι-τₕ]); zeros(Complex,ω-τ-1)];
count+=1;	
end
FFTW.fft(Ψ,1)
end

function pseudoWignerVille{T<:Number}(σ::AbstractVector{T},ω::Int,τ::Int,ξ::Int)
#Pseudo-Wigner-Ville transform with Gauss windowing.
Φ=hilbert(σ)	
σn=length(σ)
σp=div(σn,τ)
Ψ=zeros(Complex,ω,σn)
τₕ=div(τ,2)
count=1;
Ω=gaussWindow(τ,τₕ,ξ);
Ωₕ=hilbert(Ω)
for ι=1+τₕ:σn-τₕ-1
	Ψ[:,count]=[Ωₕ[:].*conj(Ωₕ[:]).*Φ[ι-τₕ:ι+τₕ-1].*conj(Φ[ι+τₕ-1:-1:ι-τₕ]); zeros(Complex,ω-τ)];
count+=1;
end
FFTW.fft(Ψ,1)
end

function smoothedPseudoWignerVille{T<:Number}(s::AbstractVector{T},samplingFrequency::Int,tsWinL::Int,fsWinL::Int)
	

end

function gaussWindow{T<:Int}(τ::T,ξ::T,Δ::T,γ=1::T)
				#This function returns a Gaussian bell 
				#curve Window. 
				#With total length τ, mid-point ξ, 
				#width Δ at ≈ 78 % of the maximum and 
				#gain γ. 
				return γ*exp(-(([1:τ]-floor(ξ))/(Δ)).^2);
end
