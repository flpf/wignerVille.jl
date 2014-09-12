#=module WignerVilleTransform=#

export wignerVille, pseudoWignerVille, smoothedPseusoWignerVille
include("hilbert.jl")
using Base.LinAlg.BLAS
function wignerVille{T<:Number}(σ::AbstractVector{T},ω::Int,τ::Int)
#The Wigner-Ville transform without windowing.
#
  Φ = hilbert(σ)	
  σn=length(σ)
  σp=div(σn,τ)
  τₕ=div(τ,2)
  Ψ=zeros(Complex,ω,σn-τₕ)
  count=1;
	gc_disable()
  for ι=1+τₕ:σn-τₕ-1
    #=Ψ[:,count]=[Φ[ι-τₕ:ι+τₕ].*conj(Φ[ι+τₕ:-1:ι-τₕ]); zeros(Complex,ω-τ-1)];=#
		Ψ[1:τ,count]=Φ[ι-τₕ:ι+τₕ].*conj(Φ[ι+τₕ:-1:ι-τₕ]);
    count+=1;	
  end
	gc_enable()
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
  Θ=gaussWindow(τ,τₕ,ξ);
  Θₕ=hilbert(Θ)
	gc_disable()
  for ι=1+τₕ:σn-τₕ-1
    #=Ψ[:,count]=[Θₕ[:].*conj(Θₕ[:]).*Φ[ι-τₕ:ι+τₕ-1].*conj(Φ[ι+τₕ-1:-1:ι-τₕ]); zeros(Complex,ω-τ)];=#
    Ψ[1:τ,count]=Θₕ[:].*conj(Θₕ[:]).*Φ[ι-τₕ:ι+τₕ-1].*conj(Φ[ι+τₕ-1:-1:ι-τₕ]);
    count+=1;
  end
	gc_enable()
  FFTW.fft(Ψ,1)
end

function smoothedPseudoWignerVille{T<:Number}(σ::AbstractVector{T},ω::Int,τ::Int,ξ::Int,β::Int,Ξ::Int)
#Smoothed-pseudo-Wigner-Ville transform
  σn=length(σ)
  σp=div(σn,τ)
	offset=β+τ;
	σ=[σ[offset:-1:1]; σ ; σ[σn:-1:σn-(offset)]];
  Φ=hilbert(σ)	
	#=Φ=[Φ[offset:-1:1]; Φ ; Φ[σn:-1:σn-(offset)]]; =#
  Ψ=zeros(Complex,div(ω,2),σn)
  τₕ=div(τ,2)
  βₕ=div(β,2)
  count=1;
  Θ=gaussWindow(τ,τₕ,ξ);
  Ω=gaussWindow(β,βₕ,Ξ);
  Ωₕ=hilbert(Ω)
  Ωₕₓ=Ωₕ[:].*conj(Ωₕ[:])
	Ψₜ=zero(Complex)
	gc_disable()
	for ι=1:σn-1
		for tt=1:βₕ-1
      #=[>Ψₜ=sum(Θ.*Φ[ι-τₕ-tt:ι+τₕ-1-tt].*conj(Φ[ι-τₕ+tt:ι+τₕ-1+tt]));<]=#
			@inbounds Ψₜ=BLAS.dotc(τ,Θ.*Φ[offset+ι-τₕ+tt:offset+ι+τₕ-1+tt],1,Θ.*Φ[offset+ι-τₕ-tt:offset+ι+τₕ-1-tt],1);
			@inbounds Ψ[tt,ι]=Ωₕₓ[tt].*Ψₜ ;
		end
	end
	gc_enable()
  FFTW.fft(Ψ,1) 	
end

function gaussWindow{T<:Int}(τ::T,ξ::T,Δ::T,γ=1.0::Float64)
#This function returns a Gaussian bell 
#curve Window. 
#With total length τ, mid-point ξ, 
#width Δ at ≈ 78 % of the maximum and 
#gain γ. 
return γ*exp(-(([1:τ]-floor(ξ))/(Δ)).^2);
end


#=end=#
