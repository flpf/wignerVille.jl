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
  Θ=gaussWindow(τ,τₕ,ξ);
  Θₕ=hilbert(Θ)
  for ι=1+τₕ:σn-τₕ-1
    Ψ[:,count]=[Θₕ[:].*conj(Θₕ[:]).*Φ[ι-τₕ:ι+τₕ-1].*conj(Φ[ι+τₕ-1:-1:ι-τₕ]); zeros(Complex,ω-τ)];
    count+=1;
  end
  FFTW.fft(Ψ,1)
end

function smoothedPseudoWignerVille{T<:Number}(σ::AbstractVector{T},ω::Int,τ::Int,ξ::Int,β::Int,Ξ::Int)
  Φ=hilbert(σ)	
  σn=length(σ)
  σp=div(σn,τ)
  Ψ=zeros(Complex,ω,σn)
  τₕ=div(τ,2)
  βₕ=div(β,2)
  count=1;
  #=Θ=gaussWindow(τ,τₕ,ξ);=#
  Θ=gaussWindow(τ,τₕ,ξ);
  Θₕ=hilbert(Θ)
  Ω=gaussWindow(β,βₕ,Ξ);
  Ωₕ=hilbert(Ω)
  for ι=1+β:σn-β-1
		for tt=1:βₕ 
    	Ψₜ=sum(Θₕ[:].*conj(Θₕ[:]).*Φ[ι-τₕ+tt:ι+τₕ-1+tt].*conj(Φ[ι-τₕ-tt:ι+τₕ-1-tt]));
			Ψ[ω+1-tt,ι]=[Ω[Ωₕ-tt+1].*Ψₜ; zeros(Complex,ω-β)];
		end
  end
  FFTW.fft(Ψ,1) 	
end

function gaussWindow{T<:Int}(τ::T,ξ::T,Δ::T,γ=1::T)
#This function returns a Gaussian bell 
#curve Window. 
#With total length τ, mid-point ξ, 
#width Δ at ≈ 78 % of the maximum and 
#gain γ. 
return γ*exp(-(([1:τ]-floor(ξ))/(Δ)).^2);
end
