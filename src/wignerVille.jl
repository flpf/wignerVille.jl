include("hilbert.jl")

function wignerVille{T<:Number}(σ::AbstractVector{T},samplingFrequency::Int,tsWinL::Int,fsWinL::Int)
Ψ = hilbert(σ);	
σn=lenght(σ);
σp=div(σn,tsWinL);

end

function pseudoWignerVille{T<:Number}(s::AbstractVector{T},samplingFrequency::Int,fsWinL::Int)

	

end

function smoothedPseudoWignerVille{T<:Number}(s::AbstractVector{T},samplingFrequency::Int,tsWinL::Int,fsWinL::Int)
	

end
