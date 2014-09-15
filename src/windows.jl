
function gaussWindow{T<:Int}(τ::T,ξ::T,Δ::T,γ=1::T)
				#This function returns a Gaussian bell 
				#curve Window. 
				#With total length τ, mid-point ξ, 
				#width Δ at ≈ 78 % of the maximum and 
				#gain γ. 
				return γ*exp(-(([1:τ]-floor(ξ))/(Δ)).^2);
end
