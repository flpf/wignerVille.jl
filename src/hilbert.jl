#Hilbert-transform. Kammeyer, K.D.: Digitale Signalverarbeitung. Teubner. 2006. pp.28 ff.
function hilbert{T<:Number}(σ::AbstractVector{T})
  σₙ=length(σ)
  σₕ=div(σₙ,2)
  σᵧ=σₙ%2
  σ=complex(2*σ)
  FFTW.fft!(σ)
  σ[1]=0
  FFTW.ifft!([σ[1:σₕ];zeros(σₕ+σᵧ)])
end
