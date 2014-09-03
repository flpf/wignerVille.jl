using WignerVille
using Base.Test

chirp=readdlm("../files/chirp.txt") # 8000 sample long chirp produced with MATLAB
results=readdlm("../files/results_wvt.txt")

@test wvt(chirp) == results
