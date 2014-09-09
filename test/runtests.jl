using WignerVille
using Base.Test

chirp=readdlm("../files/chirp.txt") # 8000 sample long chirp produced in MATLAB with the chirp()

results=readdlm("../files/results_wvt.txt") # The results of the wigner ville transform from MATLAB


@test wvt(chirp) == results
