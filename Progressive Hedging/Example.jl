#=
Input initial data (hardcoded for now) and call ProgressiveHedging function
Created by: Fabricio Oliveira on 25/11/2015
=#


c = [150 230 260]

P = [0.25 0.50 0.25]

q = [238 210 -170 -150 -36 -10]

A = [1 1 1]

b = [500]

T = [-2.5   0   0;
        0  -3   0;
        0   0 -20;
        0   0   0]

W = [-1  0 1  0 0 0;
      0 -1 0 +1 0 0;
      0  0 0  0 1 1;
      0  0 0  0 1 0]

h = [-200; -240; 0; 6000]

#Replicating the data scenario-wise (must be removed in later versions with automated input)
α = [0.8 1 1.2] #Yield Adjustment Coefficient for each scenario.

detT = T

m2 = length(h)
dimX = length(c)
totalScenarios = length(P)

T = zeros(m2, dimX, totalScenarios)
for s in 1:totalScenarios
  T[:,:,s] = α[s].*detT
end

include("ProgressiveHedging.jl")
include("FrankWolfeProgressiveHedging.jl")

tol_plot = ProgressiveHedging(c,P,q,A,b,T,W,h,1)
tol_plot2 = FrankWolfeProgressiveHedging(c,P,q,A,b,T,W,h,ρ)
