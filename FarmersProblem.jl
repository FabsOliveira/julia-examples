#(This example refers to the Farmer's problem example, available in
# Birge & Louveaux - Introduction to Stochastic Programming (1998), Chapter 1.
#
# created by: Fabricio Oliveira

using JuMP, Gurobi #I'm using Gurobi, but if not installed, just remove it from here...

#... and use this instead
#MyModel = Model()
MyModel = Model(solver = GurobiSolver())

##Sets
Prods = 1:3  #set of products i in [1-wheat, 2-corn, 3-sugar beets]
NumProds = length(Prods)
Ranges = 1:2 #set of price ranges j in [1-good, 2-bad]
NumRanges = length(Ranges)

TotalScenarios = 3
Scenarios = 1:TotalScenarios #set of scenarios s
π = ones(TotalScenarios).*(1/TotalScenarios) #Scenario probabilities
α = [0.8 1 1.2] #land production yield variation per scenario

##Parameters
L = 500 #acres of land available
Y = [2.5 3 20] #yield (T/acre)
P = [150 230 260] #planting cost ($/acre)
S = [170 150 36;
       0   0 10] #selling price ($/T) in each range
A = [238 210 0] #purchase price
M = [200 240 0] #minimum requirement
N = 6000 #maximum for selling in good price range

##Define variables generation control
@variable(MyModel, x[i in Prods] >= 0) #total area planted
@variable(MyModel, y[i in Prods, s in Scenarios] >= 0) #total purchased
@variable(MyModel, w[i in Prods,j in Ranges, s in Scenarios] >= 0) #total sold in price range

##Set obj function
@objective(MyModel, Min,
              sum{P[i]*x[i], i in Prods}
            + sum{π[s]*(
              sum{A'[i]*y[i,s], i in Prods}
            - sum{S'[i,j]*w[i,j,s], i in Prods, j in Ranges}), s in Scenarios}
              )

##State constraints
@constraint(MyModel, #1 - land availability
               sum{x[i], i in Prods} <= L)
@constraint(MyModel,  c2[i=1:2, s in Scenarios], #2 - minimum requirements
               α'[s]*Y[i]*x[i] + y[i,s] - sum{w[i,j,s],j=1:2} >= M[i])
@constraint(MyModel, c3[s in Scenarios], #3 - maximum beets sold
               sum{w[3,j,s], j=1:2} - α'[s]*Y[3]*x[3] <= 0)
@constraint(MyModel, c4[s in Scenarios], #4 - maximum for price range
               w[3,1,s] <= N)

#Call solve instruction
status = solve(MyModel)

#Show results
println("Opt. Sol.: " , getvalue(x))
println("z*= ", getobjectivevalue(MyModel))
