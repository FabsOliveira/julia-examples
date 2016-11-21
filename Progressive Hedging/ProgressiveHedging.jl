#=
Code for using Progressive Hedging in a 2-stage stochastic problem
Created by: Fabricio Oliveira on 25/11/2015
=#

using JuMP, Gurobi

#=
Problem's current standard format:

Min c'x + \sum_{s \in S}P_s q' y_s
s.t.: Ax <= b
      Tx + W y <= h_s
=#

function ProgressiveHedging(c,P,q,A,b,T,W,h,ρ)

   #Reference parameters obtained from the input data
   dimX = length(c) #could be size(A,2) or size(T,2)  as well
   dimY = length(q) #same
   m1 = length(b) #total of 1st-stage constraints
   m2 = length(h) #total of 2nd-stage constraints
   totalScenarios = length(P) #total scenarios
   scenarios = 1:totalScenarios #range of scenarios

   #Input algorithm's parameters (penalty, initial x_bar, initial multipliers)

   #Define penalty parameter
   println("Obtaining initial value for parameters...")

   ϵ = 1e-3
   kMax = 1000

   tol_plot = zeros(0)
   x_ref = zeros(dimX)

   #If validation is required, turn this on
   validate = true

   if validate ==  true
    detEqProblem = Model(solver = solver=GurobiSolver(OutputFlag = 0))

    @variable(detEqProblem, x[j in 1:dimX] >= 0) #creating x
    @variable(detEqProblem, y[j in 1:dimY, s in scenarios] >= 0) #creating y

    @constraint(detEqProblem, FirstStageConsts[i=1:m1], #Ax = b
                   sum{A[i,j]*x[j], j in 1:dimX} <= b[i])
    @constraint(detEqProblem, SecondStageConsts[i=1:m2, s in scenarios], #Tx + Wy = h
                   sum{T[i,j,s]*x[j], j in 1:dimX} + sum{W[i,j]*y[j,s], j in 1:dimY} <= h[i])

    @objective(detEqProblem, Min, sum{c[j]*x[j], j in 1:dimX} + sum{P[s]*sum{q[j]*y[j,s], j in 1:dimY}, s in scenarios})

    solve(detEqProblem);
    println("Reference problem optimal solution...")
    println("x* = ", getvalue(x))
    println("z* = ", getobjectivevalue(detEqProblem))

    x_ref = getvalue(x)
   end

   #Obtaining initial x for each scenario so we can calculate an initial x_bar

   x_scenarios = ones(dimX, totalScenarios) #we assume that x is column vector (i.e. each column is a scenario and each line a component of x)
   ω = zeros(dimX, totalScenarios)

   for s in scenarios
    scenarioProblem = Model(solver = GurobiSolver(OutputFlag = 0))

    @variable(scenarioProblem, x[1:dimX] >= 0) #creating x
    @variable(scenarioProblem, y[1:dimY] >= 0) #creating y

    @constraint(scenarioProblem, FirstStageConsts[i=1:m1], #Ax = b
                   sum{A[i,j]*x[j], j in 1:dimX} <= b[i])
    @constraint(scenarioProblem, SecondStageConsts[i=1:m2], #Tx + Wy = h
                   sum{T[i,j,s]*x[j], j in 1:dimX} + sum{W[i,j]*y[j], j in 1:dimY} <= h[i])

    @objective(scenarioProblem, Min, sum{c[j]*x[j], j in 1:dimX} + sum{q[j]*y[j], j in 1:dimY})

    solve(scenarioProblem)

    x_scenarios[:,s] = getvalue(x)
   end

   #calculating initial x_bar
   x_bar = P*x_scenarios'

   #calculating initial \omega
   for s in scenarios
      ω[:,s] = ρ.*(x_scenarios[:,s] - x_bar')
   end

   println("Starting algorithm...")

   #Start main loop
   Converged = false
   k = 0

   while (Converged == false)
      k = k + 1

      println("Iteration $k:")
      #Solve augmented subproblems for each scenario
      for s in scenarios
         augScenarioProblem = Model(solver=GurobiSolver(OutputFlag = 0))

         @variable(augScenarioProblem, x[1:dimX] >= 0) #creating x
         @variable(augScenarioProblem, y[1:dimY] >= 0) #creating y

         @constraint(augScenarioProblem, FirstStageConsts[i=1:m1], #Ax = b
                        sum{A[i,j]*x[j], j in 1:dimX} <= b[i])

         @constraint(augScenarioProblem, SecondStageConsts[i=1:m2], #Tx + Wy = h
                        sum{T[i,j,s]*x[j], j in 1:dimX} + sum{W[i,j]*y[j], j in 1:dimY} <= h[i])

         @objective(augScenarioProblem, Min, sum{c[j]*x[j], j in 1:dimX} +
                         sum{ω[j,s]*x[j], j in 1:dimX} +
                         sum{(ρ/2)*((x[j] - x_bar[j])^2), j in 1:dimX} +
                         sum{q[j]*y[j], j in 1:dimY})

         solve(augScenarioProblem)

         x_scenarios[:,s] = getvalue(x)
      end

      #Update x_bar
      x_bar = P*x_scenarios'

      #update multipliers
      old_ω = ω
      for s in scenarios
         ω[:,s] = old_ω[:,s] + ρ.*(x_scenarios[:,s] - x_bar')
      end

      #Calculating the reference tolerance parameter:
      tol = 0
      for s in scenarios
         tol = tol + P[s]*norm(x_scenarios[:,s] - x_bar')
      end

      push!(tol_plot, tol)
      println("Tolerance: $tol")
      #Checking convergence

      if tol <= ϵ || k == kMax
         Converged = true
         println("Algorithm converged in $k iterations with tolerance = ", tol)
      end
   end

   println("Final solution: ", x_bar)
   if validate == true
      println("Reference solution: ", x_ref)
   end
end
