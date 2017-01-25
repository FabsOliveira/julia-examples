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

function FrankWolfeProgressiveHedging(c,P,q,A,b,T,W,h,ρ)

   #Reference parameters obtained from the input data
   dimX = length(c) #could be size(A,2) or size(T,2)  as well
   dimY = length(q) #same
   m1 = length(b) #total of 1st-stage constraints
   m2 = length(h) #total of 2nd-stage constraints
   totalScenarios = length(P) #total scenarios
   scenarios = 1:totalScenarios #range of scenarios

   #Array Declarations
   x_scenarios = Array(Float64, dimX, totalScenarios) #we assume that x is column vector (i.e. each column is a scenario and each line a component of x)
   y_scenarios = Array(Float64,dimY, totalScenarios)
   ω = Array(Float64, dimX, totalScenarios)
   ωHat = Array(Float64, dimX, totalScenarios)
   xHat_scenarios = Array(Float64, dimX, totalScenarios)
   yHat_scenarios = Array(Float64, dimY, totalScenarios)

   ϵ = 1e-3
   kMax = 1000

   tol_plot = zeros(0)
   x_ref = zeros(dimX)

   validate = true #If validation is required, turn this on
   verbose = false #If detailed output, turn this on

   if validate ==  true
      detEqProblem = Model(solver=GurobiSolver(OutputFlag = 0))

      @variable(detEqProblem, x[j in 1:dimX] >= 0, Int) #creating x
      @variable(detEqProblem, y[j in 1:dimY, s in scenarios] >= 0, Int) #creating y

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
   for s in scenarios
      scenarioProblem = Model(solver = GurobiSolver(OutputFlag = 0))

      @variable(scenarioProblem, x[1:dimX] >= 0, Int) #creating x
      @variable(scenarioProblem, y[1:dimY] >= 0, Int) #creating y

      @constraint(scenarioProblem, FirstStageConsts[i=1:m1], #Ax = b
                   sum{A[i,j]*x[j], j in 1:dimX} <= b[i])
      @constraint(scenarioProblem, SecondStageConsts[i=1:m2], #Tx + Wy = h
                   sum{T[i,j,s]*x[j], j in 1:dimX} + sum{W[i,j]*y[j], j in 1:dimY} <= h[i])

      @objective(scenarioProblem, Min, sum{c[j]*x[j], j in 1:dimX} + sum{q[j]*y[j], j in 1:dimY})

      solve(scenarioProblem)

      x_scenarios[:,s] = getvalue(x)
      y_scenarios[:,s] = getvalue(y)
   end

   #Initialising sets Vx and Vy
   Vx = x_scenarios
   Vy = y_scenarios

   #calculating initial x_bar (z)
   x_bar = P*x_scenarios'

   #calculating initial \omega
   for s in scenarios
      ω[:,s] = ρ.*(x_scenarios[:,s] - x_bar')
   end

   println("Starting algorithm...")

   #Start main loop
   Converged = false
   k = 0

   while (Converged == false) #main loop
      k = k + 1

      #update multipliers for the linearised Lagrangian
      for s in scenarios
         ωHat[:,s] = ω[:,s] + ρ.* (x_scenarios[:,s] - x_bar')
      end

      #(x^,y^) <- Solve Linearised augmented subproblems for each scenario
      for s in scenarios
         augLinSP = Model(solver=GurobiSolver(OutputFlag = 0))

         @variable(augLinSP, x[1:dimX] >= 0, Int) #creating x
         @variable(augLinSP, y[1:dimY] >= 0, Int) #creating y

         @constraint(augLinSP, FirstStageConsts[i=1:m1], #Ax = b
                        sum{A[i,j]*x[j], j in 1:dimX} <= b[i])

         @constraint(augLinSP, SecondStageConsts[i=1:m2], #Tx + Wy = h
                        sum{T[i,j,s]*x[j], j in 1:dimX} + sum{W[i,j]*y[j], j in 1:dimY} <= h[i])

         @objective(augLinSP, Min, sum{c[j]*x[j], j in 1:dimX} +
                         sum{ωHat[j,s]*x[j], j in 1:dimX} +
                         sum{q[j]*y[j], j in 1:dimY})

         solve(augLinSP)

         xHat_scenarios[:,s] = getvalue(x)
         yHat_scenarios[:,s] = getvalue(y)
      end

      #Add (x^,y^) to the list of vertices and solve
      Vx = includeVertex(xHat_scenarios, Vx, k+1, dimX, totalScenarios)
      Vy = includeVertex(yHat_scenarios, Vy, k+1, dimY, totalScenarios)

      #(x_s,y_s) <- solve master (continuous quadratic) problem
      for s in scenarios
         convSP = Model(solver=GurobiSolver(OutputFlag = 0))

         @variable(convSP, 0 <= λ[1:k] <=1) #creating λ
         @variable(convSP, x[1:dimX] >= 0)
         @variable(convSP, y[1:dimY] >= 0)

         @constraint(convSP, definex[i = 1:dimX],
            x[i] == sum{Vx[i, s, l]*λ[l], l in 1:k})
         @constraint(convSP, definey[i = 1:dimY],
            y[i] == sum{Vy[i, s, l]*λ[l], l in 1:k})
         @constraint(convSP, ConvexComb,
            sum{λ[l], l in 1:k} == 1)

         @objective(convSP,
         Min, sum{c[j]*x[j], j in 1:dimX} +
         sum{ω[j,s]*x[j], j in 1:dimX} +
         sum{(ρ/2)*((x[j] - x_bar[j])^2), j in 1:dimX} +
         sum{q[j]*y[j], j in 1:dimY})

         solve(convSP)

         x_scenarios[:,s] = getvalue(x)
         y_scenarios[:,s] = getvalue(y)
      end

      #Update x_bar
      oldxBar = x_bar
      x_bar = P*x_scenarios'

      #update multipliers
      old_ω = ω
      for s in scenarios
         ω[:,s] = old_ω[:,s] + ρ.*(x_scenarios[:,s] - x_bar')
      end

      #Calculating the reference tolerance parameter:
      tol = 0
      for s in scenarios
         tol = tol + P[s]*norm(x_scenarios[:,s] - oldxBar')
      end

      push!(tol_plot, tol)

      if verbose println("Iteration $k. Current tolerance: $tol") end

      #Checking convergence
      if tol <= ϵ || k == kMax
         Converged = true
         println("Algorithm converged in $k iterations with tolerance = $tol")
      end
   end

   println("Final solution: ", x_bar)
   if validate == true
      println("Reference solution: ", x_ref)
   end

   return tol_plot
end

# Auxiliary functions
function includeVertex(x, Vx, k::Int64, dimX::Int64, totalScenarios::Int64)
   flatVx = Vx[:]
   flatx = x[:]

   for i in 1:length(flatx)
      push!(flatVx, x[i])
   end

   Vx = reshape(flatVx, (dimX, totalScenarios, k))
   return Vx
end
