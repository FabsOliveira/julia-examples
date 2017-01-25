#=
Code for using Progressive Hedging in a 2-stage stochastic problem
Created by: Fabricio Oliveira on 25/11/2015
=#

using JuMP, Gurobi

#sets
totalI = 20
totalJ = 4
totalV = 2
totalS = 10
totalN = totalI + totalJ

I = 1:totalI #(set I)
J = totalI+1:totalN #(set J)
V = 1:totalV
S = 1:totalS
N = 1:totalN #(set N)

#parameters
B = Array(Float64, totalI, totalS) #latest arrival time
EC = Array(Float64, totalJ) #installation cost
NV = Array(Float64, (totalI, totalS)) #asset value
ST = Array(Float64, (totalN, totalS)) #operation time
TT = Array(Float64,(totalN, totalN, totalV, totalS)) #routing cost
P = Array(Float64, totalS) #scenario probability

M = 500

# Auxiliary function
function includeVertex(x, Vx, k::Int64)
   dim = tuple(size(x)...,k...)
   flatVx = Vx[:]
   flatx = x[:]

   for i in 1:length(flatx)
      push!(flatVx, x[i])
   end

   Vx = reshape(flatVx, dim)
   return Vx
end

# Main function
function FrankWolfeProgressiveHedging(ρ)

   # control elements
   ϵ = 1e-6 # tolerance
   kMax = 1000 # total iterations

   tol_plot = zeros(0) # register the calculated tolerance
   bound_plot = zeros(0) # register the calculated Lag. bound
   #### remove later!
   of_plot = zeros(0)
   ####

   x_ref = zeros(totalJ) # if full-space solution can be calculates, we store it here

   validate = false # if validation is required, turn this on
   verbose = true # if detailed output is required, turn this on

   #Array Declarations
   x_scenarios = zeros(Float64, (totalJ, totalS))

   z_scenarios = zeros(Float64, (totalI, totalJ, totalS))
   y_scenarios = zeros(Float64, (totalI, totalS))
   ry_scenarios = zeros(Float64, (totalN, totalN, totalV, totalS))
   t_scenarios = zeros(Float64, (totalN, totalV, totalS))

   ω = zeros(totalJ, totalS)
   ωHat = zeros(totalJ, totalS)
   xHat_scenarios = zeros(Float64, totalJ, totalS)

   zHat_scenarios = zeros(Float64, (totalI, totalJ, totalS))
   yHat_scenarios = zeros(Float64, (totalI, totalS))
   ryHat_scenarios = zeros(Float64, (totalN, totalN, totalV, totalS))
   tHat_scenarios = zeros(Float64, (totalN, totalV, totalS))

   if validate ==  true
      DetEq = Model(solver = GurobiSolver(OutputFlag = 1))

      @variable(DetEq, x[J, S], Bin) #server location decision
      @variable(DetEq, z[I, J, S], Bin) #service assignment
      @variable(DetEq, y[I, S], Bin) #service performed
      @variable(DetEq, ry[N, N, V, S], Bin) #sequence of services
      @variable(DetEq, t[N, V, S] >= 0)
      @variable(DetEq, xBar[J], Bin)

      @objective(DetEq, Min,
         sum(P[s]*(
         sum(EC[j-totalI]*x[j,s] for j in J) -
         sum(NV[i,s]*y[i,s] for i in I)  +
         sum(TT[n,np,v,s]*ry[n,np,v,s] for n in N, np in N, v in V)) for s in S))

      @constraint(DetEq, const1[i in I, s in S],
         sum(ry[n,i,v,s] for n in N, v in V) == y[i,s])

      @constraint(DetEq, const2[i in I, j in J, s in S],
         z[i,j,s] <= x[j,s])

      @constraint(DetEq, const3[n in N, i in I, v in V, s in S],
         t[n,v,s] + ST[n,s] + TT[n,i,v,s] <= t[i,v,s] +  M*(1 - ry[n,i,v,s]))

      @constraint(DetEq, const4[i in I, j in J, v in V, s in S],
         sum(ry[j,np,v,s] for np in N) + sum(ry[np,i,v,s] for np in N) <= z[i,j,s] + 1)

      @constraint(DetEq, const5[i in I, s in S],
         sum(z[i,j,s] for j in J) == y[i,s])

      @constraint(DetEq, const6[i in I, v in V, s in S],
         t[i,v,s] <= B[i,s])

      @constraint(DetEq, const7[v in V, s in S],
         sum(ry[n,j,v,s] for n in N, j in J) == 1)

      @constraint(DetEq, const8[n in N, v in V, s in S],
         sum(ry[n,np,v,s] for np in N) - sum(ry[np,n,v,s] for np in N) == 0)

      @constraint(DetEq, NAC[j in J, s in S],
         x[j,s] == xBar[j])

      solve(DetEq)

      x_ref = getvalue(xBar)

   end

   if verbose == true
      println("Intitialing elements...")
   end

   #Obtaining initial x for each scenario so we can calculate an initial x_bar
   for s in S
      if verbose == true
         println("Getting inital x for scenario ", s,"/ ", totalS)
      end
      ScenProb = Model(solver = GurobiSolver(OutputFlag = 0))

      @variable(ScenProb, x[J], Bin) #server location decision
      @variable(ScenProb, z[I, J], Bin) #service assignment
      @variable(ScenProb, y[I], Bin) #service performed
      @variable(ScenProb, ry[N, N, V], Bin) #sequence of services
      @variable(ScenProb, t[N, V] >= 0)

      @objective(ScenProb, Min,
         sum(EC[j-totalI]*x[j] for j in J) -
         sum(NV[i,s]*y[i] for i in I)  +
         sum(TT[n,np,v,s]*ry[n,np,v] for n in N, np in N, v in V))

      @constraint(ScenProb, const1[i in I],
         sum(ry[n,i,v] for n in N, v in V) == y[i])

      @constraint(ScenProb, const2[i in I, j in J],
         z[i,j] <= x[j])

      @constraint(ScenProb, const3[n in N, i in I, v in V],
         t[n,v] + ST[n,s] + TT[n,i,v,s] <= t[i,v] +  M*(1 - ry[n,i,v]))

      @constraint(ScenProb, const4[i in I, j in J, v in V],
         sum(ry[j,np,v] for np in N) + sum(ry[np,i,v] for np in N) <= z[i,j] + 1)

      @constraint(ScenProb, const5[i in I],
         sum(z[i,j] for j in J) == y[i])

      @constraint(ScenProb, const6[i in I, v in V],
         t[i,v] <= B[i,s])

      @constraint(ScenProb, const7[v in V],
         sum(ry[n,j,v] for n in N, j in J) == 1)

      @constraint(ScenProb, const8[n in N, v in V],
         sum(ry[n,np,v] for np in N) - sum(ry[np,n,v] for np in N) == 0)

      solve(ScenProb)

      ####Remove later!
      push!(of_plot, getobjectivevalue(ScenProb))
      ####

      x_scenarios[:,s] = getvalue(x[:])
      z_scenarios[:,:,s] = getvalue(z[:])
      y_scenarios[:,s] = getvalue(y[:])
      ry_scenarios[:,:,:,s] = getvalue(ry[:])
      t_scenarios[:,:,s] = getvalue(t[:])
   end

   if verbose
      println("Scenario-wise solution: $x_scenarios")
      println("Scenario-wise obj. function: $of_plot")
   end

   #Initialising sets Vx and Vy
   Vx = x_scenarios
   Vz = z_scenarios
   Vy = y_scenarios
   Vry = ry_scenarios
   Vt = t_scenarios

   #calculating initial x_bar (z)
   x_bar = P'*x_scenarios'

   #calculating initial \omega
   for s in S
      ω[:,s] = ρ.*(x_scenarios[:,s] - x_bar')
   end

   println("Starting algorithm...")

   #Start main loop
   Converged = false
   k = 1

   while (Converged == false) #main loop
      k = k + 1

      #update multipliers for the linearised Lagrangian
      for s in S
         ωHat[:,s] = ω[:,s] + ρ.* (x_scenarios[:,s] - x_bar')
      end

      #(x^,y^) <- Solve Linearised augmented subproblems for each scenario
      LagBound = 0 # Auxiliary to acumulate the bound value
      for s in S
         augLinSP = Model(solver=GurobiSolver(OutputFlag = 0))

         @variable(augLinSP, x[J], Bin) #server location decision
         @variable(augLinSP, z[I, J], Bin) #service assignment
         @variable(augLinSP, y[I], Bin) #service performed
         @variable(augLinSP, ry[N, N, V], Bin) #sequence of services
         @variable(augLinSP, t[N, V] >= 0)

         @objective(augLinSP, Min,
            sum(EC[j-totalI]*x[j] for j in J) -
            sum(NV[i,s]*y[i] for i in I)  +
            sum(TT[n,np,v,s]*ry[n,np,v] for n in N, np in N, v in V) +
            sum(ωHat[j-totalI,s]*x[j] for j in J)) #workaround for j indexing

         @constraint(augLinSP, const1[i in I],
            sum(ry[n,i,v] for n in N, v in V) == y[i])

         @constraint(augLinSP, const2[i in I, j in J],
            z[i,j] <= x[j])

         @constraint(augLinSP, const3[n in N, i in I, v in V],
            t[n,v] + ST[n,s] + TT[n,i,v,s] <= t[i,v] +  M*(1 - ry[n,i,v]))

         @constraint(augLinSP, const4[i in I, j in J, v in V],
            sum(ry[j,np,v] for np in N) + sum(ry[np,i,v] for np in N) <= z[i,j] + 1)

         @constraint(augLinSP, const5[i in I],
            sum(z[i,j] for j in J) == y[i])

         @constraint(augLinSP, const6[i in I, v in V],
            t[i,v] <= B[i,s])

         @constraint(augLinSP, const7[v in V],
            sum(ry[n,j,v] for n in N, j in J) == 1)

         @constraint(augLinSP, const8[n in N, v in V],
            sum(ry[n,np,v] for np in N) - sum(ry[np,n,v] for np in N) == 0)

         solve(augLinSP)

         # Stores the partial value to calculate the bound
         LagBound = LagBound + P[s]*getobjectivevalue(augLinSP)

         xHat_scenarios[:,s] = getvalue(x[:])
         zHat_scenarios[:,:,s] = getvalue(z[:])
         yHat_scenarios[:,s] = getvalue(y[:])
         ryHat_scenarios[:,:,:,s] = getvalue(ry[:])
         tHat_scenarios[:,:,s] = getvalue(t[:])
      end
      push!(bound_plot, LagBound)

      #Add (x^,y^) to the list of vertices and solve
      Vx = includeVertex(xHat_scenarios, Vx, k)
      Vz = includeVertex(zHat_scenarios, Vz, k)
      Vy = includeVertex(yHat_scenarios, Vy, k)
      Vry = includeVertex(ryHat_scenarios, Vry, k)
      Vt = includeVertex(tHat_scenarios, Vt, k)

      #(x_s,y_s) <- solve master (continuous quadratic) problem
      for s in S
         convSP = Model(solver=GurobiSolver(OutputFlag = 0))

         @variable(convSP, 0 <= λ[1:k] <=1) #creating λ
         @variable(convSP, 0 <= x[J] <= 1)
         @variable(convSP, 0 <= z[I, J] <= 1) #service assignment
         @variable(convSP, 0 <= y[I] <= 1) #service performed
         @variable(convSP, 0 <= ry[N, N, V] <= 1) #sequence of services
         @variable(convSP, t[N, V] >= 0)

         @objective(convSP, Min, #work around to deal with j indexing
            sum(EC[j-totalI]*x[j] for j in J) -
            sum(NV[i,s]*y[i] for i in I)  +
            sum(TT[n,np,v,s]*ry[n,np,v] for n in N, np in N, v in V) +
            sum(ω[j-totalI,s]*x[j] for j in J) +
            sum((ρ/2)*(x[j] - x_bar[j-totalI])^2 for j in J))

         @constraint(convSP, definex[j in J],
            x[j] == sum(Vx[j-totalI, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definez[i in I, j in J],
            z[i, j] == sum(Vz[i, j-totalI, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definey[i in I],
            y[i] == sum(Vy[i, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definery[n in N, np in N, v in V],
            ry[n, np, v] == sum(Vry[n, np, v, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definet[i in I, v in V],
            t[i, v] == sum(Vz[i, v, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, ConvexComb,
            sum(λ[l] for l in 1:k) == 1)

         solve(convSP)

         x_scenarios[:,s] = getvalue(x[:])
         z_scenarios[:,:,s] = getvalue(z[:])
         y_scenarios[:,s] = getvalue(y[:])
         ry_scenarios[:,:,:,s] = getvalue(ry[:])
         t_scenarios[:,:,s] = getvalue(t[:])
      end

      #Update x_bar
      oldxBar = x_bar
      x_bar = P'*x_scenarios'

      #update multipliers
      old_ω = ω
      for s in S
         ω[:,s] = old_ω[:,s] + ρ.*(x_scenarios[:,s] - x_bar')
      end

      #Calculating the reference tolerance parameter:
      tol = 0
      for s in S
         tol = tol + P[s]*(norm(x_scenarios[:,s] - oldxBar'))^2
      end

      push!(tol_plot, tol)

      if verbose
            println("Iteration $k. Current tolerance: $tol. xBar = $x_bar")
      end

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

   return tol_plot, bound_plot
end

include("ReadData.jl")
B, EC, M, NV, ST, P, TT = readData(totalN, totalJ, totalV, totalS)

# Call FWPH
a, b = FrankWolfeProgressiveHedging(50)
