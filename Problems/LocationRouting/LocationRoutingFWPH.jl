#=
Code for using Progressive Hedging in a 2-stage stochastic location routing problem
Created by: Fabricio Oliveira on 27/02/2017
=#

using JuMP, Gurobi

type InstanceData
   I::UnitRange{Int64} #total clients
   J::UnitRange{Int64} #total centres
   N::UnitRange{Int64} #total nodes = clients + servers
   V::UnitRange{Int64} #total vehicles
   S::UnitRange{Int64} #total scenarios
   B::Array{Float64,2} #time limit for protection
   EC::Array{Float64,1} #instalation cost
   NV::Array{Float64,2} #node value
   ST::Array{Float64,2} #service time
   TT::Array{Float64,4} #travel time
   M::Float64 #Big-M
   P::Array{Float64,1} #probabilities
   Q::Float64 #Budget
end

function createinstance(totalI::Int64, totalJ::Int64, totalV::Int64, totalS::Int64, Q::Float64)
   I = 1:totalI
   J = (totalI+1):(totalI + totalJ)
   S = 1:totalS
   V = 1:totalV

   totalN = totalI + totalJ
   N = 1:totalN

   B = readcsv("B.csv" )
   B = B[I, S]

   EC = readcsv("EC.csv")
   EC = EC[1:totalJ]

   M = 500

   NV = readcsv("NV.csv")
   NV = NV[I, S]

   ST = readcsv("ST.csv")
   ST = ST[N, S]
   ST = round(ST, 3)

   P = 1/totalS.*ones(totalS)

   TTaux = readcsv("TT.csv")
   TTaux = TTaux[1:totalN^2, S]
   TTaux = reshape(TTaux[:, :], (totalN,totalN,1,totalS))
   TT = Array(Float64, (totalN, totalN, totalV, totalS))
   for v in V
      TT[:, :, v, :] = 1.*TTaux[:, :, 1, :]
   end
   TT = round(TT, 3)

   return InstanceData(I, J, N, V, S, B, EC, NV, ST, TT, M, P, Q)
end

function createdeteqmodel(Id::InstanceData, flag = 1)
   DetEq = Model(solver = GurobiSolver(OutputFlag = flag))

   totalI = length(Id.I)

   @variable(DetEq, x[Id.J, Id.S], Bin) #server location decision
   @variable(DetEq, z[Id.I, Id.J, Id.S], Bin) #service assignment
   @variable(DetEq, y[Id.I, Id.S], Bin) #service performed
   @variable(DetEq, ry[Id.N, Id.N, Id.V, Id.S], Bin) #sequence of services
   @variable(DetEq, t[Id.N, Id.V, Id.S] >= 0)
   @variable(DetEq, xBar[Id.J], Bin)

   @objective(DetEq, Min,
      sum(Id.P[s]*(
      -sum(Id.NV[i,s]*y[i,s] for i in Id.I)) for s in Id.S))

   @constraint(DetEq, const1[i in Id.I, s in Id.S],
      sum(ry[n,i,v,s] for n in Id.N, v in Id.V) == y[i,s])

   @constraint(DetEq, const2[i in Id.I, j in Id.J, s in Id.S],
      z[i,j,s] <= x[j,s])

   @constraint(DetEq, const3[n in Id.N, i in Id.I, v in Id.V, s in Id.S],
      t[n,v,s] + Id.ST[n,s] + Id.TT[n,i,v,s] <= t[i,v,s] +  Id.M*(1 - ry[n,i,v,s]))

   @constraint(DetEq, const4[i in Id.I, j in Id.J, v in Id.V, s in Id.S],
      sum(ry[j,np,v,s] for np in Id.N) + sum(ry[np,i,v,s] for np in Id.N) <= z[i,j,s] + 1)

   @constraint(DetEq, const5[i in Id.I, s in Id.S],
      sum(z[i,j,s] for j in Id.J) == y[i,s])

   @constraint(DetEq, const6[i in Id.I, v in Id.V, s in Id.S],
      t[i,v,s] <= Id.B[i,s])

   @constraint(DetEq, const7[v in Id.V, s in Id.S],
      sum(ry[n,j,v,s] for n in Id.N, j in Id.J) == 1)

   @constraint(DetEq, const8[n in Id.N, v in Id.V, s in Id.S],
      sum(ry[n,np,v,s] for np in Id.N) - sum(ry[np,n,v,s] for np in Id.N) == 0)

   @constraint(DetEq, const9[s in Id.S],
      sum(Id.EC[j-totalI]*x[j,s] for j in Id.J) +
      sum(Id.TT[n,np,v,s]*ry[n,np,v,s] for n in Id.N, np in Id.N, v in Id.V) <= Id.Q)

   @constraint(DetEq, NAC[j in Id.J, s in Id.S],
      x[j,s] == xBar[j])

   return DetEq
end

function solvedeteqmodel(Id::InstanceData)
   DetEq = createdeteqmodel(Id)

   soltime = @elapsed solve(DetEq)

   x_ref = getvariable(DetEq, :xBar)
   obj = getobjectivevalue(DetEq)

   return x_ref, obj, soltime
end

function createmodelindstruct(Id::InstanceData, m::Model)

   @variable(m, x[Id.J], Bin) #server location decision
   @variable(m, z[Id.I, Id.J], Bin) #service assignment
   @variable(m, y[Id.I], Bin) #service performed
   @variable(m, ry[Id.N, Id.N, Id.V], Bin) #sequence of services
   @variable(m, t[Id.N, Id.V], lowerbound = 0)

   const1 = @constraint(m, [i in Id.I],
      sum(ry[n,i,v] for n in Id.N, v in Id.V) == y[i])

   const2 = @constraint(m, [i in Id.I, j in Id.J],
      z[i,j] <= x[j])

   const4 = @constraint(m, [i in Id.I, j in Id.J, v in Id.V],
      sum(ry[j,np,v] for np in Id.N) + sum(ry[np,i,v] for np in Id.N) <= z[i,j] + 1)

   const5 = @constraint(m, [i in Id.I],
      sum(z[i,j] for j in Id.J) == y[i])

   const7 = @constraint(m, [v in Id.V],
      sum(ry[n,j,v] for n in Id.N, j in Id.J) == 1)

   const8 = @constraint(m, [n in Id.N, v in Id.V],
      sum(ry[n,np,v] for np in Id.N) - sum(ry[np,n,v] for np in Id.N) == 0)

   return m
end

function solvescenariosp(Id::InstanceData,  ω::Array{Float64}, flag = 0)

   ScenProb = Model(solver = GurobiSolver(OutputFlag = flag))

   totalI = length(Id.I)
   totalJ = length(Id.J)
   totalN = totalI + totalJ
   totalV = length(Id.V)
   totalS = length(Id.S)

   x_scenarios = zeros(Float64, (totalJ, totalS))
   z_scenarios = zeros(Float64, (totalI, totalJ, totalS))
   y_scenarios = zeros(Float64, (totalI, totalS))
   ry_scenarios = zeros(Float64, (totalN, totalN, totalV, totalS))
   t_scenarios = zeros(Float64, (totalN, totalV, totalS))
   LagBound = 0

   ScenProb = createmodelindstruct(Id,ScenProb)

   x = getvariable(ScenProb, :x)
   z = getvariable(ScenProb, :z)
   y = getvariable(ScenProb, :y)
   ry = getvariable(ScenProb, :ry)
   t = getvariable(ScenProb, :t)

   const3 = @constraint(ScenProb, [n in Id.N, i in Id.I, v in Id.V],
      t[n,v] + Id.ST[n,1] + Id.TT[n,i,v,1] <= t[i,v] +  Id.M*(1 - ry[n,i,v]))

   const6 = @constraint(ScenProb, [i in Id.I, v in Id.V],
      t[i,v] <= Id.B[i,1])

   const9 = @constraint(ScenProb,
         sum(Id.EC[j-totalI]*x[j] for j in Id.J) +
         sum(Id.TT[n,np,v,1]*ry[n,np,v] for n in Id.N, np in Id.N, v in Id.V) <= Id.Q)

   for s in Id.S
      @objective(ScenProb, Min,
         -sum(Id.NV[i,s]*y[i] for i in Id.I) +
         sum(ω[j-totalI,s]*x[j] for j in Id.J)) #workaround for j indexing

      if s == 1
         solve(ScenProb)
      else
         const9 = @constraint(ScenProb,
               sum(Id.EC[j-totalI]*x[j] for j in Id.J) +
               sum(Id.TT[n,np,v,s]*ry[n,np,v] for n in Id.N, np in Id.N, v in Id.V) <= Id.Q)

         for i in Id.I, v in Id.V
            JuMP.setRHS(const6[i,v], Id.B[i,s])
            for n in Id.N
               JuMP.setRHS(const3[n,i,v], - Id.ST[n,s] - Id.TT[n,i,v,s] + Id.M)
            end
         end
         solve(ScenProb)
      end
      x_scenarios[:,s] = getvalue(x[:])
      z_scenarios[:,:,s] = getvalue(z[:])
      y_scenarios[:,s] = getvalue(y[:])
      ry_scenarios[:,:,:,s] = getvalue(ry[:])
      t_scenarios[:,:,s] = getvalue(t[:])
      LagBound += Id.P[s]*getobjectivevalue(ScenProb)
   end

   return x_scenarios, z_scenarios, y_scenarios, ry_scenarios, t_scenarios, LagBound
end

# Solution printing...
function printroute(ry::JuMP.JuMPArray{JuMP.Variable,4})
   fullry = getvalue(ry)
   N = 1:JuMP.size(fullry,1)
   V = 1:JuMP.size(fullry,3)
   S = 1:JuMP.size(fullry,4)

   for s in S
      println("Scenario $s")
      for v in V
         route = ""
         for n in N, np in N
            if fullry[n,np,v,s] == 1
                  if route == ""
                     route = string("(", n, "," ,np ,")")
                  else
                     route = string(route, ";", "(", n, ",", np,")")
                  end
            end
         end
         println("Vehicle $v: ", route)
      end
   end
end

# Auxiliary function
function includevertex(x::Array{Float64}, Vx::Array{Float64}, k::Int64)
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
function fwph(Id::InstanceData, ρ::Float64, kMax = 100, validate=false, ϵ = 1e-6)

   Converged = false # main loop flag.
   k = 1 #iteration counter inititalisation
   tol_plot = zeros(0) # register the calculated tolerance
   bound_plot = zeros(0) # register the calculated Lag. bound

   # calculate problem dimenions
   totalI = length(Id.I)
   totalJ = length(Id.J)
   totalN = totalI + totalJ
   totalV = length(Id.V)
   totalS = length(Id.S)

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

   #Deterministic equivalent validation
   if validate
      x_ref, obj = solvedeteqmodel(Id)
   end

   # Initialisationstep for ω = 0
   x_scenarios, z_scenarios, y_scenarios, ry_scenarios, t_scenarios, LagBound =
      solvescenariosp(Id,ω)

   push!(bound_plot, LagBound)

   #Initialising sets Vx and Vy
   Vx = x_scenarios
   Vz = z_scenarios
   Vy = y_scenarios
   Vry = ry_scenarios
   Vt = t_scenarios

   #calculating initial x_bar (z)
   x_bar = Id.P'*x_scenarios'

   println("Iteration ", k, ": bound: ", round(LagBound, 6) ,"; x = ", x_bar)

   #calculating initial \omega
   for s in Id.S
      ω[:,s] = ρ.*(x_scenarios[:,s] - x_bar')
   end

   #Start main loop
   soltime = @elapsed while (Converged == false) #main loop
      k = k + 1

      # update multipliers for the linearised Lagrangian
      for s in Id.S
         ωHat[:,s] = ω[:,s] + ρ.* (x_scenarios[:,s] - x_bar')
      end

      # (x^,y^) <- Solve Linearised augmented subproblems for each scenario
      xHat_scenarios, zHat_scenarios, yHat_scenarios, ryHat_scenarios, tHat_scenarios, LagBound =
         solvescenariosp(Id,ωHat)

      push!(bound_plot, LagBound)

      #Add (x^,y^) to the list of vertices and solve
      Vx = includevertex(xHat_scenarios, Vx, k)
      Vz = includevertex(zHat_scenarios, Vz, k)
      Vy = includevertex(yHat_scenarios, Vy, k)
      Vry = includevertex(ryHat_scenarios, Vry, k)
      Vt = includevertex(tHat_scenarios, Vt, k)

      #(x_s,y_s) <- solve master (continuous quadratic) problem
      for s in Id.S
         convSP = Model(solver=GurobiSolver(OutputFlag = 0))

         λ = @variable(convSP, [1:k], lowerbound = 0, upperbound = 1) #creating λ
         x = @variable(convSP, [Id.J], lowerbound = 0, upperbound = 1)
         z = @variable(convSP, [Id.I, Id.J], lowerbound = 0, upperbound = 1)
         y = @variable(convSP, y[Id.I], lowerbound = 0, upperbound = 1)
         ry = @variable(convSP, [Id.N, Id.N, Id.V], lowerbound = 0, upperbound = 1)
         t = @variable(convSP, [Id.N, Id.V], lowerbound = 0)

         @objective(convSP, Min, #work around to deal with j indexing
            - sum(Id.NV[i,s]*y[i] for i in Id.I)  +
            sum(ω[j-totalI,s]*x[j] for j in Id.J) +
            sum((ρ/2)*(x[j] - x_bar[j-totalI])^2 for j in Id.J))

         @constraint(convSP, definex[j in Id.J],
            x[j] == sum(Vx[j-totalI, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definez[i in Id.I, j in Id.J],
            z[i, j] == sum(Vz[i, j-totalI, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definey[i in Id.I],
            y[i] == sum(Vy[i, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definery[n in Id.N, np in Id.N, v in Id.V],
            ry[n, np, v] == sum(Vry[n, np, v, s, l]*λ[l] for l in 1:k))

         @constraint(convSP, definet[i in Id.I, v in Id.V],
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

      # Update x_bar
      oldxBar = x_bar
      x_bar = Id.P'*x_scenarios'

      # update multipliers
      old_ω = ω
      for s in Id.S
         ω[:,s] = old_ω[:,s] + ρ.*(x_scenarios[:,s] - x_bar')
      end

      # Calculating the reference tolerance parameter:
      tol = sqrt(sum(Id.P[s]*(sum((x_scenarios[j-totalI,s] - oldxBar[j-totalI])^2 for j in Id.J)) for s in Id.S))
      push!(tol_plot, tol)

      println("Iteration ", k, ": bound: ", round(LagBound,6), "; tolerance: ", round(tol, 6), "; x = ", x_bar)

      # Checking convergence
      if tol <= ϵ || k > kMax
         Converged = true
         println("Algorithm converged. ", k, " iterations; bound: ", round(LagBound, 6), "; tolerance: ",
         round(tol, 6), ".")
      end
   end

   # Finalising
   println("Final solution: ", x_bar)
   if validate == true
      println("Reference solution: ", getvalue(x_ref))
   end

   return tol_plot, bound_plot, x_bar, soltime
end
