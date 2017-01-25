using JuMP, Gurobi
using LightGraphs, GraphPlot, Colors

# Sets
totalI = 20
totalJ = 4
totalV = 2
totalS = 10
totalN = totalI + totalJ

I = 1:totalI # set of assets
J = totalI+1:totalN # set of server locations
V = 1:totalV # set of vehicles
S = 1:totalS # set of scenarios
N = 1:totalN # set of nodes (I U J)

include("ReadData.jl")
B, EC, M, NV, ST, P, TT = readData(totalN, totalJ, totalV, totalS)

model = Model(solver = GurobiSolver(OutputFlag = 1))

@variable(model, x[J, S], Bin) # server location decision per scenario
@variable(model, z[I, J, S], Bin) # service assignment
@variable(model, y[I, S], Bin) # service performed
@variable(model, ry[N, N, V, S], Bin) # sequence of services
@variable(model, t[N, V, S] >= 0) # time when protection starts
@variable(model, xBar[J], Bin) # server location decision

@objective(model, Min,
   sum(P[s]*(
   sum(EC[j-totalI]*x[j,s] for j in J) -
   sum(NV[i,s]*y[i,s] for i in I)  +
   sum(TT[n,np,v,s]*ry[n,np,v,s] for n in N, np in N, v in V)) for s in S))

@constraint(model, const1[i in I, s in S],
   sum(ry[n,i,v,s] for n in N, v in V) == y[i,s])

@constraint(model, const2[i in I, j in J, s in S],
   z[i,j,s] <= x[j,s])

@constraint(model, const3[n in N, i in I, v in V, s in S],
   t[n,v,s] + ST[n,s] + TT[n,i,v,s] <= t[i,v,s] +  M*(1 - ry[n,i,v,s]))

@constraint(model, const4[i in I, j in J, v in V, s in S],
   sum(ry[j,np,v,s] for np in N) + sum(ry[np,i,v,s] for np in N) <= z[i,j,s] + 1)

@constraint(model, const5[i in I, s in S],
   sum(z[i,j,s] for j in J) == y[i,s])

@constraint(model, const6[i in I, v in V, s in S],
   t[i,v,s] <= B[i,s])

@constraint(model, const7[v in V, s in S],
   sum(ry[n,j,v,s] for n in N, j in J) == 1)

@constraint(model, const8[n in N, v in V, s in S],
   sum(ry[n,np,v,s] for np in N) - sum(ry[np,n,v,s] for np in N) == 0)

@constraint(model, NAC[j in J, s in S],
   x[j,s] == xBar[j])

solve(model)

println("Optimal solution: ", getvalue(xBar))

#Building graphical representation
printgraph = true
if printgraph
   g = Graph(totalN)
   for n in N, np in N, v in V, s in S
      if getvalue(ry[n,np,v,s]) > 0
         add_edge!(g,n,np)
      end
   end

   labels = 1:totalN
   membership = Array(Int64,totalN)
   for n in N
      if n in I
         membership[n] = 1
      else
         membership[n] = 2
      end
   end
   nodecolor = [colorant"lightseagreen", colorant"orange"]
   # membership color
   colors = nodecolor[membership]
   gplothtml(g , nodelabel = labels, nodefillc=colors)
end
