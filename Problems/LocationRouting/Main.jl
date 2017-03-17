using Plots

# This calls the functions that are going to be used
include("LocationRoutingFWPH.jl")

# Creating an instance. createinstance(TotalI, TotalJ, TotalV, TotalS, Budget) returns a instance object.
ins = createinstance(20,5,3,10,350.0)

# Calls FWPH. fwph(Instance, Penalty)
tolplot, boundplot, xbar, soltime = fwph(ins, 50.0)

#If you wanna solve the full-space problem:
xref, obj, soltime = solvedeteqmodel(ins)

#if you wanna plot them, use
#plot(boundplot)
#plot(tolplot)

#To get the last objective value, use:
#println(boundplot[end])
