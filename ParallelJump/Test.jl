@everywhere using JuMP, NLopt

function main(max :: Int)

    argin = [i for i=1:max]
    result = Array(Float64,max,2)

    for i =1:max
        result[i,:] = testjump(argin[i])
    end

    presult = pmap( testjump, argin )

    return result,presult
end

function testjump(seed :: Int)

    N = 20

    srand(seed)
    Z=[ones(N) randn(N)]
    b=[2; 3]
    y=Z*b+rand(N)

    m = Model( solver=NLoptSolver(algorithm=:LD_SLSQP))
    @variable(m, x[ i = 1 : size(Z)[2] ]  )
    @NLobjective(m, Min, sum((y[i]- sum(Z[i,j]*x[j] for j=1:size(Z)[2]))^2 for i=1:size(Z)[1]))

    print(m)

    solve(m)
    println("beta= ",getvalue(x))
    return getvalue(x)
end
