
# scaling function to constrain relevant params between 0 and 1
function phi(x)
    return 0.5 + 0.5 * erf(x / sqrt(2))
end

# reduced versions of dual-component model, for nested comparisons 
etliksepNoMF = (x,d) -> etliksep([x[1:2]; -1e99; -1e99; x[3:end]],d)
etliksepNoGlobal = (x,d) -> etliksep([1e99; -1e99; x[1:2]; -1e99; x[3:end]],d)

# dual-component model
function etliksep(params,data)
    lrT = phi(params[1])
    lrV = phi(params[2])
    lr0 = phi(params[3])
    gamma =  phi(params[4])
    lambda = phi(params[5])
    bint = params[6]
    bV = params[7]
    sig2 = exp(params[8])
    
    paramType = typeof(params[1])
    
    hexStates::Array{Int64,1} = data.pairedHexState .+ 1
    rwd::Array{Int64,1} = data.rwd
    da::Array{Float64,1} = data.DA
    atport::Array{Int64,1} = data.port
    
    V = zeros(paramType,126) .+ 0.2
    Tl = zeros(paramType,3,126)
    
    el = zeros(paramType,126)    
    
    lik = 0.
    
    for i = 1:length(hexStates)-1
        s = hexStates[i]
        snext = hexStates[i+1]


        if (atport[i] == -100)
            # Not at port
            # Gaussian log likelihood
            lik += -1/2 * log(2 * pi * sig2) - 1/(2*sig2) * (da[i] - (bint + bV * V[s]))^2

            V[s] += lr0 * (rwd[i] + gamma * V[snext] - V[s])
            el .*= lambda
            el[s] = 1

        else
            # At port
            #el .*= lambda
            #el[s] = 1

            # update paths in
            Tl[atport[i] + 1,:] = (1-lrT) * Tl[atport[i] + 1,:] + lrT * el

            # Truncated TD0 update
            V[s] += lr0 * (rwd[i] - V[s])

            # TD1 / MB update
            V += lrV .* Tl[atport[i] + 1,:] .* (rwd[i] .- V)
            
            # reset eligibility traces
            el = zeros(paramType,126)
        end
 
    end
    
    return -lik
end

# TD lambda
function tdlambdalik(params,data)
    lr = phi(params[1])
    gamma =  phi(params[2])
    lambda = phi(params[3])
    bint = params[4]
    bV = params[5]
    sig2 = exp(params[6])
   
    paramType = typeof(params[1])
   
    hexStates::Array{Int64,1} = data.pairedHexState .+ 1
    rwd::Array{Int64,1} = data.rwd
    da::Array{Float64,1} = data.DA
    atport::Array{Int64,1} = data.port
   
    V = zeros(paramType,126) .+ 0.2
   
    e = zeros(paramType,126)
   
    lik = 0.
   
    for i = 1:length(hexStates)-1
        s = hexStates[i]
        snext = hexStates[i+1]

        if (atport[i] == -100)
            # Not at port
            # Gaussian log likelihood
            lik += -1/2 * log(2 * pi * sig2) - 1/(2*sig2) * (da[i] - (bint + bV * V[s]))^2
            
            e .*= lambda * gamma
            e[s] = 1

            V += lr .* e .* (rwd[i] + gamma * V[snext] - V[s])

        else
            # At port
            #e .*= lambda * gamma
            #e[s] = 1
           
            # Truncated TD update
            V += lr .* e .* (rwd[i] - V[s])

            # reset eligibility traces
            e = zeros(paramType,126)
        end
 
    end
   
    return -lik
end

#simple 3 port state q learning
function qlik3port(params,data)
        # basic rescorla wagner learning over allocentric ports
        beta = params[1]
        lr = 0.5 + 0.5 * erf(params[2] / sqrt(2))
        bias = params[3]
        bdist = params[4]
        #biasRet = params[5]
        
        c = data[:,:port]
        Q = zeros(typeof(beta),3).+0.5

        lik = 0
    
        s = data[:,:currentport]
        r = data[:,:rwd]

        AC = data[:,:lenAC]
        BC = data[:,:lenBC]
        AB = data[:,:lenAB]
    
        sesh = data[1,:session]
        for i = 1:length(c)
                if (i>1)
                    if (data[i,:session] != sesh)
                        Q = zeros(typeof(beta),3).+0.5
                    end
                    dist = [0 AB[i] AC[i]; AB[i] 0 BC[i]; AC[i] BC[i] 0]
                    inds = mod.(s[i] .+ [0,1],3).+1
                    Q2 = beta * Q
                    Q2[inds[2]] += bias # CW/CCW bias
                    Q2 += bdist * dist[s[i],:] # distance sensitivity
                    lik += Q2[c[i]] - logsumexp(Q2[inds])
                end
                Q[c[i]] = (1-lr) * Q[c[i]] + lr * r[i]
                #Q[c[i]] = Q[c[i]] + lr * (r[i]-Q[c[i]])
                sesh = data[i,:session]
        end

    return -lik
end

#hybrid on- and off-path q learning
function qlikh(params,data)
        # hybrid allo and egocentric updates along specific paths
        beta = params[1]
        lrmf = 0.5 + 0.5 * erf(params[2] / sqrt(2))
        lrmb = 0.5 + 0.5 * erf(params[3] / sqrt(2))
        bias = params[4]
        bdist = params[5]
        #biasRet = params[6]

        Q = zeros(typeof(beta),3,2) .+ 0.5 
    
        lik = 0
    
        
        s = data[:,:currentport]
        r = data[:,:rwd]
        c = data[:,:lrchoice]
        d = data[:,:port]
        
        AC = data[:,:lenAC]
        BC = data[:,:lenBC]
        AB = data[:,:lenAB]
        
        # these are ordered in a left-right logic, unlike the ordering in the
        # allocentric version    
        #dist = [AB AC; BC AB; AC BC]
        sesh = data[1,:session]
        for i = 2:length(c)
            if (data[i,:session] != sesh)
                Q = zeros(typeof(beta),3,2) .+ 0.5 
            end
            dist = [AB[i] AC[i]; BC[i] AB[i]; AC[i] BC[i]]
            Q2 = beta * Q[s[i],:]
            Q2[1] += bias # CW/CCW bias
            #if left taken from last, inhibit right from current port
            #Q2[abs(c[i-1]-2)+1] += biasRet # inhibition of return bias.
            Q2 += bdist * dist[s[i],:] # distance sensitivity
            lik += Q2[c[i]] - logsumexp(Q2)
            Q[s[i],c[i]] = (1-lrmf) * Q[s[i],c[i]] + lrmf * r[i]
            alts = setdiff(1:3,[s[i],d[i]])[1]
            altc = 3-c[i]
        
            Q[alts,altc] = (1-lrmb) * Q[alts,altc] + lrmb * r[i]
            sesh = data[i,:session]

        end

    return -lik
end

#nested likelihood functions 
qlikmb = (x,data) -> qlikh([x[1];x[2];x[2];x[3:end]],data)
qlikmf = (x,data) -> qlikh([x[1];x[2];-1e20;x[3:end]],data)

qliknull = (x,data) -> qlikh([0;0;0;x],data)