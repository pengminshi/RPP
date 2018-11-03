
# interface for fitting rpp model
# input: a list of 
#   year; citperyear; paper title(if plot = T)
rppAnalysis <- function(paper, m.list, fit.nyears = NULL, nyear.extend = 0,
                        real.only = F, k.smoothing = 0, decay.type = 'lognormal',
                        verbose = F, max.iter = 1000, plot = F,
                        mu.init = 5, sigma.init = 1,eps = 10^(-8),
                        n.day = 365, main = NULL, awake.time = NULL,
                        alpha = 0, beta = 0){
    
    years = paper$years
    n.years = length(years)
    n.m = length(m.list)
    
    citperyear = paper$citperyear
    
    if(plot){
        legend.y = max(citperyear) * 3/2
        if (is.null(main)){
            paper.title = paper$paper.title
            if (nchar(paper.title)>40){
                paper.title = paste(substr(paper.title,1,40),'-\n',
                                    substr(paper.title,41,nchar(paper.title)),sep = '')
                if(nchar(paper.title) >85){
                    paper.title = paste(substr(paper.title,1,85),'-\n',
                                        substr(paper.title,86,nchar(paper.title)),sep = '')
                }
            }
            main = paper.title
        }
        
        plot(years, citperyear, type = 'o', pch = 1, col = 1, 
             xlim = c(min(years)-1, max(years)+ nyear.extend),
             ylim = c(0,max(citperyear)*3/2),
             main = main)
    }
    
    if (real.only)
        return(0)
    
    if ( k.smoothing > 0 ) {
        citperyear = citationSmoothing(citperyear, k.smoothing)
    }
    
    if (is.null(fit.nyears)){
        fit.nyears = n.years
    } else if( fit.nyears > n.years){
        cat('Warning: not enough years, Fit ', n.years, ' years.\n')
        fit.nyears = n.years
    }
    
    if(!is.null(awake.time)){
        year.start.n= which(years==awake.time)
        if (year.start.n > fit.nyears -5){
            cat('Warning! Less than 5 years to fit given the sleeping time.')
            return(-1)
        }
    } else {
        year.start.n = 1
    }
    
    citation.times = rep((1:(fit.nyears-year.start.n+1)-0.5) * n.day, 
                         citperyear[year.start.n:fit.nyears])
    
    cat('Fitting ...\n')
    fit.out.list = lapply(m.list, function(m){
        cat('m = ', m,'\t')
        fit.out = fitRPP(citation.times, m = m, verbose = verbose, eps = eps,
                         max.iter = max.iter, mu.init = mu.init, 
                         sigma.init = sigma.init, decay.type = decay.type,
                         alpha = alpha, beta = beta)
        cat('Converge?', fit.out$converge,', lambda=',fit.out$lambda, 
            ', mu=', fit.out$mu, ', sigma=',fit.out$sigma, '\n')
        
        return(fit.out)
    })
    
    if(plot){
        years.extend = years[year.start.n]:(max(years) + nyear.extend)
        for ( i in 1:n.m){
            new.cit = citationFit( to.date = (n.years+nyear.extend-year.start.n+1) *365, 
                                   lambda = fit.out.list[[i]]$lambda, 
                                   mu = fit.out.list[[i]]$mu, 
                                   sigma = fit.out.list[[i]]$sigma, 
                                   m = m.list[i],
                                   decay.type = decay.type)
            lines(years.extend, new.cit$yearcounts, type = 'o', col = i+1, pch = i+1)
        }
        abline(v = years[fit.nyears]+0.5, lty = 'dashed')
        legend.text = c('real data',
                        sapply(m.list, function(m) paste('fit m=', as.character(round(m,1)))))
        legend(min(years),legend.y, 
               legend = legend.text, 
               col = 1:(n.m+1), pch = 1:(n.m+1))
    }
    
    out = list(m.list = m.list, fit.parameters = fit.out.list)
    
    return(out)
}

#' fit RPP model
#' 
#' fit Reinforced poisson process model(RPP) proposed by Dashun Wnang (2013), 
#' by maximize log likelihood function
#' 
#' @param citation.times a vector indicating the arriving time of each citation {t_i}_{i=1}^n
#'                       (by number of days after the publish date).
#' @param m (optional) the global constant, suggested m = 30 by their science paper, as the defalt.
#' @param time.T (optional) the observation time [0,T], so that 0<= t1 <= t2<= ... <= tn <= T. 
#'               If not given, use the last citation arriving time as T, so that T = tn.
#' @param verbose (optional) whether to output esitmation at each step, by default is False.
#' @param mu.init (optional) the initial value of mu.
#' @param sigma.init (optional) the initial value of sigma
#' @param max.iter (optional) the maximum number of iteration in gradient descent, by default 1000.
#' @param eps (optional) topping accuracy, by default 10^(-8).
#' @param decay.type (optional) the type of decay function f_d(t) that characterize the age of citation intensity,
#'                    by default is lo gnormal suggested by science paper. Other types are to be implemented.
#' @param alpha/beta (optional) the prior of lambda, used for RPP model with prior, proposed in (Huawei Shen, etc 2014),
#'                   by default there is no prior (alpha = 0, beta = 0). 
#' 
#' @return A list containing  \describe{
#'   \item{mu}{estimation of mu.}
#'   \item{sigma}{estimation of sigma.}
#'   \item{lambda}{estimation of lambda.}
#'   \item{converge}{boolean, whether the optimization algorithm coverges within \code{max.iter} iterations.}
#' }
fitRPP <- function(citation.times,  m = 30, time.T = NULL,  verbose = FALSE,
                   mu.init = 9, sigma.init = 1, max.iter = 1000, eps = 10^(-8),
                    decay.type = 'lognormal', alpha = 0, beta = 0){
    if ( length(citation.times) == 0){
        cat('Empty citation!')
        return(0)
    } 
    nd = length(citation.times)
    
    if (is.null(time.T)){
        time.T = citation.times[nd]
    }
    
    # initialize
    mu = mu.init
    sigma = sigma.init
    
    ll.out = logLikelihood(citation.times, mu = mu, sigma = sigma, time.T = time.T, 
                           m = m, decay.type = decay.type, alpha = alpha, beta = beta)
    ll = ll.out$ll
    lambda = ll.out$lambda

    for ( i in 1:max.iter){
        
        grad.out = computeGradient(citation.times, lambda, mu, sigma, time.T, m, 
                                   decay.type = decay.type)
        
        ls.out = lineSearch(citation.times, mu = mu, sigma = sigma, ll = ll, 
                            grad.out$mu, grad.out$sigma, time.T, m, 
                            decay.type = decay.type, alpha = alpha, beta = beta)
        
        g.square = grad.out$mu^2 + grad.out$sigma^2
        
        if (verbose == T){
            cat('iter:',i, '\nLL:',ll, ' lambda:', lambda,' mu:', mu, ' sigma:', sigma,'\n')
            cat(' g_square:', g.square,'\n--------------------------------\n')
        }
        
        if ( g.square < eps) break
        
        ll = ls.out$ll
        mu = ls.out$mu
        sigma = ls.out$sigma
        lambda = ls.out$lambda
        
    }
    
    if ( i == max.iter){
        if (verbose == T){
            cat('WARNING! Reached maximum ', i, 'iterations!\n')
        }
        converge = F
    } else{
        if (verbose == T){
            cat('Converge!', i, ' iterations.','\n')
        }
        converge = T
    }
    
    return(list(mu = mu, sigma = sigma, lambda = lambda, converge = converge))
}


# compute parameter lambda
computeLambda <- function(citation.times, mu, sigma, time.T = NULL, m, 
                          defalt.lambda = 1, decay.type, alpha = 0, beta =0){
    
    nd = length(citation.times)
    
    if (is.null(time.T)){
        time.T = citation.times[nd]
    }
    
    tmp  = (nd + m) * Fd(time.T, mu, sigma, decay.type)
    for (i in 1:nd){
        Fd.ti = Fd(citation.times[i], mu, sigma, decay.type)
        tmp = tmp - Fd.ti
    }
    tmp = max(tmp, 0) # make sure tmp is positive
    
    if(tmp+beta == 0){
        lambda = defalt.lambda
    } else{
        lambda =  (alpha+nd) / (beta +tmp)
    } 
    
    return(lambda)
}

# calculate the log likelihood function
# alpha, beta are the prior of lambda defined in the paper (Huawei Shen, etc 2014), 
# which is an extension of original model.
logLikelihood <- function(citation.times, mu, sigma, time.T = NULL, m, 
                          lambda = NULL, decay.type, alpha = 0, beta = 0){
    nd = length(citation.times)
    
    
    if (is.null(lambda)){
        lambda = computeLambda(citation.times, mu = mu, sigma = sigma, time.T = time.T,
                               m = m, decay.type = decay.type,
                               alpha = alpha,beta = beta)
    }
    
    if (beta == 0){
        ll = 0
    } else {
        ll = 0 # alpha * log(beta) + log(gamma(nd+alpha)/gamma(alpha)) ??need fix
    }
    
    ll = ll + (nd + alpha)* log(lambda/(alpha + nd)) 
    ll = ll + sum(sapply(1:nd, function(i){
       
        log((m+i-1) * fd(citation.times[i], mu = mu, sigma = sigma, decay.type = decay.type))
        
    }))
    
    return(list(ll = - ll, lambda = lambda ))
}

# function used in computeGradient
gammaType <- function(t, mu, sigma, decay.type){
    
    if (decay.type == 'lognormal'){
        return((log(t) - mu) / sigma)
    } else {
        cat('Error! Wrong decaying function!')
    }
}

# compute the gradient of mu, sigma
computeGradient <- function(citation.times, lambda, mu, sigma, 
                            time.T = NULL, m, decay.type){
    nd = length(citation.times)
    if (is.null(time.T)){
        time.T = citation.times[nd]
    }
    
    gamma.d = gammaType(time.T, mu, sigma, decay.type)
    
    tmp = lambda * (nd + m) * dnorm(gamma.d)
    
    gradient.mu = tmp
    gradient.sigma = tmp * gamma.d - nd
    
    if (nd > 0){
        for ( i in 1:nd){
            gamma.i = gammaType(citation.times[i], mu, sigma, decay.type)
            tmp = gamma.i -lambda * dnorm(gamma.i)
            gradient.mu = gradient.mu + tmp
            gradient.sigma = gradient.sigma + gamma.i * tmp
        }
    }
    
    gradient.mu = gradient.mu / sigma
    gradient.sigma = gradient.sigma / sigma
    
    return(list(mu = -gradient.mu, sigma = -gradient.sigma))
}

# find the stepsize by line search
lineSearch <- function(citation.times, mu, sigma, ll, gradient.mu, gradient.sigma, 
                       time.T = NULL, m, mu.min = 0.5, sigma.min = 0.1, eta = 0.5,
                       decay.type, lambda = NULL, alpha = 0, beta = 0){
    t = 100
    nd = length(citation.times)
    
    if(is.null(time.T) ){
        time.T = citation.times[nd]
    }
    
    if(is.null(ll)){
        ll = logLikelihood(citation.times, mu = mu, sigma = sigma, 
                           time.T = time.T, m = m, alpha = alpha, beta = beta,
                           decay.type = decay.type, lambda = lambda )$ll
    }
    
    # check for maximum t it satisfy the boundary contraint, otherwise shrink it
    if ((sigma - t * gradient.sigma) <  sigma.min){
        if (sigma <  sigma.min){
            cat('sigma given for line search is outside boundary!\n')
            # return(list(mu = mu, sigma = sigma, ll = ll, lambda = lambda))
        } else {
            t = (sigma - sigma.min) / gradient.sigma
        }
    }
    
    if ((mu - t* gradient.mu) < mu.min){
        if (mu < mu.min){
            cat("mu given for line search is outside boundary!\n")
        } else {
            t = (mu - mu.min) / gradient.mu
        }
    }
    
    mu.new = mu - t * gradient.mu
    sigma.new = sigma - t * gradient.sigma
    
    # cat('linesearch, mu.new=',mu.new,'sigma=',sigma.new,'\n')
    ll.out = logLikelihood(citation.times = citation.times, mu = mu.new, sigma = sigma.new, 
                           time.T = time.T, m = m, decay.type = decay.type,
                           alpha = alpha, beta = beta)
    tmp = 1/2 * (gradient.mu^2 + gradient.sigma^2)
    ll.new = ll.out$ll

    
    if (is.na(ll.new)){
        cat('WARNING! Loglikelihood is NA.\n')
        cat(ll.new, ll-t*tmp,ll.out$lambda, mu.new, sigma.new, time.T, m,'\n')
    }
    
    while( ll.new > ll - t * tmp){
        t = t * eta
        
        mu.new = mu - t * gradient.mu
        sigma.new = sigma - t * gradient.sigma
        
        ll.out = logLikelihood(citation.times, mu.new, sigma.new, time.T, 
                               m = m, decay.type = decay.type,
                               alpha = alpha, beta = beta)
        ll.new = ll.out$ll
    }
    
    out = list(mu = mu.new, sigma = sigma.new, ll = ll.new, lambda = ll.out$lambda)
    return(out)
}

# decaying function f_d(t; mu, sigma)
Fd <- function(t, mu, sigma, decay.type){
    if (decay.type == 'lognormal'){
        return(plnorm(t, mu, sigma))
    } else {
        cat('Error! Wrong decaying function type!\n')
        return(0)
    }
}

# accumulative decay function integral of f_d(s;mu,sigma)
fd <- function(t, mu, sigma, decay.type){
    if (decay.type == 'lognormal'){
        return(dlnorm(t, mu, sigma))
    } else {
        cat('Error! Wrong decaying function type!\n')
        return(0)
    }
}


# generate the citation arriving times between [0, time.T], given the model parameters
citationGenerator <- function( time.T, lambda, mu, sigma, m, decay.type = 'lognormal'){
    
    # y = accmulativeCitation(x, lambda, mu, sigma, m, decay.type = decay.type)
    accum.T = accmulativeCitation(time.T, lambda = lambda, mu = mu, sigma = sigma, 
                                  m = m, decay.type = decay.type)
    
    n.citation = ceiling(accum.T)
    citation.times = arriveTime(n.citation, lambda = lambda, mu = mu, sigma = sigma, 
                                m = m, decay.type = decay.type)
    
    citation.times = citation.times[citation.times <= time.T] # cut the Inf time
    
    return(citation.times)
}

# calculate accumlative citations given the model parameters
accmulativeCitation <- function(times, lambda, mu, sigma, m, decay.type){
    accum = m * (exp(lambda * Fd(times, mu,sigma, decay.type)) - 1)
    return(accum)
}


# calculate the inverse of the accumulative function
inverseFunc <- function(y, lambda, mu, sigma, m, decay.type){
    quantile = 1/lambda * log(y/m + 1)
    if (quantile > 1){
        cat('inverseFunc', y,'\n')
    }
    t = fqd(quantile, mu, sigma, decay.type) 
    
    return(t)
}


# calculate arriving times of papers, ti, i = 1,...,n
arriveTime <- function(n, lambda, mu, sigma, m, decay.type){
    ys = 1:n
    times = sapply(ys, function(y) inverseFunc(y, lambda, mu, sigma, m, 
                                               decay.type = decay.type))
    
    return(times)
}

# quantile function of fd
fqd <- function(q, mu, sigma, decay.type){
    if (decay.type == 'lognormal'){
        return(qlnorm(q, mu, sigma))
    } else if(decay.type == 'normal'){
        return(qnorm(q, mu, sigma))
    } else {
        cat('Error! Wrong decaying function decay.type!\n')
        return(0)
    }
}

# count the number of citations by year
citationYearlyCount <- function(citation.times,  n.day = 365){
    
    year.max = ceiling(max(citation.times)/n.day)
    
    yearly.counts = sapply(1: year.max, function(y){
        length(which(citation.times > (y-1)*n.day & citation.times <= y*n.day))
    } )

    return(yearly.counts)
}

# given the parameters of the citation model
# output the acculative citation curve and the citation per years 
citationFit <- function( to.date, lambda, mu, sigma, m,
                         delta = 1, n.day = 365, decay.type = 'lognormal'){
    
    x = seq(0, to.date-delta, delta)
    y = accmulativeCitation(x, lambda, mu, sigma, m, decay.type = decay.type)
    
    if (max(y) > m*(exp(lambda)-1)){
        cat('WARNING! Accumulative citations exceeds the maximum count!\n')
        y[y > m*(exp(lambda)-1)] = m*(exp(lambda)-1)
    }
    
    ay.out = accumToYeartime(x, y, n.day = n.day)
    
    out = list(acc.x = x,
               acc.y = y,
               arr.yeartimes = ay.out$yeartimes,
               yearcounts = ay.out$yearcounts)
    
    return(out)
    
}

# function used by citationFit
accumToYeartime<- function(times.x, accum.y, n.day = 365){
    year.max = ceiling(max(times.x)/n.day)
    
    year.counts = rep(0, year.max)
    acc.sum = accum.y[1]
    for( i in 1:year.max){
        ind.i1 = max(which(times.x <= i * n.day))
        year.counts[i] = round(accum.y[ind.i1]) - acc.sum
        acc.sum = acc.sum + year.counts[i]
    }
    year.times = rep(1:year.max, year.counts)
    out = list(yeartimes = year.times, yearcounts = year.counts)
    
    return(out)
}

# smoothing the citation with the k nearest neigher years
citationSmoothing <- function(citperyear, k.smoothing){
    n.year = length(citperyear)
    
    if (n.year < k.smoothing + 1){
        cat('WARNING! Not enough years to smooth.')
        k.smoothing = n.year - 1
    }
    
    citperyear.new = citperyear
    for ( i in (k.smoothing+1): n.year){
        smoothing.ind = max(1, i-k.smoothing) : min(n.year, i + k.smoothing)
        citperyear.new[i] = mean(citperyear[smoothing.ind])
    }
    citperyear.new = ceiling(citperyear.new)
    
    return(citperyear.new)
}
