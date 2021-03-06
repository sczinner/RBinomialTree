#Generates a binomial lattice for stock price, returns its parameters
#assuming stock can go up to S_{t-1}u_ or down to S_{t-1}d_
#uses u_=e^{sigma*sqrt{h}}, d_=e^{-sigma*sqrt{h}},
#or  u_=e^{rh}+delta, d_=e^{rh}-delta, delta=e^{rh}*sqrt(e^(sigma^2h)-1)
#T_ is the number of years
#h_ is the period time in years
#is_rc is TRUE if recombining, otherwise false
#S_0 is the initial stock price
#r_ is the risk free rate (generally assumed>0)
#sigma is the volatility (>0)
#Also requires T_%%h_=0 (i.e number of levels should be an integer)
#method 1 is Cox-Ross-Rubinstein Tree, method 2 is the Jarrow-Rudd Risk Neutral Tree

library(zoo)#used for computing rolling average (asian options)
binomial_lattice<-function(T_,h_, is_rc, S_0, r_, sigma,method){
  model_parameters<-list()
  
  n_levels<-T_/h_#number of levels in the lattice
  if(method==1){
    u_<-exp(sigma*sqrt(h_))
    d_<-1/u_
  }else{
    delta<-exp(r_*h_)*sqrt(exp(sigma^2*h_)-1)
    u_<-exp(r_*h_)+delta
    d_<-exp(r_*h_)-delta
  }
  lattice<-vector("list",length = n_levels+1)
  
  for(ii in 1:(n_levels+1)){
    level<-ii-1# note i'm doing this because of R not being zero-indexed
    if(ii==1){
      row=list(list(price = S_0))# set initial stock price to S_0
      model_parameters<-list(T_=T_,h_=h_,is_rc=is_rc,r_=r_,sigma=sigma, method=method,u_=u_,d_=d_)
      #also store parameters here 
      #a row is a list of nodes, each node has a list of characteristics like price
    }else{
      if(is_rc){
        l1<-length(lattice[[ii-1]])#number of nodes in in previous period
        row<-vector("list",length = l1+1)#one additional node per period
        for(jj in 1:length(lattice[[ii-1]])){
          row[[jj]]$price<-d_*lattice[[ii-1]][[jj]]$price
        }
        row[[l1+1]]$price<-u_*lattice[[ii-1]][[l1]]$price#set the additional node
      }else{
        row<-vector("list",length = 2^(level))
        for(jj in 1:length(lattice[[ii-1]])){
          row[[jj*2-1]]$price<-d_*lattice[[ii-1]][[jj]]$price
          row[[jj*2]]$price<-u_*lattice[[ii-1]][[jj]]$price
        }
      }
    }
    
    lattice[[ii]]<-row
  }
  return(list(lattice,model_parameters))
}



call_payoff<-function(lattice,param){
  K<-param[[1]]
  final_period=length(lattice)
  for(ii in 1:length(lattice[[final_period]])){
    lattice[[final_period]][[ii]]$payoff<-max(lattice[[final_period]][[ii]]$price-K,0)
  }
  #the idea is to return a lattice but with an additional payoff parameter in the final layer.
  return(lattice)
}

put_payoff<-function(lattice,param){
  K<-param[[1]]
  final_period=length(lattice)
  for(ii in 1:length(lattice[[final_period]])){
    lattice[[final_period]][[ii]]$payoff<-max(K-lattice[[final_period]][[ii]]$price,0)
  }
  #the idea is to return a lattice but with an additional payoff parameter in the final layer.
  return(lattice)
}

#a recursive function which returns the index in the row of a levelth degree ancestor 
#from a node given by its index in its row
getNodeParent<-function(index,levels){
  if(levels==0){
    return(index)
  }else{
    return(getNodeParent(floor((index+1)/2),levels-1))
  }
}

#asian options are path dependent and require non-recombining lattice
asian_call_payoff<-function(lattice,param){
  K<-param[[1]]
  numtoaverage<-min(param[[2]],length(lattice)-1)
  
  final_period=length(lattice)
  for(ii in 1:length(lattice[[final_period]])){
    average_price<-lattice[[final_period]][[ii]]$price/numtoaverage
    if(numtoaverage-1>0){
      for(jj in 1:(numtoaverage-1)){
        average_price<-average_price+lattice[[final_period-jj]][[getNodeParent(ii,jj)]]$price/numtoaverage
      }
    }
    
    
    lattice[[final_period]][[ii]]$payoff<-max(average_price-K,0)
  }
  return(lattice)
}

#asian options are path dependent and require non-recombining lattice
asian_put_payoff<-function(lattice,param){
  K<-param[[1]]
  numtoaverage<-min(param[[2]],length(lattice)-1)
  
  final_period=length(lattice)
  for(ii in 1:length(lattice[[final_period]])){
    average_price<-lattice[[final_period]][[ii]]$price/numtoaverage
    if(numtoaverage-1>0){
      for(jj in 1:(numtoaverage-1)){
        average_price<-average_price+lattice[[final_period-jj]][[getNodeParent(ii,jj)]]$price/numtoaverage
      }
    }
     
    lattice[[final_period]][[ii]]$payoff<-max(K-average_price,0)
  }
  return(lattice)
}


price_binomial_lattice<-function(lattice, is_european,payoff_function,param=NULL){
  model_parameters<-lattice[[2]]
  lattice<-lattice[[1]]
  #if european then call the payoff on the final time period in the lattice, otherwise every time period
  #if european then use the fact that discounted price is a martingale under risk neutral measure
  #if american use the fact above but also consider the option of exercising at all time periods
  #the idea is to return a lattice but with an additional option price parameter in the final layer. 
  is_rc<-model_parameters$is_rc
  r_<-model_parameters$r_
  h_<-model_parameters$h_
  u_<-model_parameters$u_
  d_<-model_parameters$d_
  qu<-(exp(r_*h_)-d_)/(u_-d_)#risk neutral measure in binomial model
  qd<-1-qu
  
  if(is_european){
    lattice<-payoff_function(lattice,param)
    for(ii in length(lattice):1){#backwards for recursive pricing
      for(jj in 1:length(lattice[[ii]])){
        if(ii==length(lattice)){
          lattice[[ii]][[jj]]$optionprice<-lattice[[ii]][[jj]]$payoff
        }else{
          if(is_rc){
            lattice[[ii]][[jj]]$optionprice<-exp(-r_*h_)*(lattice[[ii+1]][[jj]]$optionprice*qd+lattice[[ii+1]][[jj+1]]$optionprice*qu)
          }else{
            lattice[[ii]][[jj]]$optionprice<-exp(-r_*h_)*(lattice[[ii+1]][[2*jj-1]]$optionprice*qd+lattice[[ii+1]][[2*jj]]$optionprice*qu)
          }
        }
      }
    }
  }else{
    for(ii in  1:length(lattice)){
      lattice[[ii]]<-payoff_function(lattice[1:ii],param)[[ii]]#computes payoffs from exercising at all times
    }
    for(ii in length(lattice):1){#backwards for recursive pricing
      for(jj in 1:length(lattice[[ii]])){
        if(ii==length(lattice)){
          lattice[[ii]][[jj]]$optionprice<-lattice[[ii]][[jj]]$payoff
        }else{
          if(is_rc){
            continuingprice<-exp(-r_*h_)*(lattice[[ii+1]][[jj]]$optionprice*qd+lattice[[ii+1]][[jj+1]]$optionprice*qu)
            lattice[[ii]][[jj]]$optionprice<-max(continuingprice,lattice[[ii]][[jj]]$payoff)
          }else{
            continuingprice<-exp(-r_*h_)*(lattice[[ii+1]][[2*jj-1]]$optionprice*qd+lattice[[ii+1]][[2*jj]]$optionprice*qu)
            lattice[[ii]][[jj]]$optionprice<-max(continuingprice,lattice[[ii]][[jj]]$payoff)
          }
        }
      }
    }
  }
  return(list(lattice,model_parameters))
}

print_binomial_lattice<-function(lattice,parameter){
  lines<-list()
  longest<-0
  for(ii in length(lattice):1){#go through the list backwards so get the longest line first
      str<-""
      
    for(jj in 1:length(lattice[[ii]])){
      num<-format(lattice[[ii]][[jj]][parameter],digits=2,nsmall=2)
      padding<-max(longest/(length(lattice[[ii]])+1)-nchar(num),0)
      str<-paste(str,paste(strrep(" ",padding)))#padding to center
      str<-paste(str,num)
    }
    
    if(ii==length(lattice)){
      longest<-nchar(str)
    }
    lines[[ii]]<-str
  }
  
  for(ii in 1:length(lines)){
    cat(lines[[ii]])
    cat("\n")
  }
}




#monte-carlo implementations - does not compute the entire lattice but rather some simulated paths
#this is useful for path-dependent options due to otherwise exponential lattice size, for n periods
#a non-recombining tree has O(2^n) nodes and a recombining tree has O(n^2) nodes. It is still useful
#for recombining trees and path-independent options but less so. 

mc_paths<-function(T_,h_, S_0, r_, sigma,method,N_){
  model_parameters<-list()
  
  n_levels<-T_/h_#number of levels in the lattice
  if(method==1){
    u_<-exp(sigma*sqrt(h_))
    d_<-1/u_
  }else{
    delta<-exp(r_*h_)*sqrt(exp(sigma^2*h_)-1)
    u_<-exp(r_*h_)+delta
    d_<-exp(r_*h_)-delta
  }
  qu<-(exp(r_*h_)-d_)/(u_-d_)
  qd<-1-qu
    
  paths<-list()
  paths$prices<-matrix(data = S_0, nrow=n_levels+1, ncol=N_)
  #movements<-sample(x = c(u_,d_),size = N_*n_levels,replace=TRUE,prob=c(qu,qd))
  movements<-exp((r_-sigma^2/2)*h_+sigma*sqrt(h_)*rnorm(N_*n_levels,mean=0,sd=1))
  movements<-matrix(data = movements,nrow=n_levels,ncol=N_)
  for(ii in 1:n_levels){
    time<-ii+1
    paths$prices[time,]=paths$prices[time-1,]*movements[time-1,]
  }
  model_parameters<-list(T_=T_,h_=h_,r_=r_,sigma=sigma, method=method,u_=u_,d_=d_,qu=qu,qd=qd)
  return(list(paths,model_parameters))    
}

call_payoff_mc<-function(paths,param){
  K<-param[[1]]
  paths$payoffs<-apply(paths$prices,1:2,function(x)max(x-K,0))
  return(paths)
}

put_payoff_mc<-function(paths,param){
  K<-param[[1]]
  paths$payoffs<-apply(paths$prices,1:2,function(x)max(K-x,0))
  return(paths)
}

asian_call_payoff_mc<-function(paths,param){
  K<-param[[1]]
  numtoaverage<-min(param[[2]],nrow(paths$prices))
  
 
  #it would be more efficient to make more matrices than use rollapply
  averageprice<-rollapply(paths$prices[2:nrow(paths$prices),],width=numtoaverage,by.column=TRUE,partial=TRUE,FUN=mean,align="right")
  averageprice<-rbind(paths$prices[1,],averageprice)
  
  paths$payoffs<-apply(averageprice,1:2,function(x)max(x-K,0))
  
  return(paths)
}

asian_put_payoff_mc<-function(paths,param){
  K<-param[[1]]
  numtoaverage<-min(param[[2]],nrow(paths$prices))
  
  
  averageprice<-rollapply(paths$prices[2:nrow(paths$prices),],width=numtoaverage,by.column=TRUE,partial=TRUE,FUN=mean,align="right")
  averageprice<-rbind(paths$prices[1,],averageprice)

  
  paths$payoffs<-apply(averageprice,1:2,function(x)max(K-x,0))
  
  return(paths)
}



#american options and early exercise https://people.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf
#Longstaff–Schwartz (LS) algorithm as described in
#Valuing American Options by Simulation: A Simple Least-Squares Approach  (2001)
#note that although im generating binomial tree paths, this pricing function
#prices any path actually
price_mc_paths<-function(paths, is_european,payoff_function,param=NULL){
  model_parameters<-paths[[2]]
  paths<-paths[[1]]
  
  is_rc<-model_parameters$is_rc
  r_<-model_parameters$r_
  h_<-model_parameters$h_

  if(is_european){
    paths<-payoff_function(paths,param)
    paths$optionprices<-matrix(0,ncol=ncol(paths$prices),nrow=nrow(paths$prices))
    for(ii in nrow(paths$prices):1){#backwards for recursive pricing
      if(ii==nrow(paths$prices)){
        paths$optionprices[ii,]<-paths$payoffs[ii,]
      }else{
        paths$optionprices[ii,]<-exp(-r_*h_)*paths$optionprices[ii+1,]
      }
    }
  }else{
    paths<-payoff_function(paths,param)
    paths$cashflows<-matrix(0,ncol=ncol(paths$prices),nrow=nrow(paths$prices))
    for(ii in nrow(paths$prices):1){#backwards for recursive pricing
      if(ii==nrow(paths$prices)){
        paths$cashflows[ii,]<-paths$payoffs[ii,]
      }else{
        
        #following the LS algorithm paper:
        itmindices<-paths$payoffs[ii,]>0
        X<-paths$prices[ii,][itmindices]
        discounted_cashflows<-apply(paths$cashflows[ii:nrow(paths$cashflows),],2,function(x) sum(x*exp(-r_*h_*(ii:nrow(paths$cashflows)-ii))))
        Y<-discounted_cashflows[itmindices]
        if(length(Y)>0){
         
          #paper:.~One possible choice of basis functions is the set of 
          #(weighted) Laguerre polynomials
          #Thus we regress discounted realized cash flows on a  constant and three nonlinear functions of the stock price
          
          # LSlm<-lm(Y~1+I(exp(-X/2))+I(exp(-X/2)*(1-X))+I(exp(-X/2)*(1-2*X+X^2/2)))
          LSlm<-lm(Y~1+X+I(X^2))
          continuation_values<-LSlm$fitted.values
          
          exercise_values<-paths$payoffs[ii,]
          oeindices<-exercise_values[itmindices]>continuation_values
          oeindices<-which(itmindices)[oeindices]
          
          #zero the columns of columns that will be updated
          for(jj in (ii+1):nrow( paths$cashflows)){
            for(kk in oeindices){
              paths$cashflows[jj,kk]<-0
            }
          }
          paths$cashflows[ii,oeindices]<-exercise_values[oeindices]
         
        }
      }
    }
    #discount all the cashflows
    optionprices<-matrix(data=0, nrow=1,ncol=ncol(paths$cashflows))
    optionprices[1,]<-apply(paths$cashflows,2,function(x) sum(x*exp(-r_*h_*(1:nrow(paths$cashflows)-1))))
    paths$optionprices<-optionprices
  }
  return(list(paths,model_parameters))
}



example<-matrix(data=c(1,1.09,1.08,1.34,
                       1,1.16,1.26,1.54,
                       1,1.22,1.07,1.03,
                       1,0.93,0.97,0.92,
                       1,1.11,1.56,1.52,
                       1,0.76,0.77,0.9,
                       1,0.92,0.84,1.01,
                       1,0.88,1.22,1.34),ncol=8,nrow=4)
T_=1
h_=1
mcpaths2<-mc_paths(T_,h_,S_0,r_,sigma,method=2,10000)
ppaths2<-price_mc_paths(mcpaths2,is_european=FALSE,put_payoff_mc,list(K=K))
ppaths2[[1]]<-list(prices = example)
ppaths2[[2]]$r_=0.06
ppaths2[[2]]$h_=1
temp<-price_mc_paths(ppaths2,is_european=FALSE,put_payoff_mc,list(K=1.1))
mean(temp[[1]]$optionprices[1,])


a<-binomial_lattice(1/4,1/12,TRUE, 100, 0.05,0.2,1)
#print_binomial_lattice(a[[1]],"price")
a<-price_binomial_lattice(a,TRUE,call_payoff,list(K=100))
#print_binomial_lattice(a[[1]],"optionprice")

a1<-mc_paths(1/4,1/12,100,0.05,0.2,1,5000)
a1<-price_mc_paths(a1,TRUE,call_payoff_mc,list(K=100))
paste("mc: ", mean(a1[[1]]$optionprices[1,]),", bintree: ",a[[1]][[1]][[1]]$optionprice)

b<-binomial_lattice(1,1/3,FALSE, 200, 0.09,0.3,1)
#print_binomial_lattice(b[[1]],"price")
b<-price_binomial_lattice(b,TRUE,asian_put_payoff,list(K=200,numtoaverage=3))
#print_binomial_lattice(b[[1]],"optionprice")

b1<-mc_paths(1,1/3,200,0.09,0.3,1,5000)
b1<-price_mc_paths(b1,TRUE,asian_put_payoff_mc,list(K=200,numtoaverage=3))
paste("mc: ", mean(b1[[1]]$optionprices[1,]),", bintree: ",b[[1]][[1]][[1]]$optionprice)


c<-binomial_lattice(1,1/3,FALSE, 200, 0.09,0.3,1)
#print_binomial_lattice(c[[1]],"price")
c<-price_binomial_lattice(c,FALSE,asian_put_payoff,list(K=200,numtoaverage=3))
#print_binomial_lattice(c[[1]],"optionprice")

c1<-mc_paths(1,1/3,200,0.09,0.3,1,5000)
c1<-price_mc_paths(c1,FALSE,asian_put_payoff_mc,list(K=200,numtoaverage=3))
paste("mc: ", mean(c1[[1]]$optionprices[1,]),", bintree: ",c[[1]][[1]][[1]]$optionprice)



########example: pricing American put option

T_=1/12#time to expiration 
h_=T_*1/4#length of period 
S_0=100#itial stock price
r_=0.09#risk-free rate
sigma=0.3#volatility of stock
K=100#strike price

#American options generally require non-recombining trees

#Nonrecombining Binomial Tree, using method 1 which converges lognormal prices
bintree1<-binomial_lattice(T_,h_, is_rc=FALSE, 100, r_,sigma,method=1)
print_binomial_lattice(bintree1[[1]],"price")

#Nonrecombining Binomial Tree, using method 2 which also converges to lognormal prices
bintree2<-binomial_lattice(T_,h_, is_rc=FALSE, 100, r_,sigma,method=2)
print_binomial_lattice(bintree2[[1]],"price")

#note the slight differences. As a result they have slightly different risk neutral measures
#For large T_/h_ the stock prices converge to lognormal

#price the first using recursive pricing of the binomial tree
ptree1<-price_binomial_lattice(bintree1,is_european=FALSE,put_payoff,list(K=K))
print_binomial_lattice(ptree1[[1]],"optionprice")
pt1<-ptree1[[1]][[1]][[1]]$optionprice


#price the second using recursive pricing of the binomial tree
ptree2<-price_binomial_lattice(bintree2,is_european=FALSE,put_payoff,list(K=K))
print_binomial_lattice(ptree2[[1]],"optionprice")
pt2<-ptree2[[1]][[1]][[1]]$optionprice

#Large T_/h_ is not computational using non-recombining binomial trees
#There are 2^(T_/h_)-1 nodes in a non-recombining binomial tree

#Instead price using monte-carlo and Longstaff-Schwartz algorithm for American options, first method
mcpaths1<-mc_paths(T_,h_,S_0,r_,sigma,method=1,50000)
ppaths1<-price_mc_paths(mcpaths1,is_european=FALSE,put_payoff_mc,list(K=K))
pp1<-mean(ppaths1[[1]]$optionprices[1,])

#second method
mcpaths2<-mc_paths(T_,h_,S_0,r_,sigma,method=2,50000)
ppaths2<-price_mc_paths(mcpaths2,is_european=FALSE,put_payoff_mc,list(K=K))
pp2<-mean(ppaths2[[1]]$optionprices[1,])

#however with monte-carlo we can increase T_/h_ and still be reasonably computationally efficient
#we can do this and compare to online calculators, 
#so redoing 
h_=T_*1/50


mcpaths12<-mc_paths(T_,h_,S_0,r_,sigma,method=1,1000)
ppaths12<-price_mc_paths(mcpaths12,is_european=FALSE,put_payoff_mc,list(K=K))
pp12<-mean(ppaths12[[1]]$optionprices[1,])


mcpaths22<-mc_paths(T_,h_,S_0,r_,sigma,method=2,1000)
ppaths22<-price_mc_paths(mcpaths22,is_european=FALSE,put_payoff_mc,list(K=K))
pp22<-mean(ppaths22[[1]]$optionprices[1,])

qqnorm(log(ppaths12[[1]]$prices[nrow(ppaths12[[1]]$prices),]))#ultimate price agrees with normal quantiles

library(LSMonteCarlo)

pLSM<-AmerPutLSM(Spot = S_0, sigma = sigma, n = 1000, m = T_/h_, Strike = K, r = r_,dr = 0, mT = T_)[[1]]

print(paste(pt1, pt2, pp1, pp2, pp12,pp22,pLSM))
