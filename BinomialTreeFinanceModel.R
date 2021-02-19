#Generates a binomial lattice for stock price
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
binomial_lattice<-function(T_,h_, is_rc, S_0, r_, sigma,method){
  n_levels<-T_/h_#number of levels in the lattice
  if(method==1){
    u_<-exp(sigma*sqrt(h_))
    d_<-1/u_
  }else{
    delta<-exp(r_*h_)*sqrt(e^(sigma^2*h_)-1)
    u_<-exp(r_*h_)+delta
    d_<-exp(r_*h_)-delta
  }
  lattice<-vector("list",length = n_levels+1)
  
  for(ii in 1:(n_levels+1)){
    level<-ii-1# note i'm doing this because of R not being zero-indexed
    if(ii==1){
      row=list(list(price = S_0, T_=T_,h_=h_,is_rc=is_rc,r_=r_,sigma=sigma, method=method))# set initial stock price to S_0,
      #also store parameters here for convenience
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
  return(lattice)
}



call_payoff<-function(lattice){
  #the idea is to return a lattice but with an additional payoff parameter in the final layer. 
}


price_binomial_lattice<-function(lattice, is_european,payoff_function){
  #if european then call the payoff on the final time period in the lattice, otherwise every time period
  #if european then use the fact that discounted price is a martingale under risk neutral measure
  #if american use the fact above but also consider the option of exercising at all time periods
  #the idea is to return a lattice but with an additional option price parameter in the final layer. 
}

#this print function leaves a lot to be desired but eh
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

a<-binomial_lattice(3,1,TRUE, 100, 0.05,0.2,1)
print_binomial_lattice(a,"price")

#monte-carlo implementations - does not compute the entire lattice but rather some simulated paths
#this is useful for path-dependent options due to otherwise exponential lattice size, for n periods
#a non-recombining tree has O(2^n) nodes and a recombining tree has O(n^2) nodes. It is still useful
#for recombining trees and path-independent options but less so. 




