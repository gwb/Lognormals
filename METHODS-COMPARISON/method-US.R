




# Numerically unstable versions of the f1, f2, f3 and their derivatives
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

f_u = function(w, mu, sigsq){
  return(-1/(2*sigsq)*(w-mu)^2 + log( log( 1+exp(w) ) ) - 1/2 * log(2*pi*sigsq) )
  }


f2_u = function(w, mu, sigsq){
  return(-1/(2*sigsq)*(w-mu)^2 + w - log( 1+exp(w) ) - 1/2 * log(2*pi*sigsq) )
  }

f3_u = function(w, mu, sigsq){
  return(-1/(2*sigsq)*(w-mu)^2 + log( log( 1+exp(w) )^2 ) - 1/2 * log(2*pi*sigsq) )
  }


df_u = function(w, mu, sigsq){
  return( -(w-mu)/sigsq + exp(w)/( (1+exp(w)) * log(1+exp(w))) )
}

df2_u = function(w, mu, sigsq){
  return( 1 - exp(w)/( (1+exp(w)) ) -(w-mu)/sigsq)
}

df3_u = function(w, mu, sigsq){
  return( -(w-mu)/sigsq + 2*exp(w)/( (1+exp(w)) * log(1+exp(w))) )
}

ddf_u = function(w, mu, sigsq){
  t1 = -1/sigsq
  t2 = - exp(2*w)/( (1+exp(w))^2 * log(1+exp(w))^2 )
  t3 = - exp(2*w)/( (1+exp(w))^2 * log(1+exp(w)) )
  t4 = exp(w)/( (1+exp(w)) * log(1+exp(w)) )
  return(t1 + t2 + t3 + t4)
}

ddf2_u = function(w, mu, sigsq){
  t1 = exp(2*w)/( (1+exp(w))^2 )
  t2 = - exp(w)/( (1+exp(w)) )
  t3 = -1/sigsq
  return(t1 + t2 + t3)
}

ddf3_u = function(w, mu, sigsq){
  t1 = -1/sigsq
  t2 = - 2*exp(2*w)/( (1+exp(w))^2 * log(1+exp(w))^2 )
  t3 = - 2*exp(2*w)/( (1+exp(w))^2 * log(1+exp(w)) )
  t4 = 2*exp(w)/( (1+exp(w)) * log(1+exp(w)) )
  return(t1 + t2 + t3 + t4)
}

dddf_u = function(w, mu, sigsq){
  t1 = 2*exp(3*w)/( (1+exp(w))^3 * log(1+exp(w))^3 )
  t2 = 3*exp(3*w)/( (1+exp(w))^3 * log(1+exp(w))^2 )
  t3 = -3*exp(2*w)/( (1+exp(w))^2 * log(1+exp(w))^2 )
  t4 = 2*exp(3*w)/( (1+exp(w))^3 * log(1+exp(w)) )
  t5 = -3*exp(2*w)/( (1+exp(w))^2 * log(1+exp(w)) )
  t6 = exp(w)/( (1+exp(w)) * log(1+exp(w)) )
  return(t1 + t2 + t3 + t4 + t5 + t6)
}

dddf2_u = function(w, mu, sigsq){
  t1 = - 2*exp(3*w)/( (1+exp(w))^3 )
  t2 = 3*exp(2*w)/( (1+exp(w))^2 )
  t3 = - exp(w)/( (1+exp(w)) )
  return(t1 + t2 + t3 )
}

dddf3_u = function(w, mu, sigsq){
  t1 = 4*exp(3*w)/( (1+exp(w))^3 * log(1+exp(w))^3 )
  t2 = 6*exp(3*w)/( (1+exp(w))^3 * log(1+exp(w))^2 )
  t3 = -6*exp(2*w)/( (1+exp(w))^2 * log(1+exp(w))^2 )
  t4 = 4*exp(3*w)/( (1+exp(w))^3 * log(1+exp(w)) )
  t5 = -6*exp(2*w)/( (1+exp(w))^2 * log(1+exp(w)) )
  t6 = 2*exp(w)/( (1+exp(w)) * log(1+exp(w)) )
  return(t1 + t2 + t3 + t4 + t5 + t6)
}


# Stable versions of the f1, f2, f3 and their derivatives
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

f = function(w, mu, sigsq){
  res = f_u(w,mu,sigsq)
  if(!is.infinite(res) && !is.nan(res)){
    return(res)
  }else{
    if(w > 0){
      res_alt_pos = -1/(2*sigsq)*(w-mu)^2 + log(w) - 1/2 * log(2*pi*sigsq)
      return(res_alt_pos)
    }else{
      res_alt_neg = -1/(2*sigsq)*(w-mu)^2 + w - 1/2 * log(2*pi*sigsq)
      return(res_alt_neg)
    }
  }
}

f2 = function(w, mu, sigsq){
  res = f2_u(w,mu,sigsq)
  res_alt = -1/(2*sigsq)*(w-mu)^2 - 1/2 * log(2*pi*sigsq)
  return( ifelse( !is.infinite(res), res, res_alt) )
  }


f3 = function(w, mu, sigsq){
  res = f3_u(w,mu,sigsq)
  res_alt_pos = -1/(2*sigsq)*(w-mu)^2 + log(w^2) - 1/2 * log(2*pi*sigsq)
  res_alt_neg = -1/(2*sigsq)*(w-mu)^2 + 2*w - 1/2 * log(2*pi*sigsq)
  if(!is.infinite(res) && !is.nan(res)){
    return(res)
  }else{
    return(ifelse(w>0, res_alt_pos,res_alt_neg))
  }
}

df = function(w, mu, sigsq){
  res = df_u(w,mu,sigsq)
  res_alt_pos = -(w-mu)/sigsq
  res_alt_neg = -(w-mu)/sigsq + 1/(1+exp(w))
  if(!is.nan(res) && !is.infinite(res)){
    return(res)
  }else{
    return(ifelse(w > 0, res_alt_pos, res_alt_neg))
  }
}

df2 = function(w, mu, sigsq){
  res = df2_u(w,mu,sigsq)
  return( ifelse( !is.nan(res), res, -(w-mu)/sigsq) )
}

df3 = function(w, mu, sigsq){
  res = df3_u(w,mu,sigsq)
  res_alt_pos = -(w-mu)/sigsq
  res_alt_neg = -(w-mu)/sigsq + 2/(1+exp(w))
  if(!is.nan(res) && !is.infinite(res)){
    return(res)
  }else{
    return(ifelse(w > 0, res_alt_pos, res_alt_neg))
  }
}

ddf_old = function(w, mu, sigsq){
  t = -1/sigsq
  res = ddf_u(w,mu,sigsq)
  return(ifelse( !is.nan(res), res, t)) 
}

ddf = function(w, mu, sigsq){
  res = ddf_u(w,mu,sigsq)
  res_alt_pos = -1/sigsq
  res_alt_neg = -1/sigsq - 1/(1+exp(w))^2 - exp(w)/(1+exp(w))^2 + 1/(1+exp(w))
  return(ifelse( !is.nan(res), res, ifelse(w > 0, res_alt_pos, res_alt_neg))) 
}

ddf2 = function(w, mu, sigsq){
  t = -1/sigsq
  res = ddf2_u(w,mu,sigsq)
  return(ifelse( !is.nan(res), res, t)) 
}

ddf3_old = function(w, mu, sigsq){
  t = -1/sigsq
  res = ddf3_u(w,mu,sigsq)
  return(ifelse( !is.nan(res), res, t)) 
}

ddf3 = function(w, mu, sigsq){
  res = ddf3_u(w,mu,sigsq)
  res_alt_pos = -1/sigsq
  res_alt_neg = -1/sigsq - 2/(1+exp(w))^2 - 2*exp(w)/(1+exp(w))^2 + 2/(1+exp(w))
  return(ifelse( !is.nan(res), res, ifelse(w > 0, res_alt_pos, res_alt_neg)))
}

dddf = function(w, mu, sigsq){
  res = dddf_u(w,mu,sigsq)
  return( ifelse( !is.nan(res), res, 0))
}

dddf2 = function(w, mu, sigsq){
  res = dddf2_u(w,mu,sigsq)
  return( ifelse( !is.nan(res), res, 0) )
}

dddf3 = function(w, mu, sigsq){
  res = dddf3_u(w,mu,sigsq)
  return(ifelse( !is.nan(res), res, 0))
}

# Convenient Structures
# # # # # # # # # # # #

fs_u = c(f_u, f2_u, f3_u)
dfs_u = c(df_u, df2_u, df3_u)
ddfs_u = c(ddf_u, ddf2_u, ddf3_u)
dddfs_u = c(dddf_u, dddf2_u, dddf3_u)

fs = c(f, f2, f3)
dfs = c(df, df2, df3)
ddfs = c(ddf, ddf2, ddf3)
dddfs = c(dddf, dddf2, dddf3)


# Halley's approximations
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

bisect = function(dfn, mu, sigsq, err=5){
  l = 0
  u = 0
  while(dfn(l,mu,sigsq) < 0){
    l = l - 10
    }
  while(dfn(u,mu,sigsq) > 0){
    u = u + 10
    }
  diff = 100000
  while(diff > err){
    mid = (u+l)/2
    if(dfn(mid, mu, sigsq) == 0){
        return(mid)
      }else{
        if(dfn(mid,mu,sigsq) > 0){
          l=mid
        }else{
          u=mid
        }
        diff = (u-l)/2
      }
  }
  return(mid)
}

# numerically unstable version
get_wargmax_u = function(w0, fn_id, mu, sigsq, eps=0.01){
  diff = eps+1
  ow = w0
  nw = ow - (2*dfs[[fn_id]](ow, mu, sigsq)*ddfs[[fn_id]](ow,mu,sigsq))/( 2*ddfs[[fn_id]](ow,mu,sigsq)^2 - dfs[[fn_id]](ow,mu,sigsq)*dddfs[[fn_id]](ow,mu,sigsq) )
  while(abs(nw - ow)  > eps){
    ow = nw
    nw = ow - (2*dfs[[fn_id]](ow, mu, sigsq)*ddfs[[fn_id]](ow,mu,sigsq))/( 2*ddfs[[fn_id]](ow,mu,sigsq)^2 - dfs[[fn_id]](ow,mu,sigsq)*dddfs[[fn_id]](ow,mu,sigsq) )
  }
  return(nw)
}

# numerically stable version
get_wargmax = function(fn_id, mu,sigsq, eps=0.01){
  v0 = bisect(dfs[[fn_id]], mu,sigsq)
  return(get_wargmax_u(v0, fn_id,mu,sigsq,eps))
}


# Approximation of the moments
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

get_ci = function(i, wargmax, mu, sigsq){
  return(sqrt( 2*pi/abs(ddfs[[i]](wargmax,mu,sigsq)) ) * exp(fs[[i]](wargmax,mu,sigsq)))
}

get_ftilde = function(fn_id, mu, sigsq){
  wargmax = get_wargmax(fn_id, mu, sigsq)
  ci = get_ci(fn_id, wargmax,mu, sigsq)
  
  ftilde = function(w){
    return(ci * (sqrt(abs(ddfs[[fn_id]](wargmax,mu,sigsq)))/sqrt(2*pi))*exp(-abs(ddfs[[fn_id]](wargmax,mu,sigsq)) *(w-wargmax)^2/2))
    }
  return(ftilde)
}

EZ = function(mu_y1, mu, sigsq){
  wargmax = get_wargmax(1, mu,sigsq)
  c1 = get_ci(1,wargmax, mu, sigsq)
  return(c1 + mu_y1)
  }

SIGSQZ = function(mu1, sigsq1, mu2, sigsq2){
  wargmax = sapply(c(1,2,3), get_wargmax, mu2-mu1, sigsq1+sigsq2)
  ci = sapply(c(1,2,3), function (x) get_ci(x,wargmax[x],mu2-mu1, sigsq2+sigsq1))
  return(sigsq1 + mu1^2 + 2*mu1*ci[1] - 2*sigsq1/(sigsq1+sigsq2)*(sigsq1+sigsq2)*ci[2] + ci[3] - (mu1 + ci[1])^2 )
}


get.US <- function(u1, u2, s1, s2){
  if (u2 > u1){
    tmp <- u1
    u1 <- u2
    u2 <- tmp

    tmp <- s1
    s1 <- s2
    s2 <- tmp
  }
  
  ez <- EZ(u1, u2 - u1, s1^2 + s2^2)
  ssq <- SIGSQZ(u1, s1^2, u2, s2^2)
  return(c(ez, sqrt(ssq)))
}
