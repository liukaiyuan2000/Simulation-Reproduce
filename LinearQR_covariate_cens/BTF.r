BTF = function(y, delta)
{
  # the weight of Bang-Tsiatis (2003) estimator
  # y: observed survival times
  # delta: censoring indicator
  
  #weight is an internal function to calculate delta/G(yi)
  ind = order(y) #对y排序的下标
  y = y[ind] #对y升序排列
  delta = delta[ind] #对delta以y的顺序排列
  n = length(y)
  
  s = numeric(length(y)) #s是n维向量
  l = survfit(Surv(y, 1 - delta) ~ 1)
  s[duplicated(y) == F] = l$surv #若y中没有重复变量，则s=l$surv
  index = seq(1, length(y), 1) #index=1,2,3,...,n
  s[duplicated(y) == T] = s[index[duplicated(y) == T] - 1] #若y中有重复变量，则为l$surv的上一个值
  
  if (s[length(y)] == 0) #若s(n)=0则赋给s(n)=1,避免出现NAN
  { s[length(y)] = 1 }   #assign any nonzero value to s[n] to avoid NaN
  
  for (j in 1:(n - 1)) 
  { 
    if (s[j + 1] == 0) 
    {
      s[j + 1] = s[j] #把s中为0的项数赋值为前一项
    }
  }
  weight = delta / s  #避免分母为0
  weight2 = rep(0, n) 
  weight2[ind] = weight #重排序，回到最初的delta的顺序以对应y的删失情况
  return(weight2)
}
