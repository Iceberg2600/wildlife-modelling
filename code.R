# 环境搭建
# --------------------------------------------------------
library(statsecol)
library(tidyverse)
library(zoo)
library(jagsUI)
library(MCMCvis)

# 1.数据预处理
# --------------------------------------------------------
# 加载数据
data("wildebeest")
str(wildebeest)

# 提取数据
years <- wildebeest$year
N_obs <- wildebeest$Nhat
se_obs <- wildebeest$sehat
rainfall <- wildebeest$rain
catch <- wildebeest$Catch

# 只对初始缺失值进行后向填充(back-filling)
N_init <- na.locf(N_obs, fromLast = TRUE)[1]

# 创建观测数据存在的年份索引
obs_years <- which(!is.na(N_obs))

# 过程模型数据（所有年份）
process_data <- list(
  rain = rainfall,
  catch = catch,
  n = length(years)
)

# 观测模型数据（只包含有观测值的年份）
observation_data <- list(
  y = N_obs[obs_years],        # 观测值
  se_y = se_obs[obs_years],    # 标准误
  n_obs = length(obs_years),   # 观测值数量
  obs_years = obs_years        # 观测年份的索引
)

# 合并数据用于JAGS
jags_data <- c(process_data, observation_data, list(N_init = N_init))

# 2.JAGS模型 - 离散时间密度无关指数增长
# --------------------------------------------------------
# I.Specify model in BUGS language
sink("wildebeestSSM.txt")
cat("
model{ #tell JAGS this is the model

  # Priors and constraints
  beta0 ~ dnorm(-2, 10)        
   beta1 ~ dlnorm(-1, 1)          
  sigma_proc ~ dunif(0.01, 0.1)      
  
  # initial state
  N[1] <- N_init
  
  # Likelihood - State process
  for(t in 1:(n-1)) {
    # Calculating growth rates r[t]
    log_r[t] <- beta0 + beta1 * rain[t]
    r[t] <- exp(log_r[t])
    
    # Expected value calculation (with protection against negative values)
    exp_N[t] <- N[t] + r[t] * N[t] - catch[t]
    
    # State equation with process error
    log_N[t+1] ~ dnorm(log(max(exp_N[t], 0.01)), pow(sigma_proc, -2))
    N[t+1] <- exp(log_N[t+1])
  }
  
  # Likelihood - Observation process
  for(i in 1:n_obs) {
    y[i] ~ dnorm(N[obs_years[i]], 1/(se_y[i]^2))
  }

}
",fill=TRUE)
sink()

# II.Collect and package data
jags_data

# III.Set initial values
inits <- function() {
  list(
    beta0 = rnorm(1, -2, 0.1),  
    beta1 = rlnorm(1, -1, 0.1),  
    sigma_proc = runif(1, 0.02, 0.05)
  )
}

# IV.Specify parameters to monitor
parms <- c("beta0","beta1","r","sigma_proc","N")

# V.Define MCMC settings
nc <- 3
nb <- 5000
ni <- 50000 + nb
nt <- 5

# VI.Model fitting
out <- jags(data = jags_data,
            inits = inits,
            parameters.to.save = parms,
            model.file = "wildebeestSSM.txt",
            n.chains = nc,
            n.iter = ni,
            n.burnin = nb,
            n.thin = nt)

# 3.评估参数收敛性
# --------------------------------------------------------
MCMCtrace(out,                
          params = parms[c(1, 2, 4)], 
          ind = TRUE,              
          pdf = FALSE)             

MCMCsummary(out,                
          params = parms[c(1, 2, 4)]
)


# 4.图像--拟合的种群轨迹，包括不确定性，并添加观测数据 
# --------------------------------------------------------
wildebeest_traj <- data.frame(Year = wildebeest$year,
                          Mean = out$mean$N,
                          Lower = out$q2.5$N,
                          Upper = out$q97.5$N,
                          Obs = wildebeest$Nhat)

ggplot(data = wildebeest_traj) + 
  geom_ribbon(aes(x=Year, y=Mean, ymin=Lower, ymax=Upper),
              fill="grey", alpha = 0.25) +
  geom_line(aes(x=Year, y=Mean), linewidth=1, color="blue") + 
  geom_line(aes(x=Year, y=Obs)) +
  geom_point(aes(x=Year, y=Obs), size=1.2) +
  theme_bw()

# 5.预测
# --------------------------------------------------------
# 计算平均降雨量
mean_rain <- mean(rainfall, na.rm = TRUE)
# 计算平均观测误差
mean_se <- mean(se_obs, na.rm = TRUE)
# 获取最后观测年份的捕获量作为未来捕获水平
last_catch <- tail(catch, 1)

# number of projection years
nproj <- 5 # data goes to 1989, we want to project to 1994

# new data
process_data_project <- list(
  rain = c(rainfall,rep(mean_rain,5)),
  catch = c(catch,rep(last_catch,5)),
  n = length(c(years,1990:1994))
)

observation_data_project <- list(
  y = c(N_obs[obs_years],rep(NA, nproj)),        
  se_y = c(se_obs[obs_years],rep(mean_se, nproj)),
  n_obs = length(c(obs_years,31:35)),  
  obs_years = c(obs_years,31:35)
)

wildebeest_project <- c(process_data_project, observation_data_project, list(N_init = N_init))

wildebeestproj <- jags(data = wildebeest_project, #<- CHANGE
                   inits = inits,
                   parameters.to.save = parms,
                   model.file = "wildebeestSSM.txt",
                   n.chains = nc,
                   n.iter = ni,
                   n.burnin = nb,
                   n.thin = nt)

wildebeestproj_traj <- data.frame(Year = c(years,1990:1994),
                              Mean = wildebeestproj$mean$N,
                              Lower = wildebeestproj$q2.5$N,
                              Upper = wildebeestproj$q97.5$N,
                              Obs = c(wildebeest$Nhat,rep(NA, nproj)))

ggplot(data = wildebeestproj_traj) + 
  geom_ribbon(aes(x=Year, y=Mean, ymin=Lower, ymax=Upper),
              fill="grey", alpha = 0.25) +
  geom_line(aes(x=Year, y=Mean), linewidth=1, color="blue") + 
  geom_line(aes(x=Year, y=Obs)) +
  geom_point(aes(x=Year, y=Obs), size=1.2) +
  theme_bw()




