rm(list = ls())
set.seed(42)

library(rjags)
library(coda)
library(parallel)
library(openxlsx)

#==================================================
# 0. 写入【纵向-多状态联合模型】JAGS代码
#==================================================
write_joint_ms_jags <- function() {
  jags_code <- "
model{
  # ========== 纵向标志物参数 ==========
  beta0_y ~ dnorm(0, 0.01)
  beta1_y ~ dnorm(0, 0.01)
  sigma_y ~ dunif(0, 5)
  tau_y <- pow(sigma_y, -2)

  # ========== 多状态转移基线风险 ==========
  lambda12 ~ dlnorm(0, 1)
  lambda13 ~ dlnorm(0, 1)
  lambda23 ~ dlnorm(0, 1)

  # ========== 协变量系数 + 纵向关联参数alpha ==========
  beta12 ~ dnorm(0, 1.0/4)
  beta13 ~ dnorm(0, 1.0/4)
  beta23 ~ dnorm(0, 1.0/4)
  alpha ~ dnorm(0, 0.01)

  # ========== 个体共享随机效应（纵向+多状态共用） ==========
  sigma_b ~ dunif(0, 5)
  tau_b <- pow(sigma_b, -2)
  for(i in 1:N){
    b[i] ~ dnorm(0, tau_b)
  }

  # ========== 纵向测量似然 ==========
  for(l in 1:n_long) {
    mu_y[l] <- beta0_y + beta1_y * time_long[l] + b[id_long[l]]
    y[l] ~ dnorm(mu_y[l], tau_y)
  }

  # ========== 多状态转移区间似然 ==========
  for(n in 1:n_trans){
    dt[n] <- end_time[n] - start_time[n]
    zi[n] <- z[id[n]]
    bi[n] <- b[id[n]]

    # Cox风险：加入纵向关联项 alpha*bi[n]
    h12[n] <- lambda12 * exp(beta12 * zi[n] + alpha * bi[n])
    h13[n] <- lambda13 * exp(beta13 * zi[n] + alpha * bi[n])
    h23[n] <- lambda23 * exp(beta23 * zi[n] + alpha * bi[n])

    is_from1[n] <- equals(from[n], 1)
    is_from2[n] <- equals(from[n], 2)
    is_to0[n]   <- equals(to[n], 0)
    is_to2[n]   <- equals(to[n], 2)
    is_to3[n]   <- equals(to[n], 3)

    # 状态1出发对数似然
    loglik1[n] <- is_from1[n] * (
      is_to2[n] * (log(h12[n]) - (h12[n] + h13[n]) * dt[n]) +
      is_to3[n] * (log(h13[n]) - (h12[n] + h13[n]) * dt[n]) +
      is_to0[n] * (-(h12[n] + h13[n]) * dt[n])
    )
    # 状态2出发对数似然
    loglik2[n] <- is_from2[n] * (
      is_to3[n] * (log(h23[n]) - h23[n] * dt[n]) +
      is_to0[n] * (-h23[n] * dt[n])
    )

    loglik[n] <- loglik1[n] + loglik2[n]
    phi[n]    <- -loglik[n] + C
    zeros[n]  ~ dpois(phi[n])
  }
}
"
writeLines(jags_code, "joint_multistate_long.jags")
}

#==================================================
# 1. 全局真实参数（纵向+多状态+关联参数整合）
#==================================================
true_par <- c(
  # 纵向标志物
  beta0_y = 2,
  beta1_y = 0.5,
  sigma_y = 0.8,
  # 多状态基线风险
  lambda12 = 1,
  lambda13 = 1,
  lambda23 = 1,
  # 转移协变量系数
  beta12 = -1,
  beta13 = 2,
  beta23 = 1,
  # 随机效应SD + 纵向-转移关联参数
  sigma_b = 0.5,
  alpha   = 0.8
)

#==================================================
# 2. 模拟函数：同时生成纵向数据 + 多状态转移数据
#==================================================
simulate_joint_ms <- function(
    N,
    par,
    censor_rate = 0.08,
    Tmax = 15,
    n_visit_max = 4
){
  # 提取真实参数
  beta0_y  <- par["beta0_y"]
  beta1_y  <- par["beta1_y"]
  sigma_y  <- par["sigma_y"]
  lambda12 <- par["lambda12"]
  lambda13 <- par["lambda13"]
  lambda23 <- par["lambda23"]
  beta12   <- par["beta12"]
  beta13   <- par["beta13"]
  beta23   <- par["beta23"]
  sigma_b  <- par["sigma_b"]
  alpha    <- par["alpha"]
  
  z <- rnorm(N)
  b <- rnorm(N, 0, sigma_b) # 个体共享随机效应
  censor_time <- rexp(N, rate = censor_rate)
  
  # 2.1 生成纵向标志物数据集
  long_list <- list()
  long_row <- 1
  for(i in 1:N){
    obs_end_i <- min(censor_time[i], Tmax)
    n_vis <- sample(2:n_visit_max, size = 1)
    visit_t <- sort(runif(n_vis, min = 0, max = obs_end_i))
    for(tv in visit_t){
      mu <- beta0_y + beta1_y * tv + b[i]
      y_obs <- mu + rnorm(1, 0, sigma_y)
      long_list[[long_row]] <- data.frame(
        id = i, time = tv, y = y_obs, z = z[i]
      )
      long_row <- long_row + 1
    }
  }
  long_data <- do.call(rbind, long_list)
  
  # 2.2 生成多状态转移区间数据集
  trans_list <- list()
  trans_row <- 1
  for(i in 1:N){
    t_cur <- 0
    state_cur <- 1
    obs_end <- min(censor_time[i], Tmax)
    while(t_cur < obs_end && state_cur != 3){
      if(state_cur == 1){
        h12 <- lambda12 * exp(beta12 * z[i] + alpha * b[i])
        h13 <- lambda13 * exp(beta13 * z[i] + alpha * b[i])
        total_h <- h12 + h13
        t_next <- t_cur - log(runif(1)) / total_h
        if(t_next > obs_end){
          trans_list[[trans_row]] <- data.frame(
            id = i, from = 1, to = 0, start = t_cur, end = obs_end, status = 0
          )
          trans_row <- trans_row + 1
          break
        }
        p_vec <- c(h12, h13) / total_h
        next_s <- sample(c(2, 3), size = 1, prob = p_vec)
      } else if(state_cur == 2){
        h23 <- lambda23 * exp(beta23 * z[i] + alpha * b[i])
        t_next <- t_cur - log(runif(1)) / h23
        if(t_next > obs_end){
          trans_list[[trans_row]] <- data.frame(
            id = i, from = 2, to = 0, start = t_cur, end = obs_end, status = 0
          )
          trans_row <- trans_row + 1
          break
        }
        next_s <- 3
      }
      # 记录有效转移
      trans_list[[trans_row]] <- data.frame(
        id = i, from = state_cur, to = next_s, start = t_cur, end = t_next, status = 1
      )
      trans_row <- trans_row + 1
      t_cur <- t_next
      state_cur <- next_s
    }
  }
  trans_data <- do.call(rbind, trans_list)
  
  return(list(
    N = N,
    z = z,
    b = b,
    long_data = long_data,
    trans_data = trans_data
  ))
}

#==================================================
# 3. JAGS模型拟合函数（同时输入纵向+多状态数据）
#==================================================
fit_joint_jags <- function(simdat, n_adapt = 100, n_burn = 300, n_iter = 800, C_const = 10000){
  long_df <- simdat$long_data
  trans_df <- simdat$trans_data
  N <- simdat$N
  
  jags_data <- list(
    # 纵向观测
    n_long = nrow(long_df),
    id_long = long_df$id,
    time_long = long_df$time,
    y = long_df$y,
    # 多状态转移区间
    N = N,
    n_trans = nrow(trans_df),
    from = trans_df$from,
    to = trans_df$to,
    id = trans_df$id,
    start_time = trans_df$start,
    end_time = trans_df$end,
    z = simdat$z,
    zeros = rep(0, nrow(trans_df)),
    C = C_const
  )
  
  # 初始值
  inits <- list(
    beta0_y = 2, beta1_y = 0.5, sigma_y = 1,
    lambda12 = 0.8, lambda13 = 0.8, lambda23 = 0.8,
    beta12 = 0, beta13 = 0, beta23 = 0,
    alpha = 0.5, sigma_b = 1,
    b = rnorm(N, 0, 0.1)
  )
  
  # 需要监控的全部参数
  monitor_pars <- c(
    "beta0_y", "beta1_y", "sigma_y",
    "lambda12", "lambda13", "lambda23",
    "beta12", "beta13", "beta23",
    "alpha", "sigma_b"
  )
  
  jm <- jags.model(
    file = "joint_multistate_long.jags",
    data = jags_data,
    inits = inits,
    n.chains = 3,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  update(jm, n_burn)
  samp <- coda.samples(jm, variable.names = monitor_pars, n.iter = n_iter, thin = 1)
  return(samp)
}

#==================================================
# 4. 提取后验均值、标准差、95%CI
#==================================================
extract_post_summary <- function(samp, rep_id, true_par){
  mat <- as.matrix(samp)
  post_mean <- colMeans(mat)
  post_sd   <- apply(mat, 2, sd)
  ci_mat    <- apply(mat, 2, quantile, probs = c(0.025, 0.975))
  
  df <- data.frame(
    Parameter = colnames(mat),
    Estimate  = post_mean,
    PostSD    = post_sd,
    Lower95   = ci_mat[1, ],
    Upper95   = ci_mat[2, ],
    Replicate = rep_id
  )
  # 绑定真实参数
  df$TrueValue <- true_par[df$Parameter]
  df$Cover <- as.integer(df$Lower95 <= df$TrueValue & df$TrueValue <= df$Upper95)
  return(df)
}

#==================================================
# 5. 单次蒙特卡洛重复：模拟+拟合+提取结果
#==================================================
one_rep <- function(rep_id, N, true_par){
  # 生成联合数据
  simdat <- simulate_joint_ms(N = N, par = true_par)
  # 拟合JAGS
  samp <- fit_joint_jags(simdat)
  # 后验汇总
  res_df <- extract_post_summary(samp, rep_id = rep_id, true_par = true_par)
  res_df$SampleSize <- N
  return(res_df)
}

#==================================================
# 6. 批量串行运行MC模拟
#==================================================
run_mc_serial <- function(n_rep, N, true_par){
  res_list <- vector("list", n_rep)
  for(r in 1:n_rep){
    cat("==== Rep", r, "/", n_rep, " N =", N, "====\n")
    res_list[[r]] <- tryCatch(
      one_rep(rep_id = r, N = N, true_par = true_par),
      error = function(e){
        cat("Rep", r, "failed, error:", e$message, "\n")
        return(NULL)
      }
    )
  }
  valid_idx <- !sapply(res_list, is.null)
  valid_df  <- do.call(rbind, res_list[valid_idx])
  return(list(
    detail = valid_df,
    ok_rep = sum(valid_idx),
    fail_rep = sum(!valid_idx)
  ))
}

#==================================================
# 7. 汇总评价指标：Bias、RBias、MSE、RMSE、覆盖率、MCSE
#==================================================
summarize_performance <- function(mc_obj){
  df_all <- mc_obj$detail
  if(nrow(df_all) == 0) stop("无有效重复结果")
  by_par <- split(df_all, df_all$Parameter)
  perf_list <- lapply(by_par, function(subdf){
    true_val <- unique(subdf$TrueValue)
    mean_est <- mean(subdf$Estimate)
    bias     <- mean_est - true_val
    rbias    <- bias / true_val
    mse      <- mean((subdf$Estimate - true_val)^2)
    rmse     <- sqrt(mse)
    cover    <- mean(subdf$Cover)
    mean_psd <- mean(subdf$PostSD)
    mcse     <- sd(subdf$Estimate) / sqrt(nrow(subdf))
    
    data.frame(
      Parameter     = unique(subdf$Parameter),
      True_Value    = round(true_val, 4),
      Mean_Estimate = round(mean_est, 4),
      Mean_PostSD   = round(mean_psd, 4),
      Bias          = round(bias, 4),
      Relative_Bias = round(rbias, 4),
      MSE           = round(mse, 4),
      RMSE          = round(rmse, 4),
      Coverage      = round(cover, 4),
      MCSE          = round(mcse, 5),
      SampleSize    = unique(subdf$SampleSize)
    )
  })
  perf_df <- do.call(rbind, perf_list)
  rownames(perf_df) <- NULL
  return(perf_df)
}

#==================================================
# 8. 主执行流程
#==================================================
# 1）写入JAGS模型文件
write_joint_ms_jags()

# 2）设置模拟参数
sample_sizes <- c(100, 300, 500) # 样本量可自行修改
n_rep_total   <- 100         # 重复次数，正式模拟改为100

all_perf_tables <- list()

# 3）循环不同样本量运行模拟
for(n in sample_sizes){
  cat("\n=====================\n")
  cat("开始 N =", n, "\n")
  cat("=====================\n")
  mc_out <- run_mc_serial(n_rep = n_rep_total, N = n, true_par = true_par)
  cat("成功重复数：", mc_out$ok_rep, "失败：", mc_out$fail_rep, "\n")
  
  # 保存该样本量原始明细
  write.csv(mc_out$detail, paste0("JointMS_Detail_N", n, ".csv"), row.names = F)
  
  # 计算性能汇总表
  perf_tab <- summarize_performance(mc_out)
  all_perf_tables[[as.character(n)]] <- perf_tab
  print(perf_tab)
  write.csv(perf_tab, paste0("JointMS_Perf_N", n, ".csv"), row.names = F)
}

# 合并所有样本量汇总结果
final_perf <- do.call(rbind, all_perf_tables)

# 导出全部汇总结果
write.csv(final_perf, "JointMS_All_Performance.csv", row.names = F)
saveRDS(final_perf, "JointMS_All_Performance.rds")

# Excel输出
wb <- createWorkbook()
addWorksheet(wb, "Performance Summary")
writeData(wb, 1, final_perf)
# 写入真实参数表
par_df <- data.frame(Param = names(true_par), True = as.numeric(true_par))
addWorksheet(wb, "True Parameters")
writeData(wb, 2, par_df)
saveWorkbook(wb, "Joint_MultiState_Long_Sim_Results.xlsx", overwrite = TRUE)

cat("\n模拟全部完成！文件已保存至工作目录\n")
print(final_perf)
