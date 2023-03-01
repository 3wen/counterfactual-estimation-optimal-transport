# Univariate case----

## Models----

### GAM----

library(splines)

#' Returns $\hat{m}_0()$ and $\hat{m}_1()$, using GAM
#' @param target name of the target variable
#' @param treatment_name name of the treatment variable (column in `data`)
#' @param x_m_name name of the mediator variable
#' @param data data frame with the observations
#' @param treatment_0 value for non treated
#' @param treatment_1 value for treated
#' @param df degrees of freedom (in the `bs()` function, deefault to `df=3`)
models_spline_univ <- function(target, treatment_name, x_m_name, data, treatment_0, treatment_1, df = 3){
  # \hat{m}_0()
  reg_0 <-
    bquote(glm(.(target)~bs(.(x_m_name), df = df), data=base,
               family=binomial, subset = (.(treatment_name) == .(treatment_0))),
           list(target = as.name(target), x_m_name=as.name(x_m_name),
                treatment_name = as.name(treatment_name),
                treatment_0 = treatment_0)) %>% 
    eval()
  # \hat{m}_1()
  reg_1 <- bquote(glm(.(target)~bs(.(x_m_name), df = df), data=base,
                      family=binomial, subset = (.(treatment_name) == .(treatment_1))),
                  list(target = as.name(target), x_m_name=as.name(x_m_name),
                       treatment_name = as.name(treatment_name),
                       treatment_1 = treatment_1)) %>% 
    eval()
  
  list(reg_0 = reg_0, reg_1 = reg_1)
}

### Kernels----

#' @param target name of the response variable
#' @param x_m_name name of the mediator variable
#' @param x single value which will be transported
#' @param h bandwidth
#' @param data data frame, subset of data with t=0 or t=1
pred_kernel_single <- function(target, x_m_name, x, h=10, data){
  w <- dnorm(pull(data, x_m_name), x, h)
  sum(w*pull(data, target))/sum(w)
}

## Prediction Functions----

### GAM----

#' @param object regression model (GAM)
#' @param newdata data frame in which to look for the mediator variable used to predict the target
model_spline_predict <- function(object, newdata){
  predict(object, newdata = newdata, type="response")
}

### Kernels----

#' @param target name of the response variable
#' @param x_m_name name of the mediator variable
#' @param x single value which will be transported
#' @param h bandwidth
#' @param data data frame, subset of data with t=0 or t=1
pred_kernel_single <- function(target, x_m_name, x, h=10, data){
  w <- dnorm(pull(data, x_m_name), x, h)
  sum(w*pull(data, target))/sum(w)
}


#' @param object regression model (GAM)
#' @param newdata data frame in which to look for the mediator variable used to predict the target
pred_kernel <- function(object, newdata){
  target <- object$target
  x_m_name <- object$x_m_name
  h <- object$h
  data <- object$data
  pull(newdata, 1) %>% 
    map_dbl(~pred_kernel_single(target = target, x_m_name = x_m_name, x = ., h = h,
                                data = data))
}

## Transport----

### Quantile-based----

#' Quantile-based transport (Univariate)
#' @param x_0 vector of values to be transported
#' @param x_m_name name of the mediator variable
#' @param treatment_name name of the treatment variable (column in `data`)
#' @param treatment_0 value for non treated
#' @param treatment_1 value for treated
#' @param data data frame with both T=0 or T=1
transport_quantile <- function(x_0, x_m_name, treatment_name, treatment_0, treatment_1, data){
  ind_0 <- pull(data, treatment_name) == treatment_0
  x_val <- pull(data, x_m_name)
  # Empirical Cumulative Distribution Function for values of the mediator variable for non-treated
  Fn <- ecdf(x_val)
  # Probability associated in the non-treated
  u <- Fn(x_0)
  # Transported_values
  x_1 <- pull(data, x_m_name)[pull(data, treatment_name) == treatment_1]
  x_t_quantile <- quantile(x_1, u)
  list(x_0 = x_0, u = u, x_t_quantile = x_t_quantile)
}

### Gaussian assumption----

#' @param x_0 vector of values to be transported
#' @param x_m_name name of the mediator variable
#' @param treatment_name name of the treatment variable (column in `data`)
#' @param treatment_0 value for non treated
#' @param treatment_1value for treated
#' @param data data frame with both T=0 or T=1
transport_univ_gaussian <- function(x_0, x_m_name, treatment_name, treatment_0, treatment_1, data){
  x0 <- mean(pull(data, x_m_name)[pull(data, treatment_name) == treatment_0])
  x1 <- mean(pull(data, x_m_name)[pull(data, treatment_name) == treatment_1])
  s0 <- sd(pull(data, x_m_name)[pull(data, treatment_name) == treatment_0])
  s1 <- sd(pull(data, x_m_name)[pull(data, treatment_name) == treatment_1])
  u_N <- pnorm(x_0,x0,s0)
  x_t_N <- qnorm(u_N, x1, s1)
  list(x_0 = x_0, u_N = u_N, x_t_N = x_t_N)
}

# Multivariate case----

## Models----

### GAM----

library(splines)

#' Returns $\hat{m}_0()$ and $\hat{m}_1()$, using GAM
#' @param target name of the target variable
#' @param treatment_name name of the treatment variable (column in `data`)
#' @param x_m_names names of the mediator variable
#' @param data data frame with the observations
#' @param treatment_0 value for non treated
#' @param treatment_1 value for treated
#' @param df degrees of freedom (in the `bs()` function, deefault to `df=3`)
models_spline <- function(target, treatment_name, x_m_names, data, treatment_0, treatment_1, df = 3){
  
  # \hat{m}_0()
  formula_glm <- paste0(target, "~", paste0("bs(", x_m_names, ", df = ",df, ")", collapse = " + "))
  reg_0 <- bquote(
    glm(formula = .(formula_glm), data=data, family = binomial,
        subset = (.(treatment_name) == .(treatment_0))),
    list(
      formula_glm = formula_glm,
      treatment_name = as.name(treatment_name),
      treatment_0 = treatment_0)
  ) %>% eval()
  # \hat{m}_1()
  reg_1 <- bquote(
    glm(formula = .(formula_glm), data=data, family = binomial,
        subset = (.(treatment_name) == .(treatment_1))),
    list(
      formula_glm = formula_glm,
      treatment_name = as.name(treatment_name),
      treatment_1 = treatment_1)
  ) %>% eval()
  
  list(reg_0 = reg_0, reg_1 = reg_1)
}

## Prediction Functions----

### GAM----

#' @param object regression model (GAM)
#' @param newdata data frame in which to look for the mediator variable used to predict the target
model_spline_predict <- function(object, newdata){
  predict(object, newdata = newdata, type="response")
}

## Transport----

#' Optimal Transport assuming Gaussian distribution for the mediator variables (helper function)
#' @return A list with the mean and variance of the mediator in each subset, and the symmetric matrix A.
#' @param target name of the target variable
#' @param x_m_names vector of names of the mediator variables
#' @param scale vector of scaling to apply to each `x_m_names` variable to transport (default to 1)
#' @param treatment_name name of the treatment variable (column in `data`)
#' @param treatment_0 value for non treated
#' @param treatment_1 value for treated
transport_gaussian_param <- function(target, x_m_names, scale=1, treatment_name, treatment_0, treatment_1){
  base_0_unscaled <- 
    base %>% 
    select(!!c(x_m_names, treatment_name)) %>% 
    filter(!!sym(treatment_name) ==  treatment_0)
  base_1_unscaled <- 
    base %>% 
    select(!!c(x_m_names, treatment_name)) %>% 
    filter(!!sym(treatment_name) ==  treatment_1)
  
  for(i in 1:length(scale)){
    base_0_scaled <- 
      base_0_unscaled %>% 
      mutate(!!sym(x_m_names[i]) := !!sym(x_m_names[i]) * scale[i])
    base_1_scaled <- 
      base_1_unscaled %>% 
      mutate(!!sym(x_m_names[i]) := !!sym(x_m_names[i]) * scale[i])
  }
  
  # Mean in each subset (i.e., T=0, T=1)
  m_0 <- base_0_scaled %>% summarise(across(!!x_m_names, mean)) %>% as_vector()
  m_1 <- base_1_scaled %>% summarise(across(!!x_m_names, mean)) %>% as_vector()
  # Variance
  S_0 <- base_0_scaled %>% select(!!x_m_names) %>% var()
  S_1 <- base_1_scaled %>% select(!!x_m_names) %>% var()
  # Matrix A
  A <- (solve(sqrtm(S_0))) %*% sqrtm( sqrtm(S_0) %*% S_1 %*% (sqrtm(S_0)) ) %*% solve(sqrtm(S_0))
  
  list(m_0 = m_0, m_1 = m_1, S_0 = S_0, S_1 = S_1, A = A, scale = scale, x_m_names = x_m_names)
}


#' Gaussian transport for a single observation (helper function)
#' @param z vector of variables to transport (for a single observation)
#' @param A symmetric positive matrix that satisfies \(A\sigma_0A=\sigma_1\)
#' @param m_0,m_1 vectors of mean values in the subsets \(\mathcal{D}_0\) and \(\mathcal{D}_1\)
#' @param scale vector of scaling to apply to each variable to transport (default to 1)
T_C_single <- function(z, A, m_0, m_1, scale = 1){
  z <- z*scale
  as.vector(m_1 + A %*% (z-m_0))*(1/scale)
}

#' Gaussian transport
#' @param z data frame of variables to transport
#' @param params (optional) parameters to use for the transport (result of `transport_gaussian_param()`, default to NULL)
#' @param target name of the target variable
#' @param x_m_names vector of names of the mediator variables
#' @param scale vector of scaling to apply to each `x_m_names` variable to transport (default to 1)
#' @param treatment_name name of the treatment variable (column in `data`)
#' @param treatment_0 value for non treated
#' @param treatment_1 value for treated
gaussian_transport <- function(z, params = NULL, ..., x_m_names, scale, treatment_name, treatment_0, treatment_1){
  if(is.null(params)){
    # If the parameters of the transport function are not provided
    # they need to be computed first
    params <- 
      transport_gaussian_param(target = target, x_m_names = x_m_names, 
                               scale = scale, treatment_name = treatment_name, 
                               treatment_0 = treatment_0, treatment_1 = treatment_1)
  }
  A <- params$A
  m_0 <- params$m_0 ; m_1 <- params$m_1
  scale <- params$scale
  x_m_names <- params$x_m_names
  
  values_to_transport <- z %>% select(!!x_m_names)
  transported_val <- 
    apply(values_to_transport, 1, T_C_single, A = A, m_0 = m_0, m_1 = m_1,
          scale = scale, simplify = FALSE) %>% 
    do.call("rbind", .)
  colnames(transported_val) <- colnames(z)
  transported_val <- as_tibble(transported_val)
  
  structure(.Data = transported_val, params = params)
}

# Estimation of CATE and SCATE----

#' Computes the Sample Average Treatment Effect with and without transport
#' @param x_0 vector of values at which to compute $\hat{m}_0(x_0)$, $\hat{m}_1(x_0)$
#' @param x_t vector of values at which to compute $\hat{m}_1(\mathcal{T}(x_t))$
#' @param x_m_name vector of names of the mediator variables
#' @param mod_0 model $\hat{m}_0$
#' @param mod_1 model $\hat{m}_1$
#' @param pred_mod_0 prediction function for model $\hat{m}_0(\cdot)$
#' @param pred_mod_1 prediction function for model $\hat{m}_1(\cdot)$
#' @param return_x if `TRUE` (default) the mediator variables are returned in the table, as well as their transported values
sate <- function(x_0, x_t, x_m_names, mod_0, mod_1, pred_mod_0, pred_mod_1, return_x = TRUE){
  if(is.vector(x_0)){
    # Univariate case
    new_data <- tibble(!!x_m_names := x_0)
  }else{
    new_data <- x_0
  }
  if(is.vector(x_t)){
    # Univariate case
    new_data_t <- tibble(!!x_m_names := x_t)
  }else{
    new_data_t <- x_t
  }
  
  
  # $\hat{m}_0(x_0)$
  y_0 <- pred_mod_0(object = mod_0, newdata = new_data)
  # $\hat{m}_1(x_0)$
  y_1 <- pred_mod_1(object = mod_1, newdata = new_data)
  # $\hat{m}_1(\mathcal{T}(x_0))$
  y_1_t <- pred_mod_1(object = mod_1, newdata = new_data_t)
  
  
  scate_tab <- 
    tibble(y_0 = y_0, y_1 = y_1, y_1_t = y_1_t,
           CATE = y_1-y_0, SCATE = y_1_t-y_0)
  
  if(return_x){
    new_data_t <- new_data_t %>% rename_all(~paste0(.x, "_t"))
    scate_tab <- bind_cols(new_data, new_data_t, scate_tab)
  }
  scate_tab
}