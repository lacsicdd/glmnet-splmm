############################################################
# ESTUDOS DE SIMULACAO MONTE CARLO
# 12 CONFIGURACOES
#
# ALGORITMO PROPOSTO: glmnet-splmm
#
# M = 30 replicacoes Monte Carlo
#
# A: n = 30, ni = 5, N = 150
# B: n = 60, ni = 10, N = 600
#
# beta*  = (1, 1, 0, ..., 0)
# beta** = (0.1, 0.1, 0, ..., 0)
#
# X*    = Normal multivariada, correlacao 0.90
# X**   = covariaveis mistas
# X***  = covariaveis mistas permutadas
#
# Modelo:
# y_ij = beta_0 + X_ij beta +
#        b_0i + b_1i tempo_ij + epsilon_ij
#
# INTEGRACAO ALGORITMICA glmnet-splmm (Secao 3.3)
#
#  ETAPA 1 (glmnet)
#    O glmnet e usado EXCLUSIVAMENTE para definir,
#    de forma automatica, o valor de referencia lambda_max
#    e a sequencia decrescente de valores candidatos do
#    parametro de regularizacao dos efeitos fixos (lam1).
#    O glmnet ignora a estrutura de correlacao do modelo
#    misto (nao usa efeitos aleatorios).
#
#  ETAPA 2 (splmmTuning)
#    A sequencia do glmnet e fornecida ao splmmTuning, que
#    ajusta o modelo misto penalizado para cada lambda e
#    seleciona o lambda de menor BIC.
#
#  ETAPA 3 (splmm)
#    Ajuste final do modelo misto penalizado com o lambda
#    escolhido pelo BIC. A dependencia dos dados e considerada
#    somente nas etapas 2 e 3.
#
# Intercepto NAO entra nas metricas.
############################################################

rm(list = ls())
gc()

# Execute com a RAIZ DO REPOSITORIO como diretorio de trabalho
# (abra o projeto no RStudio ou use setwd() para a raiz).
# Os arquivos de saida sao gravados na pasta results/.
dir.create("results", showWarnings = FALSE)

############################################################
# 1. PACOTES
############################################################

if (!requireNamespace("splmm", quietly = TRUE)) {
  install.packages("splmm")
}

library(splmm)

if (!requireNamespace("glmnet", quietly = TRUE)) {
  install.packages("glmnet")
}

library(glmnet)

############################################################
# 2. PARAMETROS GERAIS
############################################################

set.seed(2022)

M <- 30
p <- 9

############################################################
# IMPORTANTE:
# Substitua D e sigma2 pelos valores EXATOS da imagem
# da sua metodologia, caso sejam diferentes.
############################################################

D <- matrix(
  c(1.0, 0.25,
    0.25, 1.0),
  nrow = 2,
  byrow = TRUE
)

sigma2 <- 1

############################################################
# 3. VETORES DE BETAS
############################################################

beta_forte <- c(
  1, 1, 0, 0, 0, 0, 0, 0, 0
)

beta_fraco <- c(
  0.1, 0.1, 0, 0, 0, 0, 0, 0, 0
)

############################################################
# 4. FUNCAO PARA GERAR NORMAL MULTIVARIADA
############################################################

rmvnorm_base <- function(n, mu, Sigma) {
  
  p <- length(mu)
  
  Z <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )
  
  R <- chol(Sigma)
  
  X <- Z %*% R
  
  X <- sweep(
    X,
    2,
    mu,
    "+"
  )
  
  return(X)
}

############################################################
# 5. GERACAO DAS COVARIAVEIS X*
############################################################

gerar_X_star <- function(N) {
  
  Sigma <- matrix(
    0.90,
    nrow = p,
    ncol = p
  )
  
  diag(Sigma) <- 1
  
  X <- rmvnorm_base(
    n = N,
    mu = rep(0, p),
    Sigma = Sigma
  )
  
  colnames(X) <- paste0("X", 1:p)
  
  return(X)
}

############################################################
# 6. GERACAO DAS COVARIAVEIS X**
############################################################

gerar_X_star2 <- function(N) {
  
  X <- matrix(
    NA_real_,
    nrow = N,
    ncol = p
  )
  
  # X1 e X2 ~ Bernoulli(0.5)
  X[, 1] <- rbinom(N, 1, 0.5)
  X[, 2] <- rbinom(N, 1, 0.5)
  
  # X3 e X4 ~ Normal com correlacao 0.9
  Sigma_34 <- matrix(
    c(1, 0.9,
      0.9, 1),
    nrow = 2
  )
  
  X34 <- rmvnorm_base(
    N,
    c(0, 0),
    Sigma_34
  )
  
  X[, 3] <- X34[, 1]
  X[, 4] <- X34[, 2]
  
  # X5 e X6 ~ Normal com correlacao -0.5
  Sigma_56 <- matrix(
    c(1, -0.5,
      -0.5, 1),
    nrow = 2
  )
  
  X56 <- rmvnorm_base(
    N,
    c(0, 0),
    Sigma_56
  )
  
  X[, 5] <- X56[, 1]
  X[, 6] <- X56[, 2]
  
  # X7 e X8 ~ Beta(2,10)
  X[, 7] <- rbeta(N, 2, 10)
  X[, 8] <- rbeta(N, 2, 10)
  
  # X9 ~ Normal(6,1)
  X[, 9] <- rnorm(N, mean = 6, sd = 1)
  
  colnames(X) <- paste0("X", 1:p)
  
  return(X)
}

############################################################
# 7. GERACAO DAS COVARIAVEIS X***
############################################################

gerar_X_star3 <- function(N) {
  
  X <- matrix(
    NA_real_,
    nrow = N,
    ncol = p
  )
  
  # X1 e X2 ~ Normal com correlacao 0.9
  Sigma_12 <- matrix(
    c(1, 0.9,
      0.9, 1),
    nrow = 2
  )
  
  X12 <- rmvnorm_base(
    N,
    c(0, 0),
    Sigma_12
  )
  
  X[, 1] <- X12[, 1]
  X[, 2] <- X12[, 2]
  
  # X3 e X4 ~ Bernoulli(0.5)
  X[, 3] <- rbinom(N, 1, 0.5)
  X[, 4] <- rbinom(N, 1, 0.5)
  
  # X5 e X6 ~ Normal com correlacao -0.5
  Sigma_56 <- matrix(
    c(1, -0.5,
      -0.5, 1),
    nrow = 2
  )
  
  X56 <- rmvnorm_base(
    N,
    c(0, 0),
    Sigma_56
  )
  
  X[, 5] <- X56[, 1]
  X[, 6] <- X56[, 2]
  
  # X7 e X8 ~ Beta(2,10)
  X[, 7] <- rbeta(N, 2, 10)
  X[, 8] <- rbeta(N, 2, 10)
  
  # X9 ~ Normal(6,1)
  X[, 9] <- rnorm(N, mean = 6, sd = 1)
  
  colnames(X) <- paste0("X", 1:p)
  
  return(X)
}

############################################################
# 8. GERACAO DOS EFEITOS ALEATORIOS
############################################################

gerar_efeitos_aleatorios <- function(n, D) {
  
  b <- rmvnorm_base(
    n = n,
    mu = c(0, 0),
    Sigma = D
  )
  
  colnames(b) <- c(
    "b0",
    "b1"
  )
  
  return(b)
}

############################################################
# 9. GERACAO DE UM BANCO COMPLETO
############################################################

gerar_dados <- function(
    n,
    ni,
    beta,
    tipo_X,
    D,
    sigma2) {
  
  N <- n * ni
  
  ##########################################################
  # ID
  ##########################################################
  
  id <- rep(
    seq_len(n),
    each = ni
  )
  
  ##########################################################
  # TEMPO
  ##########################################################
  
  tempo <- rep(
    seq_len(ni) - 1,
    times = n
  )
  
  ##########################################################
  # MATRIZ X
  ##########################################################
  
  if (tipo_X == "X*") {
    
    X <- gerar_X_star(N)
    
  } else if (tipo_X == "X**") {
    
    X <- gerar_X_star2(N)
    
  } else if (tipo_X == "X***") {
    
    X <- gerar_X_star3(N)
    
  } else {
    
    stop("Tipo de X desconhecido.")
    
  }
  
  ##########################################################
  # EFEITOS ALEATORIOS
  ##########################################################
  
  b <- gerar_efeitos_aleatorios(
    n = n,
    D = D
  )
  
  ##########################################################
  # EFEITO ALEATORIO DE CADA OBSERVACAO
  ##########################################################
  
  b0 <- b[id, 1]
  b1 <- b[id, 2]
  
  ##########################################################
  # ERRO RESIDUAL
  ##########################################################
  
  epsilon <- rnorm(
    N,
    mean = 0,
    sd = sqrt(sigma2)
  )
  
  ##########################################################
  # RESPOSTA
  ##########################################################
  
  y <- as.numeric(
    X %*% beta +
      b0 +
      b1 * tempo +
      epsilon
  )
  
  ##########################################################
  # MATRIZ X DO SPLMM
  #
  # PRIMEIRA COLUNA = INTERCEPTO
  ##########################################################
  
  X_fit <- cbind(
    Intercepto = rep(1, N),
    X
  )
  
  X_fit <- as.matrix(X_fit)
  
  ##########################################################
  # MATRIZ Z DOS EFEITOS ALEATORIOS
  #
  # INTERCEPTO + TEMPO
  ##########################################################
  
  Z_fit <- cbind(
    Intercepto = rep(1, N),
    Tempo = tempo
  )
  
  Z_fit <- as.matrix(Z_fit)
  
  ##########################################################
  # GARANTIAS CONTRA O ERRO:
  # "x is not a data.frame or matrix"
  ##########################################################
  
  storage.mode(X_fit) <- "double"
  storage.mode(Z_fit) <- "double"
  storage.mode(y) <- "double"
  
  ##########################################################
  # RETORNO
  ##########################################################
  
  return(
    list(
      y = y,
      X = X,
      X_fit = X_fit,
      Z_fit = Z_fit,
      id = id,
      tempo = tempo,
      beta = beta
    )
  )
}

############################################################
# 10. CALCULO DAS METRICAS
############################################################

calcular_metricas <- function(
    beta_hat,
    beta_true,
    tol = 1e-6) {
  
  ##########################################################
  # beta_hat deve ter SOMENTE os 9 efeitos de X.
  #
  # O intercepto esta em beta_hat[1] e NAO entra.
  ##########################################################
  
  if (length(beta_hat) >= 10) {
    
    beta_hat <- beta_hat[2:10]
    
  }
  
  ##########################################################
  # Garantia de tamanho
  ##########################################################
  
  if (length(beta_hat) != 9) {
    
    stop(
      paste(
        "Numero incorreto de coeficientes:",
        length(beta_hat)
      )
    )
    
  }
  
  beta_hat <- as.numeric(beta_hat)
  beta_true <- as.numeric(beta_true)
  
  ##########################################################
  # Se houver NA, a replicacao nao deve produzir
  # NA silenciosamente.
  ##########################################################
  
  if (any(!is.finite(beta_hat))) {
    
    return(
      list(
        sensibilidade = NA_real_,
        especificidade = NA_real_,
        EQM_rep = NA_real_,
        RMSE_rep = NA_real_
      )
    )
    
  }
  
  ##########################################################
  # Variaveis verdadeiramente ativas
  ##########################################################
  
  verdadeiras <- abs(beta_true) > tol
  
  ##########################################################
  # Variaveis estimadas como ativas
  ##########################################################
  
  selecionadas <- abs(beta_hat) > tol
  
  ##########################################################
  # SENSIBILIDADE
  ##########################################################
  
  if (sum(verdadeiras) > 0) {
    
    sensibilidade <- mean(
      selecionadas[verdadeiras]
    )
    
  } else {
    
    sensibilidade <- 0
    
  }
  
  ##########################################################
  # ESPECIFICIDADE
  ##########################################################
  
  if (sum(!verdadeiras) > 0) {
    
    especificidade <- mean(
      !selecionadas[!verdadeiras]
    )
    
  } else {
    
    especificidade <- 0
    
  }
  
  ##########################################################
  # EQM / RMSE
  #
  # Conforme a expressao definida:
  #
  # EQM =
  # sqrt[
  #   (beta_hat-beta)'(beta_hat-beta) / M
  # ]
  #
  # ATENCAO:
  # nesta funcao calculamos o erro da replicacao.
  # A divisao por M sera feita posteriormente,
  # depois das 30 replicacoes.
  ##########################################################
  
  erro <- beta_hat - beta_true
  
  soma_quadrados <- sum(
    erro^2
  )
  
  ##########################################################
  # Retornamos o erro quadratico da replicacao.
  ##########################################################
  
  return(
    list(
      sensibilidade = sensibilidade,
      especificidade = especificidade,
      EQM_rep = soma_quadrados,
      RMSE_rep = sqrt(soma_quadrados)
    )
  )
}

############################################################
# 11. FUNCAO PARA EXTRAIR COEFICIENTES DO SPLMM
############################################################

extrair_beta <- function(fit) {
  
  ##########################################################
  # O splmm fornece os coeficientes em $coefficients
  ##########################################################
  
  beta_hat <- fit$coefficients
  
  beta_hat <- as.numeric(
    beta_hat
  )
  
  ##########################################################
  # Deve haver 10 coeficientes:
  #
  # 1 intercepto
  # 9 covariaveis
  ##########################################################
  
  if (length(beta_hat) != 10) {
    
    stop(
      paste(
        "O modelo retornou",
        length(beta_hat),
        "coeficientes; eram esperados 10."
      )
    )
    
  }
  
  return(beta_hat)
}

############################################################
# 11B. ETAPA 1 - SEQUENCIA DE LAMBDAS VIA GLMNET
#
# O glmnet (lasso, alpha = 1) e ajustado SOMENTE com as
# covariaveis (sem o intercepto, que o glmnet ja inclui
# por conta propria) e SEM estrutura de efeitos aleatorios.
#
# Do glmnet aproveitamos apenas:
#   - lambda_max : menor lambda que zera todos os
#                  coeficientes (max da sequencia)
#   - sequencia decrescente de lambdas candidatos
#
# nlambda e lambda.min.ratio controlam o tamanho e a
# amplitude da grade (ver secao 15).
############################################################

sequencia_lambda_glmnet <- function(
    X,
    y,
    nlambda,
    lambda_min_ratio) {
  
  fit_glmnet <- glmnet(
    x = X,
    y = y,
    family = "gaussian",
    alpha = 1,
    nlambda = nlambda,
    lambda.min.ratio = lambda_min_ratio,
    standardize = TRUE
  )
  
  ##########################################################
  # fit_glmnet$lambda ja e decrescente
  ##########################################################
  
  lambdas <- as.numeric(fit_glmnet$lambda)
  
  lambdas <- lambdas[
    is.finite(lambdas) & lambdas > 0
  ]
  
  lambdas <- sort(
    unique(lambdas),
    decreasing = TRUE
  )
  
  if (length(lambdas) < 2) {
    
    stop("glmnet retornou menos de 2 valores de lambda.")
    
  }
  
  return(
    list(
      lambda_max = max(lambdas),
      lambdas = lambdas
    )
  )
}

############################################################
# 12. UMA REPLICACAO (glmnet -> splmmTuning -> splmm)
############################################################

uma_replicacao <- function(
    n,
    ni,
    beta,
    tipo_X,
    D,
    sigma2,
    lambda2,
    nlambda,
    lambda_min_ratio) {
  
  dados <- gerar_dados(
    n = n,
    ni = ni,
    beta = beta,
    tipo_X = tipo_X,
    D = D,
    sigma2 = sigma2
  )
  
  ##########################################################
  # Checagens explicitas
  ##########################################################
  
  stopifnot(
    is.matrix(dados$X_fit),
    is.matrix(dados$Z_fit),
    is.numeric(dados$y),
    nrow(dados$X_fit) == length(dados$y),
    nrow(dados$Z_fit) == length(dados$y),
    ncol(dados$X_fit) == 10,
    ncol(dados$Z_fit) == 2
  )
  
  ##########################################################
  # ETAPA 1: GLMNET
  # Define lambda_max e a sequencia decrescente de lam1
  ##########################################################
  
  seq_glmnet <- sequencia_lambda_glmnet(
    X = dados$X,
    y = dados$y,
    nlambda = nlambda,
    lambda_min_ratio = lambda_min_ratio
  )
  
  lam1_seq <- seq_glmnet$lambdas
  
  ##########################################################
  # ETAPA 2: SPLMMTUNING
  # Ajusta o modelo misto penalizado para cada lambda da
  # sequencia do glmnet e calcula o BIC de cada um.
  ##########################################################
  
  tuning <- splmmTuning(
    x = dados$X_fit,
    y = dados$y,
    z = dados$Z_fit,
    grp = dados$id,
    lam1.seq = lam1_seq,
    lam2.seq = lambda2,
    nonpen.b = 1,
    nonpen.L = 1,
    penalty.b = "lasso",
    penalty.L = "lasso",
    CovOpt = "nlminb",
    standardize = TRUE,
    control = splmmControl(
      tol = 1e-4,
      trace = 1,
      maxIter = 30,
      maxArmijo = 10,
      number = 5
    )
  )
  
  bic_seq <- as.numeric(tuning$BIC.lam1)
  
  if (length(bic_seq) != length(lam1_seq)) {
    
    stop(
      paste(
        "BIC.lam1 tem",
        length(bic_seq),
        "valores; esperados",
        length(lam1_seq)
      )
    )
    
  }
  
  bic_seq[!is.finite(bic_seq)] <- NA_real_
  
  if (all(is.na(bic_seq))) {
    
    stop("Nenhum lambda da sequencia produziu BIC valido.")
    
  }
  
  ##########################################################
  # Lambda escolhido pelo criterio de BIC
  ##########################################################
  
  lambda1_bic <- lam1_seq[
    which.min(bic_seq)
  ]
  
  ##########################################################
  # ETAPA 3: AJUSTE FINAL SPLMM COM O LAMBDA ESCOLHIDO
  ##########################################################
  
  fit <- splmm(
    x = dados$X_fit,
    y = dados$y,
    z = dados$Z_fit,
    grp = dados$id,
    lam1 = lambda1_bic,
    lam2 = lambda2,
    nonpen.b = 1,
    nonpen.L = 1,
    penalty.b = "lasso",
    penalty.L = "lasso",
    CovOpt = "nlminb",
    standardize = TRUE,
    control = splmmControl(
      tol = 1e-4,
      trace = 1,
      maxIter = 30,
      maxArmijo = 10,
      number = 5
    )
  )
  
  ##########################################################
  # Coeficientes
  ##########################################################
  
  beta_hat <- extrair_beta(
    fit
  )
  
  ##########################################################
  # Metricas
  ##########################################################
  
  met <- calcular_metricas(
    beta_hat = beta_hat,
    beta_true = beta
  )
  
  ##########################################################
  # Retorno
  ##########################################################
  
  return(
    list(
      beta_hat = beta_hat,
      sensibilidade = met$sensibilidade,
      especificidade = met$especificidade,
      EQM_rep = met$EQM_rep,
      RMSE_rep = met$RMSE_rep,
      converged = fit$converged,
      bic = fit$bic,
      lambda_max = seq_glmnet$lambda_max,
      lambda1_bic = lambda1_bic
    )
  )
}

############################################################
# 13. AS 12 CONFIGURACOES
############################################################

configuracoes <- expand.grid(
  amostra = c("A", "B"),
  beta = c("forte", "fraco"),
  X = c("X*", "X**", "X***"),
  stringsAsFactors = FALSE
)

############################################################
# 14. PARAMETROS DOS CENARIOS
############################################################

configuracoes$n <- ifelse(
  configuracoes$amostra == "A",
  30,
  60
)

configuracoes$ni <- ifelse(
  configuracoes$amostra == "A",
  5,
  10
)

configuracoes$N <- configuracoes$n * configuracoes$ni

configuracoes$beta_vec <- lapply(
  configuracoes$beta,
  function(x) {
    if (x == "forte") {
      beta_forte
    } else {
      beta_fraco
    }
  }
)

############################################################
# 15. PARAMETROS DA REGULARIZACAO
#
# lam1 (efeitos fixos):
#   NAO e fixo. Vem da sequencia do glmnet e e escolhido
#   pelo BIC dentro do splmmTuning, em cada replicacao.
#
# nlambda:
#   numero de valores candidatos pedidos ao glmnet.
#
# lambda_min_ratio:
#   razao lambda_min / lambda_max da sequencia do glmnet.
#
# lam2 (covariancia dos efeitos aleatorios):
#   mantido fixo, como no estudo anterior.
############################################################

nlambda <- 20
lambda_min_ratio <- 0.01

lambda2 <- 0.10

############################################################
# 16. OBJETO PARA ARMAZENAR RESULTADOS
############################################################

resultados <- list()

indice_resultado <- 1

############################################################
# 17. LOOP DAS 12 CONFIGURACOES
############################################################

inicio_total <- Sys.time()

for (c in seq_len(nrow(configuracoes))) {
  
  cat("\n")
  cat("====================================================\n")
  cat("CONFIGURACAO:", c, "de", nrow(configuracoes), "\n")
  cat("Amostra:", configuracoes$amostra[c], "\n")
  cat("n =", configuracoes$n[c],
      " ni =", configuracoes$ni[c],
      " N =", configuracoes$N[c], "\n")
  cat("Beta:", configuracoes$beta[c], "\n")
  cat("X:", configuracoes$X[c], "\n")
  cat("====================================================\n")
  
  n_atual <- configuracoes$n[c]
  ni_atual <- configuracoes$ni[c]
  
  beta_atual <- configuracoes$beta_vec[[c]]
  
  X_atual <- configuracoes$X[c]
  
  ########################################################
  # MATRIZ PARA AS 30 REPLICACOES
  ########################################################
  
  res_cenario <- data.frame(
    replicacao = seq_len(M),
    sensibilidade = NA_real_,
    especificidade = NA_real_,
    EQM_rep = NA_real_,
    RMSE_rep = NA_real_,
    converged = NA_real_,
    BIC = NA_real_,
    lambda_max = NA_real_,
    lambda1_bic = NA_real_,
    erro = NA_character_
  )
  
  ########################################################
  # LOOP MONTE CARLO
  ########################################################
  
  for (m in seq_len(M)) {
    
    cat(
      "  Replicacao",
      m,
      "de",
      M,
      "\n"
    )
    
    resultado_m <- tryCatch(
      
      {
        
        uma_replicacao(
          n = n_atual,
          ni = ni_atual,
          beta = beta_atual,
          tipo_X = X_atual,
          D = D,
          sigma2 = sigma2,
          lambda2 = lambda2,
          nlambda = nlambda,
          lambda_min_ratio = lambda_min_ratio
        )
        
      },
      
      error = function(e) {
        
        list(
          beta_hat = rep(NA_real_, 10),
          sensibilidade = NA_real_,
          especificidade = NA_real_,
          EQM_rep = NA_real_,
          RMSE_rep = NA_real_,
          converged = NA_real_,
          bic = NA_real_,
          lambda_max = NA_real_,
          lambda1_bic = NA_real_,
          erro = conditionMessage(e)
        )
        
      }
      
    )
    
    ######################################################
    # Armazenamento
    ######################################################
    
    res_cenario$sensibilidade[m] <-
      resultado_m$sensibilidade
    
    res_cenario$especificidade[m] <-
      resultado_m$especificidade
    
    res_cenario$EQM_rep[m] <-
      resultado_m$EQM_rep
    
    res_cenario$RMSE_rep[m] <-
      resultado_m$RMSE_rep
    
    res_cenario$converged[m] <-
      resultado_m$converged
    
    res_cenario$BIC[m] <-
      resultado_m$bic
    
    res_cenario$lambda_max[m] <-
      resultado_m$lambda_max
    
    res_cenario$lambda1_bic[m] <-
      resultado_m$lambda1_bic
    
    if (!is.null(resultado_m$erro)) {
      
      res_cenario$erro[m] <-
        resultado_m$erro
      
    }
    
  }
  
  ########################################################
  # RESULTADO RESUMIDO DO CENARIO
  ########################################################
  
  n_validas <- sum(
    is.finite(res_cenario$EQM_rep)
  )
  
  ########################################################
  # IMPORTANTE:
  #
  # NAO usamos simplesmente mean(..., na.rm=TRUE)
  # sem informar quantas replicacoes foram validas.
  ########################################################
  
  if (n_validas > 0) {
    
    sens_media <- mean(
      res_cenario$sensibilidade,
      na.rm = TRUE
    )
    
    esp_media <- mean(
      res_cenario$especificidade,
      na.rm = TRUE
    )
    
    ######################################################
    # EQM FINAL
    #
    # sqrt[
    #   sum_m ||beta_hat_m-beta||^2 / M
    # ]
    #
    # Se todas as M=30 replicacoes forem validas:
    # divisor = 30.
    ######################################################
    
    EQM_final <- sqrt(
      sum(
        res_cenario$EQM_rep,
        na.rm = TRUE
      ) / n_validas
    )
    
    ######################################################
    # RMSE medio por replicacao
    ######################################################
    
    RMSE_medio <- mean(
      res_cenario$RMSE_rep,
      na.rm = TRUE
    )
    
  } else {
    
    sens_media <- NA_real_
    esp_media <- NA_real_
    EQM_final <- NA_real_
    RMSE_medio <- NA_real_
    
  }
  
  ########################################################
  # SALVA RESULTADO
  ########################################################
  
  resultados[[indice_resultado]] <- data.frame(
    
    configuracao = c,
    
    amostra = configuracoes$amostra[c],
    
    n = n_atual,
    
    ni = ni_atual,
    
    N = configuracoes$N[c],
    
    beta = configuracoes$beta[c],
    
    X = X_atual,
    
    M = M,
    
    replicacoes_validas = n_validas,
    
    replicacoes_com_erro = M - n_validas,
    
    sensibilidade = sens_media,
    
    especificidade = esp_media,
    
    EQM = EQM_final,
    
    RMSE = RMSE_medio,
    
    lambda_max_medio = mean(
      res_cenario$lambda_max,
      na.rm = TRUE
    ),
    
    lambda1_bic_medio = mean(
      res_cenario$lambda1_bic,
      na.rm = TRUE
    )
    
  )
  
  indice_resultado <- indice_resultado + 1
  
  ########################################################
  # SALVA TAMBEM AS REPLICACOES
  ########################################################
  
  resultados[[paste0(
    "replicacoes_",
    c
  )]] <- res_cenario
  
}

############################################################
# 18. JUNTA OS RESULTADOS DOS 12 CENARIOS
############################################################

resultados_finais <- do.call(
  rbind,
  resultados[
    sapply(
      resultados,
      function(x) {
        is.data.frame(x) &&
          "configuracao" %in% names(x)
      }
    )
  ]
)

rownames(resultados_finais) <- NULL

############################################################
# 19. MOSTRA RESULTADOS
############################################################

cat("\n\n")
cat("====================================================\n")
cat("RESULTADOS FINAIS\n")
cat("====================================================\n")

print(
  resultados_finais
)

############################################################
# 20. SALVA RESULTADOS
############################################################

write.csv(
  resultados_finais,
  "results/resultados_12_cenarios_M30_glmnet_splmm.csv",
  row.names = FALSE
)

############################################################
# 21. SALVA RESULTADOS DAS REPLICACOES
############################################################

replicacoes_finais <- do.call(
  rbind,
  resultados[
    sapply(
      resultados,
      function(x) {
        is.data.frame(x) &&
          "replicacao" %in% names(x)
      }
    )
  ]
)

write.csv(
  replicacoes_finais,
  "results/resultados_replicacoes_M30_glmnet_splmm.csv",
  row.names = FALSE
)

############################################################
# 22. TEMPO TOTAL
############################################################

fim_total <- Sys.time()

cat("\n")
cat("====================================================\n")
cat("TEMPO TOTAL\n")
cat("====================================================\n")

print(
  fim_total - inicio_total
)

cat("\nArquivos gerados:\n")
cat(" - results/resultados_12_cenarios_M30_glmnet_splmm.csv\n")
cat(" - results/resultados_replicacoes_M30_glmnet_splmm.csv\n")
