############################################################
# APLICACAO: DADOS DE RIBOFLAVINA (Buhlmann et al., 2014)
#
# Selecao de variaveis com o algoritmo glmnet-splmm
#
# Modelo:
#   y = X*beta + Zb + epsilon
#
#   y     = log da taxa de producao de riboflavina
#   X     = tempo + expressao (log) de 100 genes
#   28 cepas = grupos
#   medidas repetidas da cepa (71 observacoes)
#
# Efeitos aleatorios: intercepto + tempo (Z = [1, tempo])
#
# ALGORITMO glmnet-splmm (tres etapas)
#
#  ETAPA 1 - glmnet
#    Usado EXCLUSIVAMENTE para definir, de forma
#    automatica, o valor de referencia lambda_max e a
#    sequencia decrescente de valores candidatos do
#    parametro de regularizacao dos efeitos fixos (lam1).
#    O glmnet NAO considera a correlacao entre medidas
#    repetidas (nao tem efeitos aleatorios).
#
#  ETAPA 2 - splmmTuning
#    Recebe a sequencia do glmnet, ajusta o modelo misto
#    penalizado para cada lambda e calcula o BIC.
#    Escolhe-se o lambda de MENOR BIC.
#
#  ETAPA 3 - splmm
#    Ajuste final do modelo misto penalizado com o lambda
#    escolhido. A estrutura de dependencia dos dados entra
#    nas etapas 2 e 3.
#
# Variavel selecionada: |beta_hat| > tol
# O intercepto NAO entra na selecao.
#
# Arquivos necessarios (pasta data/ do repositorio):
#   - data/riboflavinv100.csv
#   - data/riboflavinv100_structure.txt
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
# 2. PARAMETROS DO ALGORITMO
#
# Mesmos valores da simulacao:
#   nlambda          : numero de lambdas candidatos
#                      pedidos ao glmnet
#   lambda_min_ratio : razao lambda_min / lambda_max
#                      da sequencia do glmnet
#   lambda2          : penalidade da matriz de covariancia
#                      dos efeitos aleatorios (fixa)
#   tol              : |beta_hat| <= tol => variavel
#                      NAO selecionada
############################################################

nlambda <- 20
lambda_min_ratio <- 0.01

lambda2 <- 0.10

tol <- 1e-6

# Maximo de iteracoes do splmm nesta aplicacao
# (na simulacao: 30; aqui os dados reais exigem mais)
maxIter_aplic <- 100

############################################################
# 3. LEITURA DOS DADOS
#
# Estrutura do data/riboflavinv100.csv:
#   - Coluna 1     : nome da linha
#   - Linha 1      : resposta y (q_RIBFLV)
#   - Linhas 2-101 : 100 genes
#   - Colunas 2-72 : 71 observacoes (amostras)
#
# Os genes estao nas LINHAS e as amostras nas COLUNAS.
# Por isso a matriz de genes e transposta: queremos
# 71 linhas (observacoes) x 100 colunas (genes).
############################################################

dados_brutos <- read.csv(
  "data/riboflavinv100.csv",
  header = TRUE,
  check.names = FALSE
)

cat("Dimensao do arquivo:", dim(dados_brutos), "\n")

# Resposta: primeira linha, sem a primeira coluna
y_original <- as.numeric(
  unlist(dados_brutos[1, -1])
)

# Nomes dos genes: primeira coluna, a partir da linha 2
nomes_genes <- as.character(
  dados_brutos[2:nrow(dados_brutos), 1]
)

# "YQGN-P_i_at" tem hifen, que atrapalha em alguns
# comandos do R. Trocamos por "_".
nomes_genes <- gsub("-", "_", nomes_genes)

# Matriz de genes transposta (observacoes x genes)
genes <- t(
  as.matrix(
    dados_brutos[2:nrow(dados_brutos), -1]
  )
)

colnames(genes) <- nomes_genes
rownames(genes) <- NULL

storage.mode(genes) <- "double"

cat("Numero de observacoes:", nrow(genes), "\n")
cat("Numero de genes      :", ncol(genes), "\n")

############################################################
# 4. ESTRUTURA DE AGRUPAMENTO (CEPAS E TEMPO)
#
# data/riboflavinv100_structure.txt, na mesma ordem das colunas
# do CSV:
#   id    : nome da cepa
#   time  : tempo da medida
#   newid : numero do grupo (1 a 28)
############################################################

estrutura <- read.table(
  "data/riboflavinv100_structure.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(estrutura) == length(y_original),
  nrow(estrutura) == nrow(genes)
)

cat("Numero de cepas (grupos):",
    length(unique(estrutura$newid)), "\n")

############################################################
# 5. ORDENA AS OBSERVACOES POR CEPA
#
# No arquivo as medidas de uma mesma cepa nao sao
# consecutivas (ex.: 1374_ aparece nas linhas 25, 28 e 29).
# Na simulacao os dados chegam ao splmm com as observacoes
# de cada grupo juntas; aqui reordenamos as linhas para
# reproduzir essa organizacao.
#
# Isso NAO altera o modelo: y, X, Z e grupo sao
# reordenados juntos.
############################################################

ordem <- order(estrutura$newid)

y_original <- y_original[ordem]
genes      <- genes[ordem, , drop = FALSE]
estrutura  <- estrutura[ordem, ]

grupo <- as.integer(estrutura$newid)

############################################################
# 6. PADRONIZACAO
#
# Resposta, tempo e genes sao centrados e escalados
# (media 0, desvio-padrao 1), para que todos os
# coeficientes fiquem na mesma escala e a penalidade
# trate as covariaveis igualmente.
############################################################

y     <- as.numeric(scale(y_original))
tempo <- as.numeric(scale(estrutura$time))
genes <- scale(genes)

############################################################
# 7. MATRIZES DO MODELO
#
# Mesmo formato da simulacao:
#   X     : covariaveis SEM intercepto (tempo + 100 genes);
#           e o que o glmnet recebe (ele ja inclui
#           intercepto por conta propria)
#   X_fit : intercepto + covariaveis; e o que o splmm e o
#           splmmTuning recebem
#   Z_fit : intercepto + tempo (efeitos aleatorios)
############################################################

X <- cbind(
  tempo = tempo,
  genes
)

X <- as.matrix(X)

X_fit <- cbind(
  Intercepto = rep(1, length(y)),
  X
)

Z_fit <- cbind(
  Intercepto = rep(1, length(y)),
  Tempo = tempo
)

X_fit <- as.matrix(X_fit)
Z_fit <- as.matrix(Z_fit)

storage.mode(X)     <- "double"
storage.mode(X_fit) <- "double"
storage.mode(Z_fit) <- "double"
storage.mode(y)     <- "double"

p <- ncol(X)   # numero de covariaveis (sem intercepto)

############################################################
# 8. CHECAGENS EXPLICITAS (como na simulacao)
############################################################

stopifnot(
  is.matrix(X),
  is.matrix(X_fit),
  is.matrix(Z_fit),
  is.numeric(y),
  nrow(X_fit) == length(y),
  nrow(Z_fit) == length(y),
  length(grupo) == length(y),
  ncol(X_fit) == p + 1,
  ncol(Z_fit) == 2,
  !anyNA(X_fit), !anyNA(Z_fit), !anyNA(y)
)

cat("\nResumo dos dados para o ajuste:\n")
cat("  Observacoes (N)      :", length(y), "\n")
cat("  Cepas (n)            :", length(unique(grupo)), "\n")
cat("  Covariaveis (p)      :", p, "(tempo + genes)\n")
cat("  Colunas de X_fit     :", ncol(X_fit), "(com intercepto)\n")
cat("  Colunas de Z_fit     :", ncol(Z_fit), "\n")

############################################################
# 9. FUNCAO PARA EXTRAIR OS COEFICIENTES DO SPLMM
#
# Igual a da simulacao: os coeficientes estao em
# fit$coefficients (intercepto + covariaveis).
############################################################

extrair_beta <- function(fit, n_esperado) {
  
  beta_hat <- as.numeric(fit$coefficients)
  
  # (se o objeto nao tiver $coefficients, usa $fixef,
  #  como no Aplication.R original)
  if (length(beta_hat) == 0) {
    
    beta_hat <- as.numeric(fit$fixef)
    
  }
  
  if (length(beta_hat) != n_esperado) {
    
    stop(
      paste(
        "O modelo retornou",
        length(beta_hat),
        "coeficientes; eram esperados",
        n_esperado
      )
    )
    
  }
  
  return(beta_hat)
}

############################################################
# 10. ETAPA 1 - SEQUENCIA DE LAMBDAS VIA GLMNET
#
# O glmnet (lasso, alpha = 1) e ajustado somente com as
# covariaveis e SEM efeitos aleatorios.
#
# Dele aproveitamos apenas:
#   - lambda_max : menor lambda que zera todos os
#                  coeficientes (maximo da sequencia)
#   - sequencia decrescente de lambdas candidatos
############################################################

etapa1_glmnet <- function(
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
  
  # fit_glmnet$lambda ja e decrescente
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

inicio <- Sys.time()

seq_glmnet <- etapa1_glmnet(
  X = X,
  y = y,
  nlambda = nlambda,
  lambda_min_ratio = lambda_min_ratio
)

lam1_glmnet <- seq_glmnet$lambdas

cat("\n--- ETAPA 1 (glmnet) ---\n")
cat("lambda_max                :", seq_glmnet$lambda_max, "\n")
cat("Numero de lambdas na grade:", length(lam1_glmnet), "\n")
cat("Menor lambda da grade     :", min(lam1_glmnet), "\n")

############################################################
# 11. ETAPA 2 - SPLMMTUNING (BIC PARA CADA LAMBDA DA GRADE)
#
# O splmmTuning e chamado para cada lambda da sequencia do
# glmnet. Com um unico lambda ele roda o splmm e devolve o
# ajuste, do qual lemos:
#   - o BIC                        (fit$bic)
#   - o numero de coef. ativos     (|beta| > tol)
#   - o indicador de convergencia  (fit$converged)
#
# Um ajuste e VALIDO quando:
#   (i)  o BIC e finito; e
#   (ii) o numero de coeficientes ativos e menor que
#        min(p, N): quando chega nesse limite o splmm
#        interrompe o ajuste (modelo saturado).
#
# Escala da grade:
#   1 : grade do glmnet como esta (padrao do algoritmo)
#   N : grade do glmnet x N. O glmnet usa a perda dividida
#       por N; o splmm parece usar a perda sem essa divisao,
#       entao N * lambda_glmnet e o valor equivalente.
#       (So e usada se a escala 1 nao gerar nenhum ajuste
#       valido.)
############################################################

controle <- splmmControl(
  tol = 1e-4,
  trace = 0,
  maxIter = maxIter_aplic,
  maxArmijo = 10,
  number = 5
)

n_max_ativas <- min(p, length(y))

ajuste_um_lambda <- function(lam1) {
  
  saida <- tryCatch(
    {
      
      tun <- splmmTuning(
        x = X_fit,
        y = y,
        z = Z_fit,
        grp = grupo,
        lam1.seq = lam1,
        lam2.seq = lambda2,
        nonpen.b = 1,
        nonpen.L = 1,
        penalty.b = "lasso",
        penalty.L = "lasso",
        CovOpt = "nlminb",
        standardize = TRUE,
        control = controle
      )
      
      # BIC: fit$bic (como na simulacao); alternativas caso
      # o objeto tenha outro nome
      bic <- tun[["bic"]]
      
      if (is.null(bic)) bic <- tun[["bicc"]]
      
      if (is.null(bic) && !is.null(tun[["BIC.lam1"]])) {
        bic <- min(tun[["BIC.lam1"]], na.rm = TRUE)
      }
      
      bic <- suppressWarnings(as.numeric(bic)[1])
      
      if (length(bic) == 0 || !is.finite(bic)) bic <- NA_real_
      
      beta_try <- try(
        extrair_beta(tun, p + 1),
        silent = TRUE
      )
      
      n_ativas <- if (inherits(beta_try, "try-error")) {
        NA_integer_
      } else {
        sum(abs(beta_try[-1]) > tol)
      }
      
      conv <- tun[["converged"]]
      conv <- if (is.null(conv)) NA_character_ else as.character(conv)[1]
      
      list(
        bic = bic,
        n_ativas = n_ativas,
        convergiu = conv,
        erro = NA_character_
      )
      
    },
    error = function(e) {
      list(
        bic = NA_real_,
        n_ativas = NA_integer_,
        convergiu = NA_character_,
        erro = conditionMessage(e)
      )
    }
  )
  
  return(saida)
}

cat("\n--- ETAPA 2 (splmmTuning) ---\n")
cat("Limite de coeficientes ativos (min(p, N)):",
    n_max_ativas, "\n")

escalas <- c(1, length(y))

encontrou <- FALSE

for (esc in escalas) {
  
  lam1_seq <- lam1_glmnet * esc
  
  cat("\nEscala da grade =", esc,
      "| lambda1 de", format(max(lam1_seq), digits = 4),
      "a", format(min(lam1_seq), digits = 4), "\n")
  
  bic_bruto <- rep(NA_real_, length(lam1_seq))
  n_ativas  <- rep(NA_integer_, length(lam1_seq))
  convergiu <- rep(NA_character_, length(lam1_seq))
  erros     <- rep(NA_character_, length(lam1_seq))
  
  for (k in seq_along(lam1_seq)) {
    
    r <- ajuste_um_lambda(lam1_seq[k])
    
    bic_bruto[k] <- r$bic
    n_ativas[k]  <- r$n_ativas
    convergiu[k] <- r$convergiu
    erros[k]     <- r$erro
    
  }
  
  valido <- is.finite(bic_bruto) &
    (is.na(n_ativas) | n_ativas < n_max_ativas)
  
  tabela_bic <- data.frame(
    lambda1 = lam1_seq,
    BIC = bic_bruto,
    n_ativas = n_ativas,
    convergiu = convergiu,
    valido = valido
  )
  
  print(tabela_bic, digits = 6)
  
  cat("Ajustes validos:", sum(valido),
      "de", length(lam1_seq), "\n")
  
  if (any(valido)) {
    
    encontrou    <- TRUE
    escala_usada <- esc
    break
    
  } else if (any(!is.na(erros))) {
    
    cat("Primeiro erro registrado:",
        na.omit(erros)[1], "\n")
    
  }
  
}

if (!encontrou) {
  
  stop(
    paste(
      "Nenhum lambda produziu ajuste valido (BIC finito e",
      "menos de", n_max_ativas, "coeficientes ativos).",
      "Veja as tabelas acima."
    )
  )
  
}

# BIC usado na escolha: so entre os ajustes validos
bic_seq <- ifelse(valido, bic_bruto, NA_real_)

# Lambda escolhido pelo criterio de BIC
posicao_bic <- which.min(bic_seq)
lambda1_bic <- lam1_seq[posicao_bic]

cat("\nEscala usada:", escala_usada, "\n")
cat("Lambda escolhido pelo BIC:", lambda1_bic,
    "(posicao", posicao_bic, "de", length(lam1_seq), ")\n")

# Avisos sobre a posicao do minimo
if (posicao_bic == which(valido)[1]) {
  
  cat("ATENCAO: o menor BIC esta no lambda valido mais alto.",
      "Considere ajustar nlambda / lambda_min_ratio.\n")
  
}

if (posicao_bic == max(which(valido))) {
  
  cat("ATENCAO: o menor BIC esta no lambda valido mais baixo",
      "(proximo da saturacao do modelo).\n")
  
}

if (!isTRUE(as.logical(convergiu[posicao_bic]))) {
  
  cat("ATENCAO: indicador de convergencia do lambda escolhido:",
      convergiu[posicao_bic], "\n")
  
}

############################################################
# 12. ETAPA 3 - AJUSTE FINAL SPLMM
#
# Chamada identica a da simulacao (exceto maxIter, ver
# secao 2), com lam1 = lambda escolhido pelo BIC.
############################################################

cat("\n--- ETAPA 3 (splmm final) ---\n")

fit <- splmm(
  x = X_fit,
  y = y,
  z = Z_fit,
  grp = grupo,
  lam1 = lambda1_bic,
  lam2 = lambda2,
  nonpen.b = 1,
  nonpen.L = 1,
  penalty.b = "lasso",
  penalty.L = "lasso",
  CovOpt = "nlminb",
  standardize = TRUE,
  control = controle
)

fim <- Sys.time()

############################################################
# 13. COEFICIENTES ESTIMADOS E VARIAVEIS SELECIONADAS
#
# Regra da simulacao: variavel "selecionada" quando
# |beta_hat| > tol. O intercepto (posicao 1) fica de fora.
############################################################

beta_hat <- extrair_beta(
  fit,
  n_esperado = p + 1
)

intercepto <- beta_hat[1]
beta_covar <- beta_hat[2:(p + 1)]

selecionada <- abs(beta_covar) > tol

resultado <- data.frame(
  covariavel  = colnames(X),
  beta_hat    = beta_covar,
  selecionada = selecionada,
  stringsAsFactors = FALSE
)

# Somente as selecionadas, da maior para a menor |beta|
selecionadas <- resultado[resultado$selecionada, ]
selecionadas <- selecionadas[
  order(abs(selecionadas$beta_hat), decreasing = TRUE),
]
rownames(selecionadas) <- NULL

############################################################
# 14. RESULTADOS
############################################################

cat("\n====================================================\n")
cat("RESULTADOS DO ALGORITMO glmnet-splmm\n")
cat("====================================================\n")

cat("lambda_max (glmnet)      :", seq_glmnet$lambda_max, "\n")
cat("Escala da grade         :", escala_usada, "\n")
cat("lambda1 escolhido (BIC)  :", lambda1_bic, "\n")
cat("lambda2 (fixo)           :", lambda2, "\n")
cat("Convergiu (ajuste final) :", fit$converged, "\n")
cat("BIC (ajuste final)       :", fit$bic, "\n")
cat("Intercepto estimado      :", round(intercepto, 4), "\n")
cat("Variaveis selecionadas   :", sum(selecionada),
    "de", p, "\n")
cat("Variaveis zeradas        :", sum(!selecionada),
    "de", p, "\n")

cat("\nVariaveis selecionadas (ordem decrescente de |beta|):\n")
print(selecionadas, digits = 4)

cat("\nTempo de execucao:\n")
print(fim - inicio)

############################################################
# 15. GRAFICOS
#
# (a) BIC em funcao de log(lambda1): o ponto vermelho e
#     o lambda escolhido.
# (b) Coeficientes estimados (azul = selecionada).
############################################################

par(mfrow = c(2, 1))

plot(
  log(lam1_seq),
  bic_seq,
  type = "b",
  pch = 19,
  xlab = expression(log(lambda[1])),
  ylab = "BIC",
  main = "BIC (ajustes validos) ao longo da grade de lambdas"
)

points(
  log(lambda1_bic),
  bic_seq[posicao_bic],
  col = "red",
  pch = 19,
  cex = 1.8
)

barplot(
  beta_covar,
  names.arg = colnames(X),
  las = 2,
  cex.names = 0.4,
  col = ifelse(selecionada, "steelblue", "grey80"),
  border = NA,
  ylab = expression(hat(beta)),
  main = "Coeficientes estimados (azul = selecionada)"
)
abline(h = 0)

par(mfrow = c(1, 1))

############################################################
# 16. SALVA RESULTADOS
############################################################

write.csv(
  resultado,
  "results/resultado_glmnet_splmm_riboflavina_todos.csv",
  row.names = FALSE
)

write.csv(
  selecionadas,
  "results/resultado_glmnet_splmm_riboflavina_selecionadas.csv",
  row.names = FALSE
)

write.csv(
  tabela_bic,
  "results/resultado_glmnet_splmm_riboflavina_BIC_lambdas.csv",
  row.names = FALSE
)

cat("\nArquivos gerados:\n")
cat(" - results/resultado_glmnet_splmm_riboflavina_todos.csv\n")
cat(" - results/resultado_glmnet_splmm_riboflavina_selecionadas.csv\n")
cat(" - results/resultado_glmnet_splmm_riboflavina_BIC_lambdas.csv\n")
