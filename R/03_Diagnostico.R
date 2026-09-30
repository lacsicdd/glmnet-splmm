############################################################
# DIAGNOSTICO DO MODELO MISTO (lme4) - DADOS DE RIBOFLAVINA
#
# Modelo com as 15 covariaveis selecionadas
# (tempo + 14 genes) e intercepto e inclinacao (tempo)
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
# 1. PACOTES (somente os necessarios)
#
#  lme4     : ajuste do modelo misto
#  lmerTest : carrega o lme4 e acrescenta p-valores ao
#             summary(); a estimacao e a mesma do lme4
#  ggplot2 / patchwork : graficos de diagnostico
############################################################

pacotes <- c("lme4", "lmerTest", "ggplot2", "patchwork")

for (pk in pacotes) {
  if (!requireNamespace(pk, quietly = TRUE)) install.packages(pk)
}

library(lme4)
library(lmerTest)
library(ggplot2)
library(patchwork)

############################################################
# 2. VARIAVEIS SELECIONADAS (14 genes; o tempo entra a parte)
############################################################

genes_sel <- c(
  "YHZA_at", "YRPE_at", "YCGN_at", "YXLE_at", "ARGF_at",
  "XLYA_at", "YTCF_at", "YTGA_at", "YTGB_at", "PCKA_at",
  "YCKE_at", "YHFH_r_at", "GUAB_at", "TRXA_at"
)

############################################################
# 3. LEITURA E PREPARO DOS DADOS
#
# CSV: linha 1 = y (q_RIBFLV); linhas 2-101 = 100 genes;
# colunas 2-72 = 71 observacoes. Os genes estao nas linhas,
# por isso a matriz e transposta.
############################################################

bruto <- read.csv("data/riboflavinv100.csv", header = TRUE)

y_orig <- as.numeric(unlist(bruto[1, -1]))

x <- t(as.matrix(bruto[2:nrow(bruto), -1]))
colnames(x) <- as.character(bruto[2:nrow(bruto), 1])

stopifnot(all(genes_sel %in% colnames(x)))

# Estrutura: cepa (id) e tempo de cada observacao
estrutura <- read.table("data/riboflavinv100_structure.txt", header = TRUE)

stopifnot(nrow(estrutura) == length(y_orig),
          nrow(estrutura) == nrow(x))

# Covariaveis padronizadas (media 0, dp 1); y na escala original
dat <- data.frame(
  yorig      = y_orig,
  time       = as.numeric(scale(estrutura$time)),
  scale(x)[, genes_sel],
  namestrain = as.factor(estrutura$id)
)

cat("Observacoes:", nrow(dat),
    "| Cepas:", nlevels(dat$namestrain), "\n")

############################################################
# 4. AJUSTE DO MODELO MISTO
############################################################

formula_mod <- reformulate(
  termlabels = c("time", genes_sel, "(time | namestrain)"),
  response   = "yorig"
)

mod <- lmer(formula_mod, data = dat)

summary(mod)

## Variancia total estimada
SQT <- sum((dat$yorig - mean(dat$yorig))^2)
SQT
VTE <- 1 - SQE/SQT
VTE

# Ajuste singular (variancia ~0 ou correlacao +-1 nos efeitos
# aleatorios)? TRUE indica que a estrutura pode estar
# superparametrizada.
cat("\nAjuste singular:", isSingular(mod), "\n")

############################################################
# 5. QUALIDADE DO AJUSTE
############################################################

yhat <- fitted(mod)                    # inclui efeitos aleatorios
SQE  <- sum((dat$yorig - yhat)^2)      # soma dos quadrados dos erros
EQM  <- SQE / nrow(dat)                # erro quadratico medio

cat("SQE:", round(SQE, 4), "| EQM:", round(EQM, 4), "\n")

############################################################
# 6. QUANTIDADES PARA O DIAGNOSTICO
############################################################

resid_pearson <- residuals(mod, type = "pearson")

ef <- ranef(mod)$namestrain
efeitos_df <- data.frame(
  intercept = ef[, "(Intercept)"],
  time      = ef[, "time"]
)

res_df <- data.frame(
  ajustado = as.numeric(yhat),
  residuo  = as.numeric(resid_pearson),
  cepa     = dat$namestrain
)

tema <- theme_minimal()

############################################################
# 7. GRAFICOS DE DIAGNOSTICO
############################################################

# (a) Residuos x valores ajustados (homocedasticidade)
g1 <- ggplot(res_df, aes(ajustado, residuo)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted",
       x = "Fitted values", y = "Pearson residuals") + tema

# (b) QQ-plot dos residuos (normalidade dos erros)
g2 <- ggplot(res_df, aes(sample = residuo)) +
  stat_qq() + stat_qq_line() +
  labs(title = "Normal QQ-plot: Residuals",
       x = "Theoretical Quantiles", y = "Sample Quantiles") + tema

# (c) Histograma dos residuos
g3 <- ggplot(res_df, aes(residuo)) +
  geom_histogram(aes(y = after_stat(density)), bins = 12,
                 fill = "grey80", colour = "black") +
  geom_density() +
  labs(title = "Histogram of Residuals",
       x = "Pearson residuals", y = "Density") + tema

# (d) Residuos por cepa
g4 <- ggplot(res_df, aes(cepa, residuo)) +
  geom_boxplot() +
  labs(title = "Residuals by Strain",
       x = "Strain", y = "Residuals") +
  tema + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# (e) QQ-plot do intercepto aleatorio
g5 <- ggplot(efeitos_df, aes(sample = intercept)) +
  stat_qq() + stat_qq_line() +
  labs(title = "Normal QQ-plot: Random Intercept",
       x = "Theoretical Quantiles", y = "Sample Quantiles") + tema

# (f) QQ-plot da inclinacao aleatoria (tempo)
g6 <- ggplot(efeitos_df, aes(sample = time)) +
  stat_qq() + stat_qq_line() +
  labs(title = "Normal QQ-plot: Random Slope",
       x = "Theoretical Quantiles", y = "Sample Quantiles") + tema

# (g) Contorno da normalidade bivariada dos efeitos aleatorios
g7 <- ggplot(efeitos_df, aes(intercept, time)) +
  geom_point(alpha = 0.6) +
  geom_density_2d(linewidth = 1) +
  labs(title = "Bivariate Kernel Density Contour Plot",
       x = "Random Effect: Intercept",
       y = "Random Effect: Time") + tema

painel <- (g1 | g2 | g3) / (g4 | g5 | g6) / (g7 | plot_spacer() | plot_spacer())

print(painel)

ggsave("results/diagnostico_lme4.jpg", painel,
       width = 12, height = 12, dpi = 300)

############################################################
# 8. TESTES DE NORMALIDADE (Shapiro-Wilk)
############################################################

cat("\n--- Shapiro-Wilk ---\n")
cat("\nResiduos:\n");             print(shapiro.test(resid_pearson))
cat("\nIntercepto aleatorio:\n"); print(shapiro.test(efeitos_df$intercept))
cat("\nInclinacao aleatoria:\n"); print(shapiro.test(efeitos_df$time))
