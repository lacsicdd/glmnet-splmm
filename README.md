# Descobrindo genes que regulam a produção de vitamina B2: uma abordagem estatística e computacional integrada e inovadora

**Autoras**

- Larissa Castro Silva, Graduanda do Curso de Bacharelado em Matemática, DEMAT/UFSJ
- Daniela Carine Ramires de Oliveira, Professora Titular, DEMAT/UFSJ

Este repositório reúne os códigos em R do trabalho: estudo de simulação Monte Carlo, aplicação aos dados reais de expressão gênica da riboflavina e diagnóstico do modelo ajustado com o pacote `lme4`.

## Resumo

A riboflavina (vitamina B2) desempenha papel essencial no metabolismo energético e em processos celulares fundamentais. A identificação de genes associados à sua produção é estratégica para otimizar a síntese biotecnológica por microrganismos como a *Bacillus subtilis*. Contudo, dados genômicos caracterizam-se pela alta dimensionalidade, impondo desafios estatísticos e computacionais à seleção de variáveis. Este trabalho teve como objetivo aprimorar os procedimentos de seleção de efeitos fixos em Modelos Lineares de Efeitos Mistos (LMM), focando na automação e definição eficiente da grade do parâmetro de regularização para redução do custo computacional e aumento da precisão preditiva. Integrou-se o algoritmo `splmm`, baseado em penalização Lasso, e técnicas do pacote `glmnet`, para estabelecer um procedimento automatizado para otimizar o parâmetro de regularização. Estudos de simulação Monte Carlo avaliaram o desempenho da abordagem sob cenários com covariáveis correlacionadas, diferentes magnitudes de efeito e variação no tamanho amostral. Como medidas de desempenho, analisaram-se sensibilidade, especificidade e Erro Quadrático Médio para diversos cenários. Em seguida, a metodologia foi aplicada ao conjunto de dados reais de expressão gênica da riboflavina, utilizando o pacote `lme4` para reajuste, obtenção dos componentes de variância e diagnóstico do modelo. Os resultados confirmam a eficácia da integração glmnet–splmm, apontando caminhos para a redução do tempo computacional e alinhamento com práticas de *green computing*.

**Palavras-chave:** Riboflavina; Seleção de Variáveis; Modelos Mistos; Alta Dimensionalidade; Lasso; Green Computing.

## O algoritmo glmnet–splmm

1. **glmnet:** usado exclusivamente para definir, de forma automática, o valor de referência `lambda_max` e a sequência decrescente de valores candidatos do parâmetro de regularização dos efeitos fixos.
2. **splmmTuning:** recebe essa sequência, ajusta o modelo misto penalizado para cada valor e calcula o BIC.
3. **splmm:** ajuste final com o λ de menor BIC.

O `glmnet` não incorpora a correlação entre medidas repetidas; a estrutura de dependência dos dados entra apenas nas etapas 2 e 3.

## Estrutura do repositório

```
.
├── R/
│   ├── Simulacao.R       # Simulação Monte Carlo (12 cenários, M = 30)
│   ├── Aplicacao.R       # Seleção de variáveis nos dados da riboflavina
│   └── Diagnostico.R     # Ajuste com lme4 e diagnóstico do modelo
├── data/
│   ├── riboflavinv100.csv             # Resposta + 100 genes (71 observações)
│   └── riboflavinv100_structure.txt   # Cepa, tempo e grupo de cada observação
├── results/              # Saídas geradas pelos scripts
├── CITATION.cff
├── LICENSE
└── README.md
```

## Como reproduzir

1. Clone o repositório e abra a pasta como projeto no RStudio (ou use `setwd()` para a raiz do repositório). Os scripts usam caminhos relativos a essa pasta.
2. Instale os pacotes: `splmm`, `glmnet`, `lme4`, `lmerTest`, `ggplot2` e `patchwork`. Os scripts instalam automaticamente os que estiverem faltando.
3. Execute os scripts na ordem:
   - `R/Simulacao.R`: estudo de simulação (semente 2022). Grava em `results/` as métricas por cenário e por réplica. Demora, pois roda o `splmmTuning` em cada réplica.
   - `R/Aplicacao.R`: aplica o algoritmo glmnet–splmm aos dados reais e grava em `results/` os coeficientes, as variáveis selecionadas e o BIC por λ.
   - `R/Diagnostico.R`: reajusta o modelo com as variáveis selecionadas via `lme4` e produz os gráficos e testes de diagnóstico (`results/diagnostico_lme4.jpg`).

## Configurações da simulação

- **Covariáveis:** X\*, X\*\* e X\*\*\* (normal multivariada correlacionada, covariáveis mistas e covariáveis mistas permutadas).
- **Coeficientes:** β\* e β\*\* (efeitos de magnitudes diferentes).
- **Tamanhos amostrais:** n = 30 (nᵢ = 5) e n = 60 (nᵢ = 10).
- **Medidas de desempenho:** sensibilidade, especificidade e Erro Quadrático Médio.

## Observações sobre a aplicação aos dados reais

Com 101 covariáveis e 71 observações, o `Aplicacao.R` difere da simulação em dois pontos, descritos nos comentários do script: `maxIter = 100` (30 na simulação) e a possibilidade de usar a grade do `glmnet` multiplicada pelo número de observações, caso a grade original produza apenas ajustes saturados. O script informa a escala utilizada.

## Dados

Os dados de expressão gênica da riboflavina são de Bühlmann, Kalisch e Meier (2014). O arquivo utilizado contém a resposta (log da taxa de produção de riboflavina) e 100 genes, com 71 observações em 28 cepas. Consulte as condições de uso dos dados originais antes de reutilizá-los.

## Como citar

Use o botão "Cite this repository" do GitHub (arquivo `CITATION.cff`) ou:

> Silva, L. C.; Oliveira, D. C. R. *Descobrindo genes que regulam a produção de vitamina B2: uma abordagem estatística e computacional integrada e inovadora.* Universidade Federal de São João del-Rei, 2026. Código disponível em: https://github.com/SEU-USUARIO/riboflavina-glmnet-splmm

## Referência dos dados

Bühlmann, P.; Kalisch, M.; Meier, L. High-dimensional statistics with a view toward applications in biology. *Annual Review of Statistics and Its Application*, v. 1, p. 255–278, 2014.

## Licença

Código sob licença MIT (ver `LICENSE`).

## Contato

Daniela Carine Ramires de Oliveira, DEMAT/UFSJ: *(daniela@ufsj.edu.br)*
