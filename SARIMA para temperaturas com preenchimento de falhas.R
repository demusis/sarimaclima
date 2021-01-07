library(boot)
library(hexbin)
library(imputeTS)
library(robust)
library(forecast)
library(tseries)
library(fUnitRoots)
library(nortest)
library(openxlsx)
library(portes)
library(tsoutliers)
library(TSPred)

# Carga dos dados
dados <- read.csv("DADOS DAS ESTACOES COM 2018.csv")
head(dados)
sapply(dados, class) # Verificando os tipos

# Análise da temperatura
# Cria objeto do timo Time Series para o INMET
T_INMET <- ts(dados$T_INMET,start=c(1961,1),f=12)
plot(T_INMET)

# Cria objeto do timo Time Series para o ICEA
T_ICEA <- ts(dados$T_ICEA,start=c(1961,1),f=12)
plot(T_ICEA)

# Observe que as séries possuem falhas. 
# Inicialmente as preencheremos com uma regressão robusta entre INMET e ICEA.
mod1 = lmRob(T_INMET ~ T_ICEA, data = dados)
summary(mod1)

bin <- hexbin(T_ICEA, T_INMET, xbins=30)
plot(bin)

# Sem dados do INMET com dados do ICEA
aux <- dados[is.na(dados$T_INMET)&!is.na(dados$T_ICEA),]
head(aux)

# Estima INMET com base no ICEA
aux <- predict(mod1, newdata=dados[is.na(dados$T_INMET)&!is.na(dados$T_ICEA),]) 
head(aux)

# SUbstitui as estimativas na base de dados original
dados[is.na(dados$T_INMET)&!is.na(dados$T_ICEA),'T_INMET'] <- aux

# Cria novamente a TS
T_INMET <- ts(dados$T_INMET,start=c(1961,1),f=12)
plot(T_INMET)

# Preenche as falhas que sobraram pelo estimador de Kalman.
T_INMET <- na_kalman(T_INMET, model = "auto.arima")
plot(T_INMET)


# Testes da estacionariedade da série
# KPSS
urkpssTest(T_INMET, type = c("tau"), lags = c("short"),use.lag = NULL, doplot = TRUE)

# Dickey-Fuller
adf.test(T_INMET, alternative="stationary")

# ACF e ACF parcial
acf(T_INMET)
pacf(T_INMET)

# Se rejeitarmos a hipótese da estacionariedade, busca-se o valor de d.
# Calculando as correlações para uma diferença de ordem 1
acf(diff(T_INMET, differences=1))
pacf(diff(T_INMET, differences=1))

# Como o ACF cai depois do primeiro lag, podemos partir de p = 1-3. 
# Para o PACF o q poderia ser igual a 2, temos, então, um 
# ARIMA(p, d, q) = ARIMA (1, 1, 2) com os mesmos parâmetros para sazonalidade.

# Buscando melhor modelo via auto.arima
l <- BoxCox.lambda(T_INMET)
ajuste_auto <- auto.arima(T_INMET, 
                          lambda = l,
                          stepwise=FALSE, 
                          approximation=FALSE, 
                          trace=TRUE)
ajuste_auto


# Pesquisa de Outliers
outliers <- tso(BoxCox(T_INMET,l),
                tsmethod = "arima", 
                args.tsmethod = list(order=c(2, 0, 1), 
                                     seasonal=list(order=c(0, 1, 1), 
                                                   period=12)), 
                types = c("AO","LS","TC"), maxit.iloop=50)
outliers

# Substitui o outlier por NA
dados[outliers$outliers$ind, 'T_INMET'] <- NA

# Substituindo os NAs por estimativas...
dados$T_INMET <- na_kalman(dados$T_INMET, model = "auto.arima")

# Repete o procedimento desde o começo até estabilizar (+- 3 iterações)...


# Bootstrap
AIC_Boot <- function(ts) {
  aux <- Arima(ts, 
               lambda=l, 
               order=c(2, 0, 1), 
               seasonal=list(order=c(0, 1, 1), 
                             period=12))
  return(aux$aic)
}

aux <- Arima(T_INMET, 
             lambda=l, 
             order=c(2, 0, 1), 
             seasonal=list(order=c(0, 1, 1), 
                           period=12))
AIC_Boot(T_INMET)

# Bootstrap estacionário com bloco de 20
boot_ajuste_auto_bc_outliers <- tsboot(T_INMET, 
                                       AIC_Boot, 
                                       R = 1000, 
                                       l = 20, 
                                       sim = "geom")
boot_ajuste_auto_bc_outliers
boot.ci(boot_ajuste_auto_bc_outliers)