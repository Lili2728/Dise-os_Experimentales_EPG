--
  title: "Diatraea saccharalis "
author: "Lili"
format:
  html:
  toc: true
toc-expand: true
toc-depth: 4
toc-location: left
number-sections: true
self-contained: true
code-fold: true
output-file: "ESM"
editor_options:
  chunk_output_type: console
execute:
  warning: false
message: false
echo: true



library(readxl)
library(factoextra)
library(tidyverse)
library(googlesheets4)
library(cowplot)

library(emmeans)
library(multcomp)
library(multcompView)


gs4_deauth()

url <- "https://docs.google.com/spreadsheets/d/1qxamj0Dxgs6epPolsQNWnBK6mXqO20YL/export?format=xlsx"
tmp <- tempfile(fileext = ".xlsx")
download.file(url, tmp, mode = "wb")
fb <- read_excel(tmp, sheet = "datos_organizados")
str(fb)



# Datos por estadio 


huevos <- fb %>% dplyr::filter(Diatraea == "Egg")  %>% droplevels()
larvas <- fb %>% dplyr::filter(Diatraea != "Egg")  %>% droplevels()

#  ANOVA de un factor + diagnósticos + Tuky + barra
analizar <- function(data, factor, color = TRUE){
  
  frm <- reformulate(factor, response = "prob_mortality")
  md <- aov(frm, data = data)
  
  print(anova(md))
  
  mc <- emmeans::emmeans(md, specs = factor) |>
    cld(Letters = letters, reversed = TRUE) |>
    as.data.frame()
  
  print(mc)
  
  p <- ggplot2::ggplot(
    mc,
    ggplot2::aes(x = .data[[factor]], y = emmean, fill = .data[[factor]])
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .group, y = upper.CL + 0.02),
      size = 5
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = "Prob. de mortalidad") +
    ggplot2::guides(fill = "none")
  
  return(p)
}

# ---   tablas tipo   (factores x zafra con Tukey) ---
tabla_factor_zafra <- function(data) {
  data <- droplevels(data)
  md2 <- aov(prob_mortality ~ mort_factor * season, data = data)
  
  emmeans(md2, ~ mort_factor | season) %>%
    cld(Letters = letters, reversed = TRUE) %>%
    as.data.frame() %>%
    dplyr::transmute(season, mort_factor
                     , media = round(emmean, 3)
                     , grupo = gsub(" ", "", .group)) %>%
    dplyr::mutate(valor = paste0(media, " ", grupo)) %>%
    dplyr::select(mort_factor, season, valor) %>%
    tidyr::pivot_wider(names_from = season, values_from = valor)
}


# Análisis univariado para : Huevos

## Factor de mortalidad

plot_h_a <- analizar(huevos, "mort_factor")

## Tipo de bosque

plot_h_b <- analizar(huevos, "type_forest")


## Distancia al bosque


plot_h_c <- analizar(huevos, "distan_Forest", color = FALSE)


## Figura — Huevos


fig_huevos <- list(plot_h_a, plot_h_b, plot_h_c) %>%
  plot_grid(plotlist = ., ncol = 2, labels = "auto")

fig_huevos
print(fig_larvas)


ggsave("Figure_huevos.jpg", plot = fig_huevos
       , units = "cm", height = 20, width = 20)


# Análisis univariado para ; Larvas

## Factor de mortalidad

plot_l_a <- analizar(larvas, "mort_factor")


## Tipo de bosque

plot_l_b <- analizar(larvas, "type_forest")


## Distancia al bosque

plot_l_c <- analizar(larvas, "distan_Forest", color = FALSE)


## Figura — Larvas

fig_larvas <- list(plot_l_a, plot_l_b, plot_l_c) %>%
  plot_grid(plotlist = ., ncol = 2, labels = "auto")

fig_larvas

ggsave("Figure_larvas.jpg", plot = fig_larvas
       , units = "cm", height = 20, width = 20)


# TABLAS


tabla_factor_zafra(huevos) %>%
  knitr::kable(caption = "Tabla - Huevos (media + letras de Tukey)")


## Comparación de factores por zafra (larvas)

#| label: tabla-larvas
tabla_factor_zafra(larvas) %>%
  knitr::kable(caption = "Larvas (media + letras de Tukey)")


## Tablas — Detalle por bosque × distancia × factor

resumen <- fb %>%
  dplyr::group_by(diatraea, type_forest, distan_forest, mort_factor, season) %>%
  dplyr::summarise(
    media = mean(prob_mortality, na.rm = TRUE)
    , epm   = sd(prob_mortality, na.rm = TRUE) / sqrt(dplyr::n())
    , .groups = "drop"
  ) %>%
  dplyr::mutate(valor = sprintf("%.4f \u00b1 %.4f", media, epm)) %>%
  dplyr::select(diatraea, type_forest, distan_forest, mort_factor, season, valor) %>%
  tidyr::pivot_wider(names_from = season, values_from = valor) %>%
  dplyr::arrange(diatraea, type_forest, distan_forest, mort_factor)

resumen %>% dplyr::filter(diatraea == "Egg") %>%
  knitr::kable(caption = "Tabla 3: media \u00b1 EPM (huevos)")

resumen %>% dplyr::filter(diatraea != "Egg") %>%
  knitr::kable(caption = "Tabla 5: media \u00b1 EPM (larvas)")



# PCA  por estadio 

pca_estadio <- function(data) {
  data     <- droplevels(data)
  factores <- levels(data$mort_factor)
  
  dtx <- data %>%
    tidyr::pivot_wider(names_from = mort_factor, values_from = prob_mortality) %>%
    tidyr::drop_na(dplyr::all_of(factores))
  
  keep <- factores[sapply(dtx[factores], function(x) sd(x, na.rm = TRUE) > 0)]
  
  dtx <- dtx %>%
    dplyr::select(type_forest, distan_forest, dplyr::all_of(keep)) %>%
    as.data.frame()
  
  PCA(dtx, scale.unit = TRUE, graph = FALSE, quali.sup = 1:2)
}

# Reporte gráfico del PCA
pca_reporte <- function(mv, titulo) {
  print(get_eigenvalue(mv) %>% round(2) %>%
          knitr::kable(caption = paste("Autovalores -", titulo)))
  
  print(fviz_eig(mv, addlabels = TRUE, barfill = "steelblue", barcolor = "steelblue") +
          labs(title = paste("Scree -", titulo)))
  
  g_var <- fviz_pca_var(mv, col.var = "cos2"
                        , gradient.cols = c("red", "orange", "blue")
                        , repel = TRUE) + labs(title = "Variables (cos2)")
  g_c12 <- fviz_contrib(mv, choice = "var", axes = 1:2, fill = "steelblue") +
    labs(title = "Contribución Dim 1+2")
  print(plot_grid(g_var, g_c12, ncol = 2, labels = "auto"))
  
  g_for <- fviz_pca_ind(mv, habillage = "type_forest", addEllipses = TRUE
                        , label = "none", palette = c("forestgreen", "darkorange")) +
    labs(title = "Por tipo de bosque")
  g_dis <- fviz_pca_ind(mv, habillage = "distan_forest", addEllipses = TRUE
                        , label = "none", palette = c("royalblue", "firebrick")) +
    labs(title = "Por distancia")
  print(plot_grid(g_for, g_dis, ncol = 2, labels = "auto"))
  
  print(fviz_pca_biplot(mv, habillage = "type_forest", addEllipses = TRUE
                        , col.var = "black", repel = TRUE, label = "var"
                        , palette = c("forestgreen", "darkorange")) +
          labs(title = paste("Biplot -", titulo)))
  
  print(dimdesc(mv, axes = 1:2, proba = 0.05))
  
  fig <- plot_grid(g_var, g_for, ncol = 2, labels = "auto")
  ggsave(paste0("PCA_", titulo, ".jpg"), plot = fig
         , units = "cm", width = 30, height = 15)
  fig
}


## PCA ; Huevos


#| label: pca-huevos
#| fig-width: 12
#| fig-height: 5
mv_h <- pca_estadio(huevos)
pca_reporte(mv_h, "huevos")


## PCA - Larvas


#| label: pca-larvas
#| fig-width: 12
#| fig-height: 5
mv_l <- pca_estadio(larvas)
pca_reporte(mv_l, "larvas")

