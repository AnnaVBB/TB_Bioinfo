# Gráfico de barras para o estágio patológico
library(ggplot2)


#GRAFICO DE ESTAGIO DO CANCER PULMAO
ggplot(clinical_luad_filtered, aes(x = ajcc_pathologic_stage)) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Distribuição dos Estágios Patológicos",
    x = "Estágio Patológico (AJCC)",
    y = "Frequência"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
-----------------------------------------------------------------------------
# Gráfico de pizza para gênero
  # Criar o gráfico com barras mostrando a contagem
  gender_plot <- ggplot(clinical_coad_final, aes(x = gender)) +
  geom_bar(aes(fill = gender), stat = "count") +  # Usa 'stat = "count"' para mostrar as contagens
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +  # Adiciona a contagem acima das barras
  scale_fill_manual(values = c("pink", "skyblue")) +  # Cores personalizadas
  labs(
    title = "Distribuição por Gênero no TCGA_COAD",
    x = "Gênero",
    y = "Contagem",
    fill = "Gênero"
  ) +
  theme_minimal(base_size = 14)

# Exibir o gráfico
print(gender_plot) #TEMOS 33 mulheres e 24 homens na amostra
-------------------------------------------
# Gráfico de barras com porcentagem
# Calcular porcentagens para status vital
vital_table <- table(clinical_luad_filtered$vital_status)
vital_df <- data.frame(
  vital_status = names(vital_table),
  count = as.vector(vital_table)
)
vital_df$percentage <- round(100 * vital_df$count / sum(vital_df$count), 1)

# Exibir tabela de porcentagens (opcional)
print(vital_df)
#TEM 54,4% vivos e 45,6 mortos 

#Gráfico de barras com porcentagem

library(ggplot2)

ggplot(vital_df, aes(x = vital_status, y = percentage, fill = vital_status)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_text(aes(label = paste0(percentage, "%")), vjust = -0.5) + # Adiciona os valores percentuais
  labs(
    title = "Distribuição do Status Vital (Porcentagem)",
    x = "Status Vital",
    y = "Porcentagem",
    fill = "Status Vital"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("darkgreen", "red"))



---------------------------------------------------------------------------------------------------------------------
#ANALISES RELACIONADAS DOS DADOS CLINICOS
  
  # Gráfico de barras empilhadas para status vital por estágio
  ggplot(clinical_luad_filtered, aes(x = ajcc_pathologic_stage, fill = vital_status)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Status Vital por Estágio Patológico",
    x = "Estágio Patológico",
    y = "Proporção",
    fill = "Status Vital"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#COAD
# Gráfico de barras empilhadas para status vital por estágio
ggplot(clinical_coad_final, aes(x = ajcc_pathologic_stage, fill = vital_status)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Status Vital por Estágio Patológico",
    x = "Estágio Patológico",
    y = "Proporção",
    fill = "Status Vital"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
------------------------------------------------------------
  # Carregar pacotes necessários
library(ggplot2)
library(dplyr)

# Suponha que os dados estão no data frame chamado "clinical_data"
# As colunas relevantes são "ajcc_pathologic_stage" e "vital_status"

# Verificar e preparar os dados
stage_survival_counts <- clinical_luad_filtered %>%
  filter(!is.na(ajcc_pathologic_stage), !is.na(vital_status)) %>%  # Remover valores ausentes
  group_by(ajcc_pathologic_stage, vital_status) %>%  # Agrupar por estágio e status de sobrevida
  summarise(count = n(), .groups = 'drop')  # Contar pacientes

# Exibir a tabela com a contagem
print(stage_survival_counts)

# Criar o gráfico de barras empilhadas
ggplot(stage_survival_counts, aes(x = reorder(ajcc_pathologic_stage, -count), y = count, fill = vital_status)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Distribuição de Pacientes por Estágio Patológico e Sobrevida",
       x = "Estágio Patológico",
       y = "Número de Pacientes",
       fill = "Status de Sobrevida") +
  scale_fill_manual(values = c("Alive" = "#56B4E9", "Dead" = "#D55E00")) +  # Cores customizadas
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


--------------------------------------------
  # Boxplot para idade por gênero
  ggplot(clinical_luad_filtered, aes(x = gender, y = age_at_index, fill = gender)) +
  geom_boxplot() +
  labs(
    title = "Idade por Gênero",
    x = "Gênero",
    y = "Idade no Diagnóstico (anos)"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


