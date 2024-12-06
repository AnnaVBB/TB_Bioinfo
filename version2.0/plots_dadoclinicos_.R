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
# GRÁFICO DE PIZZA PARA GÊNERO 
  # Criar o gráfico com barras mostrando a contagem para coad
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
print(gender_plot) #TEMOS 21 mulheres e 20 homens na amostra
-------------------------------------------
 # Criar o gráfico com barras mostrando a contagem para luad
  gender_plot <- ggplot(clinical_luad_filtered, aes(x = gender)) +
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

---------------------------------------------------------------------------------------------------------------------
#ANALISES RELACIONADAS DOS DADOS CLINICOS
  
  # Gráfico de barras empilhadas para status vital por estágio para LUAD
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
-----------------------------------------------------------------------
 
  # Gráfico de barras empilhadas para status vital por estágio para COAD
  # Carregar pacotes necessários
library(ggplot2)
library(dplyr)

# Suponha que os dados estão no data frame chamado "clinical_data"
# As colunas relevantes são "ajcc_pathologic_stage" e "vital_status"

# Verificar e preparar os dados
stage_survival_counts <- clinical_coad_final %>%
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


