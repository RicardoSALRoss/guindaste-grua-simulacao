# Simulação de Guindaste de Torre (Grua)

Este repositório contém o código de modelagem e simulação de um guindaste de torre (Grua), modelado como um **pêndulo duplo em 3D**.

O projeto foi desenvolvido como parte da disciplina **PME3380 - Modelagem de Sistemas Dinâmicos** (2025) da Escola Politécnica da USP.

## Descrição do Projeto

O objetivo deste projeto é analisar o comportamento dinâmico de uma grua submetida a movimentações na base e distúrbios externos. O código abrange:

* **Estipulação dos parâmetros** físicos do sistema.
* **Modelagem Algébrica** utilizando o *Symbolic Math Toolbox* do MATLAB para obtenção das equações de movimento (Lagrange).
* **Linearização** do sistema em torno de pontos de equilíbrio.
* **Análise no Domínio da Frequência** (Diagramas de Bode) para avaliar a resposta dos elos ($\theta$ e $\phi$).
* **Simulação** da resposta temporal do sistema.

## Tecnologias Utilizadas

* MATLAB (R202x)
* Symbolic Math Toolbox
* Control System Toolbox (para os diagramas de Bode)

## Como Executar

1.  Clone este repositório:
    ```bash
    git clone [https://github.com/seu-usuario/guindaste-grua-simulacao.git](https://github.com/seu-usuario/guindaste-grua-simulacao.git)
    ```
2.  Abra o MATLAB na pasta do projeto.
3.  Execute o script principal (ex: `main.m` ou `run_simulation.m`) para gerar os gráficos e resultados.

##  Autores

Grupo de Engenharia Mecânica:

* Jun Moritani
* Pedro Henrique
* Ricardo Ross
* Vinícius Melo

---
*Escola Politécnica da Universidade de São Paulo (Poli-USP) - 2025*
