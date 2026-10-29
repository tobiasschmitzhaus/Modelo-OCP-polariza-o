# Como executar e ver os gráficos

Você não precisa instalar o VS Code para usar o simulador. Existem três opções.
Comece pelo teste curto: ele confere a instalação e gera exemplos de todos os
gráficos, mas usa tempos, composição e malha demonstrativos.

## 1. No navegador, pelo GitHub Actions

A página comum do GitHub e seu editor não executam Python. Este repositório
agora inclui uma execução manual no **GitHub Actions**, que roda em uma máquina
na nuvem e disponibiliza um ZIP com os resultados.

1. Abra o [repositório](https://github.com/tobiasschmitzhaus/Modelo-OCP-polariza-o).
2. Clique na aba **Actions**.
3. Se o GitHub pedir para ativar Actions, habilite os workflows do repositório.
4. Na lista à esquerda, escolha **Executar simulador e gerar graficos**.
5. Clique em **Run workflow**, escolha a branch `main` e mantenha o modo `teste`.
6. Clique no botão verde **Run workflow** para iniciar.
7. Abra a execução que aparecer e aguarde sua conclusão.
8. Na seção **Artifacts**, baixe `resultados-NUMERO`, extraia o ZIP e abra
   `execucao/index.html` no navegador. Também pode abrir cada PNG diretamente.

Nenhum pacote Python precisa ser instalado no seu computador para essa opção.
É necessário estar conectado ao GitHub e ter acesso para executar Actions.
O uso está sujeito à disponibilidade/cota de Actions da sua conta. Esta
configuração foi adicionada ao repositório; uma execução real em Actions ainda
precisa ser iniciada pelo botão acima. Os comandos Docker equivalentes foram
exercitados no ambiente de desenvolvimento.

Depois do teste, o modo `simulacao` usa a malha padrão e permite informar o tempo
de OCP em segundos (`3600` significa uma hora de tempo físico). Não confunda
esse valor com duração computacional: o cálculo pode levar muitas horas.
O job está limitado a 6 horas; uma interrupção pode deixar somente resultados
parciais e o último checkpoint, identificados pelo `status.json`. Os artefatos
ficam disponíveis por 14 dias. Resultados científicos exigem refinamento e
calibração conforme [MODELO.md](MODELO.md).

## 2. Neste ambiente de nuvem do Codex

O ambiente desta conversa já tem as dependências. Você pode pedir no chat
"rode o teste curto e mostre os gráficos", ou usar o terminal:

```sh
/workspace/ocp-environment/run.sh python3 /project/main_2d.py --smoke
```

O programa imprime o diretório da execução. No disco, ele fica dentro de
`/workspace/ocp-environment/results/resultados_2d/`. Abra seus PNG ou `index.html`.
O caminho `/workspace/ocp-environment/run.sh` pertence a este ambiente de nuvem;
ele não existe automaticamente no seu computador ou no GitHub Actions.

Para a simulação com os tempos físicos padrão e uma malha maior:

```sh
/workspace/ocp-environment/run.sh python3 /project/main_2d.py --nx 40 --ny 40 --grading 3
```

Esta execução longa não foi validada até 3600 s. Use o teste primeiro e consulte
as limitações físicas e numéricas no README e em MODELO.md.

## 3. No seu computador, com VS Code e Docker

O VS Code é o editor; quem executa o programa é o Python dentro do contêiner
Docker com FEniCS e as outras bibliotecas. Em Windows/macOS, instale e abra o
**Docker Desktop**. No Linux, use Docker Engine. Somente instalar Python no
VS Code não instala FEniCS/DOLFIN.

1. No GitHub, clique em **Code → Download ZIP** e extraia o conteúdo.
2. No VS Code, abra a pasta extraída em **File → Open Folder**.
3. Abra **Terminal → New Terminal**. Confirme que o terminal está na pasta que
   contém `Dockerfile`, `main_2d.py` e `README.md`.
4. Construa o ambiente. Na primeira vez isso baixa as bibliotecas:

   ```sh
   docker build -t ocp-fenics .
   ```

5. Crie uma pasta de resultados:

   ```sh
   mkdir resultados_2d
   ```

6. Execute o teste. O comando abaixo, em uma linha, funciona no PowerShell e
   em terminais usuais de Linux/macOS:

   ```sh
   docker run --rm -e HOME=/tmp -v "${PWD}:/project:ro" -v "${PWD}/resultados_2d:/results" -w /results ocp-fenics python3 /project/main_2d.py --smoke
   ```

7. Abra `resultados_2d` no VS Code. Dentro da pasta da execução, clique nos PNG
   ou abra `index.html` em seu navegador.

Para o caso longo, substitua `--smoke` por `--nx 40 --ny 40 --grading 3`.
No Linux, o contêiner pode criar arquivos pertencentes a root; você pode
adicionar `--user "$(id -u):$(id -g)"` após `--rm` para usar seu usuário.

## Quais arquivos são gerados?

| Arquivo | Conteúdo |
|---|---|
| `index.html` | Galeria local de todos os gráficos e links para os CSV |
| `ocp_vs_tempo.png` e `ocp.csv` | OCP vs. tempo, somente na fase de circuito aberto |
| `potencial.png` | Potencial vs. tempo nas três fases |
| `log_i_vs_E.png` | log10 da magnitude da corrente vs. potencial; ramos anódico/catódico identificados |
| `E_vs_log_i.png` | Apresentação usual com E no eixo vertical e log10\|j\| no horizontal |
| `polarizacao.png` e `polarizacao.csv` | Corrente assinada e dados da varredura; também inclui corrente faradaica/capacitiva no CSV |
| `correntes_vs_tempo.png` | Correntes faradaica, capacitiva e total vs. tempo |
| `pH_vs_tempo.png` | pH mínimo e máximo no domínio vs. tempo |
| `perfis_concentracao.png` | Fe(II) total, O2, DIC, Na+, NO3- e nitrito total vs. distância, em vários instantes |
| `perfis_especies.png` | Especiação final: ferro livre/complexado, H+/OH-, carbonatos e nitrito/ácido nitroso |
| `perfis_pH_potencial.png` | Perfis temporais de pH e potencial da solução |
| `mapa_pH.png`, `mapa_potencial.png` | Mapas 2D no último estado salvo |
| `mapas_concentracao.png`, `mapas_solidos.png` | Mapas 2D de componentes e precipitados no último estado salvo |
| `perfil.csv`, `perfil_especies.csv` | Valores finais dos perfis na meia altura, em unidades declaradas |
| `snapshots/`, `perfis_temporais/` | Campos em NPZ e perfis em CSV nos instantes armazenados |
| `history.csv`, `config.json`, `status.json`, `checkpoint.npz` | Histórico, parâmetros, conclusão/falha e retomada |

Corrente catódica é negativa e logaritmo real de número negativo não existe.
Por isso os gráficos logarítmicos usam **log10(|j|)**, com j em A/m², e excluem
os zeros em vez de inventar um piso. Eles não ajustam automaticamente retas
nem extraem parâmetros de Tafel: a existência do gráfico não comprova regime
cinético controlado por ativação.

Por padrão, os perfis são salvos no início, em 25%, 50% e 100% do tempo de OCP,
e no fim da varredura. Em uma simulação nova, personalize com, por exemplo:

```sh
python3 main_2d.py --profile-times 0 60 300 600 1800 3600
```

Esse último comando deve ser executado **dentro do ambiente com as bibliotecas**;
no Docker, substitua o final `--smoke` do comando de execução por essas opções.
O início e o fim da execução sempre são incluídos. O integrador termina um
passo em cada instante solicitado, sem confundir um estado posterior com o
instante nominal. A quantidade de arquivos cresce com os instantes escolhidos.

## Refazer as figuras sem repetir o cálculo

Se os resultados já existem, basta ler os arquivos salvos:

```sh
python3 graficos.py resultados_2d/PASTA_DA_EXECUCAO
```

No ambiente desta conversa:

```sh
/workspace/ocp-environment/run.sh python3 /project/graficos.py /results/resultados_2d/PASTA_DA_EXECUCAO
```

No computador local, use o mesmo comando `docker run`, substituindo o final por
`python3 /project/graficos.py /results/PASTA_DA_EXECUCAO`.
Resultados antigos sem snapshots permitem refazer as figuras do estado final;
para obter perfis históricos que não foram salvos, é necessário executar de novo.
