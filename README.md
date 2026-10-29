# Modelo OCP e polarização 2D

Modelo de pesquisa para eletrodo de ferro em solução de NaNO3, com transporte,
especiação, precipitados e dupla camada. **Os parâmetros cinéticos herdados não
foram validados experimentalmente.** A correção dos balanços altera as curvas;
resultados antigos não devem ser usados como referência de precisão.

## Primeira vez? Rodar e visualizar pelo navegador

Leia o [guia passo a passo para GitHub Actions, nuvem e VS Code](docs/COMO_RODAR.md).
A aba **Actions → Executar simulador e gerar graficos → Run workflow** permite
executar o teste pelo navegador e baixar um ZIP com PNG, CSV e uma galeria HTML.
A página comum do GitHub não executa Python; o workflow manual faz esse trabalho.

As saídas agora incluem **OCP vs. tempo, log10(|j|) vs. E, E vs. log10(|j|),
correntes vs. tempo, perfis de concentração em vários instantes, perfis de
especiação e mapas 2D de pH, potencial, concentrações e sólidos**. Abra o
`index.html` da pasta de resultados para navegar pelas figuras. Dados brutos
permanecem disponíveis em CSV/NPZ, e `graficos.py` permite refazer os gráficos
sem resolver o modelo novamente.

## Executar no ambiente preparado

```sh
/workspace/ocp-environment/run.sh python3 /project/main_2d.py --smoke
```

Esse teste pequeno percorre OCP, condicionamento e polarização. Usa composição,
malha e tempos demonstrativos, imprime um aviso de resolução e não representa
uma experiência de corrosão. Resultados ficam em `/workspace/ocp-environment/results`.

Para o domínio original de 1 × 2 cm, 3600 s de OCP e varredura de ±300 mV a 1 mV/s:

```sh
/workspace/ocp-environment/run.sh python3 /project/main_2d.py --nx 40 --ny 40 --grading 3
```

Essa execução longa exige estudo de convergência e pode custar bastante tempo.
A malha geométrica concentra células no eletrodo. Ajuste `--dt`, `--dt-max` e
`--rtol` e compare as curvas; o aviso inicial de resolução não substitui esse estudo.
O padrão usa 16 × 16 volumes e os mesmos tempos físicos, como ponto de partida.
`--help` lista as opções. Um condicionamento de 1 s conecta suavemente o OCP ao
início da varredura, incluindo a corrente capacitiva dessa transição.

A convecção solutal só é ativada ao fornecer coeficientes de densidade medidos
ou justificados, por exemplo `--density-fe VALOR` em kg/mol. Sem eles o modelo
representa solução quiescente. O escoamento usa FEniCS P2/P1 com Oseen transiente,
parede no metal e tração zero nos outros lados. Há limites explícitos para Davies
(I ≤ 0,5 M), Boussinesq (|Δρ|/ρ ≤ 1%) e para a representação exponencial.

## Saídas e retomada

Cada execução usa um diretório novo. `history.csv` contém corrente faradaica,
capacitiva e total, potencial vs SHE, pH de atividade, erro temporal e balanços
incluindo os fluxos do banho. `fields.npz`, `perfil.csv`, mapas e gráficos
preservam os resultados; `config.json` registra hipóteses/parâmetros e
`status.json` distingue conclusão e falha. Quantidades integrais são por metro
de profundidade fora do plano; multiplicar pela profundidade real para obter mol.

O checkpoint é atômico e inclui química, minerais, troca com o banho, potencial,
controle temporal e estado hidrodinâmico. Após uma interrupção:

```sh
/workspace/ocp-environment/run.sh python3 /project/main_2d.py --resume /results/SEU_DIRETORIO/checkpoint.npz
```

A retomada usa a configuração salva. Após falha por dt mínimo, reduza `dt_min`
em `config.json` antes de retomar; mantenha composição e parâmetros físicos.
Não sobrescreva saídas de outra execução. Falhas retornam código não zero e
mantêm histórico parcial. MPI usa um coordenador determinístico: dois processos
não aceleram o solve; somente o rank 0 calcula e escreve. FEniCS usa COMM_SELF.

## Instalação reproduzível fora deste ambiente

```sh
docker build -t ocp-fenics .
docker run --rm -v "$PWD:/project:ro" -v "$PWD/resultados_2d:/results" -w /results ocp-fenics python3 /project/main_2d.py --smoke
```

Crie `resultados_2d` antes do comando. Use a configuração de proxy/CA apropriada
à sua rede; não desabilite TLS ou assinaturas. No ambiente Codex, o script
`/workspace/ocp-environment/install.sh` já fornece a CA da plataforma.

## Testes

```sh
/workspace/ocp-environment/run.sh python3 -m unittest discover -s /project/tests -v
```

Os testes cobrem equilíbrio e atividades, conservação de elementos/carga,
Scharfetter–Gummel, convergência separada em malha e tempo, solução conhecida de
escoamento, execução das fases, MPI e recuperação de checkpoint. Veja
[formulação e hipóteses](docs/MODELO.md) e a [resposta aos 27 achados](docs/CORRECOES.md).
