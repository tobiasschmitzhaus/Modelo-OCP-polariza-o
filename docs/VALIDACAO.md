# Evidências de verificação

Foram verificados **20 casos de teste**, todos aprovados após os ajustes:

- 13 testes de química, cinética, transporte, corrente e rede RC, em
  `tests/test_physics.py`.
- 2 testes de escoamento: repouso/limite Boussinesq e solução polinomial
  incompressível conhecida, em `tests/test_flow.py`.
- 2 testes de convergência, em `tests/test_convergence.py`: redução do erro
  espacial compatível com ordem 2 no caso difusivo neutro e ordem temporal 1
  de Euler implícito contra a evolução exata do operador discreto.
- 3 testes de execução, em `tests/test_workflow.py`: três fases com CSV,
  mapas e checkpoint, falha controlada/retomada, e igualdade entre serial e
  MPI com dois processos. A falha é provocada por exigir erro temporal muito
  menor que o admissível pelo dt mínimo, depois retoma com tolerâncias
  numéricas adequadas. As constantes físicas permanecem iguais.

Uma rodada integrada inicialmente detectou que um teste dependia de uma
rejeição espontânea, que deixou de ocorrer após melhorar a solução. O teste
foi corrigido para provocar a falha explicitamente, e seus três casos foram
reexecutados e aprovados. Não se contabilizou essa rodada falha como aprovada.
Os demais 17 casos e os três casos reexecutados cobrem os 20 testes acima.

Também foram exercitados:

- Reconstrução do Dockerfile versionado e execução do script de instalação,
  com Poisson/MUMPS de solução conhecida e verificações de química.
- Demonstração das fases com acoplamento químico/hidrodinâmico FEniCS e
  encerramento limpo de MPI. O coeficiente de densidade usado nessa
  verificação é apenas um dado de teste, não calibração física do eletrólito.
- Composição traço original: o primeiro passo implícito foi verificado e o
  controlador detectou a necessidade de reduzir dt. Uma tentativa com tempo
  OCP de 2e-4 s foi interrompida pelo limite de execução de 180 s após avançar
  aproximadamente 4,6e-5 s. Isso é uma limitação de custo do diagnóstico,
  não uma simulação completa validada. A inicialização muito diluída é rígida.
  Um diagnóstico subsequente, com 10 µs de OCP, 10 µs de condicionamento e
  2 µs de polarização, completou as três fases com a composição traço original
  e controle adaptativo; log `trace-short-final.log`, saída `trace-short-final`.

As verificações executam em Debian 12, Python 3.11, DOLFIN 2019.2.0.dev0,
NumPy 1.24.2, SciPy 1.10.1, Matplotlib 3.6.3 e PETSc 3.18.5. Cada nova execução
registra as versões efetivas em config.json. Arquivos de logs e resultados desta
sessão estão em `/workspace/ocp-environment/validation` e `results`.

**Não foi realizada validação experimental nem simulação convergida de 3600 s.**
A malha da demonstração é deliberadamente pequena e aciona um aviso de resolução.
Aprovar os testes não determina que a malha padrão resolve a camada difusiva ou
que as constantes empíricas predizem o sistema experimental do usuário.
