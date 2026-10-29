# Resposta aos 27 achados

A revisão original se refere ao commit 872afb6. Esta alteração reestrutura a
formulação para tornar seus balanços explícitos. “Tratado” abaixo significa
correção de implementação ou substituição de uma hipótese inconsistente por
uma hipótese declarada; **não** significa validação experimental de todas as
constantes. Os limites científicos estão em [MODELO.md](MODELO.md).

| Nº | Tratamento no código e evidência |
|---|---|
| 1 | O solve químico foi substituído por um único solve acoplado escalado; não existe sobrescrita de NewtonSolver. O escoamento usa a API variacional FEniCS correta. Execução completa curta cobre a transição. |
| 2 | Um único controlador percorre fases com avanço de tempo e término exato. Não há laço duplicado. |
| 3 | Histórico incremental único alimenta gráficos; listas globais indefinidas foram eliminadas. Dados brutos são persistidos. |
| 4 | Falha por resíduo/erro/dt mínimo impede avançar de fase, retorna erro e salva checkpoint/histórico. Teste provoca falha e retoma. |
| 5 | Ambos os ramos do ferro compartilham i0(aFe) e eta relativa ao mesmo Nernst. Teste exige corrente Fe zero em eta=0. |
| 6 | Neutralidade algébrica e fluxos por todas as espécies determinam Phi e a corrente total. Difusão e gradiente de atividade entram; teste Boltzmann verifica equilíbrio eletrodifusivo. |
| 7 | NO3 → NO2 conserva N. HNO2/NO2 são especiados e seu Nernst usa atividades; AT inclui a base nitrito. Testes verificam N/carga. |
| 8 | Uma HER efetiva substitui os dois ramos duplicados. É uma cinética irreversível empírica com referência gasosa declarada, não uma soma de H+ e água sem reagentes. |
| 9 | Expressão injustificada de goethita/magnetita removida. Filme deriva dos inventários minerais adjacentes ao metal, com volumes molares explícitos. A redução mineral exige uma extensão futura com dados. |
| 10 | Inventário Fe(OH)2 e dissolução devolvendo Fe/AT; mesmo Q de atividades e Ksp em toda parte. Testes de massa/estequiometria. |
| 11 | Oxidação gera Fe(OH)3 sólido em vez de Fe3+ aquoso incompatível com AT. Balanço 4:1:8:4 explícito e testado; produto é irreversível no escopo atual. |
| 12 | q=Cdl(E-Phi), jC por diferença temporal, continuidade jF+jC e inventário de contraíon na dupla camada. Teste de rede RC no limite linear e balanço global de N. |
| 13 | Uma raiz escalar resolve a queda Rf*jF compartilhada. Teste verifica a soma das correntes e a mesma queda nos ramos. |
| 14 | pH usa atividade, assim como Nernst e solubilidade. Teste de água a I=0,1 M retorna pH=7. Pressões/referências gasosas explícitas. |
| 15 | Fe total é separado de Fe2+ livre e FeOH+. Cinética/solubilidade usam a espécie livre; o componente transportado é a soma. Teste verifica massa e hidrólise. |
| 16 | Fluxos de DIC/AT são combinações estequiométricas dos fluxos de espécies, sem carga/D fictícios. Cancelamento interno de fluxos é testado. |
| 17 | I usa soma de todas as espécies, resolvida com a especiação. Davies fora de 0–0,5 M é erro explícito, não clipping. |
| 18 | Bisseção verifica bracket e resíduo; casos fora do domínio falham, sem retornar silenciosamente pH limite. Testes incluem a falha antiga em AT=-1000. |
| 19 | Oseen transiente inclui inércia advectiva e tensão simétrica; contornos abertos têm tração zero e dispensam âncora de pressão. Teste usa solução incompressível conhecida; não se afirma turbulência/backflow geral. |
| 20 | Densidade depende de coeficientes solutais fornecidos, com estado de referência e limite Boussinesq. Sem coeficientes o banho é quiescente. Nenhuma constante empírica de densidade foi inventada. |
| 21 | Hidrodinâmica atualizada em cada passo de ambas as fases; tentativas não alteram estado aceito e checkpoint inclui fluido. Acoplamento Picard verifica convergência. |
| 22 | Concentrações logarítmicas, SG e advecção upwind eliminam o remendo c_safe. Malha graduada, aviso de resolução e teste espacial de segunda ordem no caso difusivo. A malha da experiência ainda precisa de estudo próprio. |
| 23 | Step doubling estima erro temporal; rejeição/restauração conjunta, dt mínimo e final de fase exato. Teste temporal isola a ordem de Euler implícito. |
| 24 | Escalas separadas para cada componente, neutralidade, força iônica, minerais e controle elétrico; aceitação exige resíduo e balanços integrais. |
| 25 | Sem tanh exponencial ou corte de Stumm–Lee por AT. Exponenciais fora de domínio falham; reações implícitas preservam disponibilidade/positividade até tolerâncias verificadas. A lei de precipitação empírica é declarada. |
| 26 | Um coordenador MPI escreve saídas; demais ranks aguardam resultado. Não se concatenam DOFs/ghosts. CSV incremental, mapas, perfis, checkpoint atômico e status. Teste exige equivalência serial/MPI e retomada sem duplicação de tempos. |
| 27 | pH/I/Nernst e cinética são avaliados no mesmo estado acoplado. Correntes vêm dos dois meios passos efetivamente aceitos; não são recalculadas após alterar a química. |

## Validação e limites

Os testes incluem execução curta das três fases, soluções conhecidas de
transporte/escoamento, equilíbrio, balanços, refinamentos separados em tempo e
espaço, MPI e recuperação. Não foram validadas contra experimento as constantes
cinéticas, volumes molares, espessuras, referência de gás H2 ou densidade solutal.
A simulação padrão de 3600 s não foi executada até convergência física nesta
alteração. Para uso quantitativo, são necessários dados experimentais e estudos
de malha/tempo do caso real, incluindo sensibilidade à composição traço inicial.
