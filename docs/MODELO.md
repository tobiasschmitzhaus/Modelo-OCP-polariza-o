# Formulação e limites

## Escopo científico

A formulação original não fechava componentes e carga. A nova química usa
volumes finitos conservativos; FEniCS continua responsável pelo escoamento.
Essa mudança evita tratar alcalinidade como um íon artificial com carga -1.
Não há uma correspondência numérica automática com curvas do código anterior.
Constantes herdadas e volumes molares indicativos constam em `cinetica.py` e são
salvos em cada `config.json`. Não há dados suficientes para afirmar precisão
experimental de OCP, corrente de corrosão ou espessura de filme.

Temperatura fixa: 298,15 K. Pressões de referência em atm relativas a 1 atm;
P_O2=0,21 atm determina o banho via Henry. A referência de H2 é 1 atm. Água tem
atividade unitária; H2 formado escapa para um reservatório gasoso idealizado,
sem balanço de bolhas. A atividade de oxigênio local determina sua fugacidade
Henry equivalente no Nernst. Todos os potenciais são vs SHE e Phi=0 no banho.

## Componentes, espécies e alcalinidade

Componentes transportados: Fe(II) dissolvido total, O2, carbono inorgânico total,
Na, nitrato, nitrito total e AT. Espécies locais: H+, OH-, CO2*, HCO3-, CO3²-,
Fe²+, FeOH+, O2, Na+, NO3-, NO2- e HNO2. A atividade é gamma*c/(1000 mol/m³).
A especiação impõe Kw, Ka1, Ka2, hidrólise FeOH+ e dissociação HNO2. pH=-log10(aH).
Davies usa I=0,5 sum(z² c)/1000, resolvida com a especiação a cada passo.
Fora de I≤0,5 M, o modelo falha explicitamente: não aplica um teto artificial.
O solve acoplado admite pH de atividade entre -2 e 16; cruzar limites não é uma
previsão válida. A busca independente de alcalinidade verifica bracket e resíduo.

AT = OH - H + HCO3 + 2 CO3 + FeOH + NO2, em mol de equivalentes/m³.
Nessa convenção HNO2 é a espécie de referência do nitrito. Por isso a redução
de um nitrato aumenta AT em **três** equivalentes: dois OH- e um NO2-.
A identidade de carga é q/F = 2 Fe_total + Na - NO3 - AT = 0.
O Na inicial é ajustado por essa identidade, em vez de manter a pequena carga
residual do banho original. Fe(II) inicial padrão=1e-12 mol/m³. O nitrito traço
padrão=1e-12 mol/m³ permite Nernst finito, é explícito/configurável e requer
estudo de sensibilidade; não é dado experimental. O modo smoke usa 1e-6 para
ambos e está identificado como demonstração. Não são incluídos outros complexos
de ferro, fases poliméricas, atividade da água não unitária ou química de NOx
além de nitrito. Acrescentá-los exige constantes e validação adicionais.

## Transporte e eletricidade

Para cada espécie, discretiza-se J=-D c grad(ln a + z F Phi/RT) + u c com
Scharfetter–Gummel e upwind de advecção. Cada face é calculada uma única vez e
entra com sinais opostos nas duas células. Fluxos de Fe, DIC, nitrito e AT são
combinações estequiométricas desses mesmos fluxos de espécies. A corrente é
F sum(z J), incluindo difusão, atividade e migração. A restrição de
neutralidade e as continuidades determinam Phi: não é um modelo de Poisson
com dupla camada espacialmente resolvida.

Três lados trocam massa com banho de composição e potencial fixos, inclusive
entrada/saída advectiva; o eletrodo fica em x=0. O banho é infinito, não uma
célula fechada. Os balanços globais incluem essa troca. Os testes fechados
sem eletrodo reativo usam paredes impermeáveis e uma referência de potencial;
célula fechada reativa exige contraeletrodo e é rejeitada pela API.

## Interface

Fe/Fe²+: uma única corrente de troca dependente de aFe multiplica ambos os
ramos exponenciais. Com eta=E_interface-E_Nernst a corrente de ferro zera no
seu equilíbrio. Os betas são alpha*n efetivos associados às inclinações de
Tafel herdadas, não coeficientes elementares que precisem somar um.

ORR, HER e nitrato são ramos **catódicos irreversíveis empíricos**, não pares
Butler–Volmer reversíveis completos. Não se adiciona uma segunda HER para H+
ao mesmo processo de formação de H2. O Nernst de nitrato inclui NO3, NO2 e OH
com atividades adimensionais. Não há amônia nem outras vias de redução.

A corrente faradaica é a soma das reações calculadas com uma **mesma** queda
Rf*jF, resolvida implicitamente por raiz escalar monotônica. Rf=h/kappa_filme.
A difusão de O2 pelo filme poroso tem resistência adicional explícita, com
D_eff=D*porosidade^1,5. Os expoentes não são achatados por tanh ou cortados;
iterados fora da representação numérica são rejeitados.

A capacitância compacta é paralela ao ramo faradaico com filme: q=Cdl(E-Phi),
e jC é sua diferença temporal. A continuidade elétrica usa jF+jC. A compensação
é uma adsorção reversível **idealizada de nitrato**: fluxo de NO3 para o líquido
=-jC/F e Γ=Γ_ref+(q-q_inicial)/F. Γ_ref=2*Cdl/F fornece o inventário de referência
em torno do estado inicial; o modelo rejeita Γ<0. Γ é somado ao balanço de N.
Não é uma descrição molecular de adsorção específica nem Gouy–Chapman;
excursões grandes ou troca de contraíon exigem outro modelo de dupla camada.
E_pzc constante seria absorvido na referência de q; apenas diferenças de carga
entram na dinâmica. A corrente capacitiva não é aproximada por Cdl*v_scan.

O eletrodo flutuante satisfaz integral(jF+jC)=0. No condicionamento e na
polarização, a equação global impõe E(t), e a corrente externa resulta dos
mesmos estados e fluxos que entraram no passo aceito.

## Minerais e reações volumétricas

- Fe²+ + 2OH- ⇌ Fe(OH)2(s): precipitação/dissolução cinética usa Q/Ksp de
  atividades, um único Ksp e inventário não negativo do sólido. A dissolução
  devolve Fe e AT à solução. A lei cinética é empírica, explicitamente
  dependente da supersaturação; não impõe equilíbrio instantâneo.
- 4Fe²+ + O2 + 10H2O → 4Fe(OH)3(s) + 8H+: a taxa de extensão consome Fe e O2,
  reduz AT e gera um inventário sólido. Não gera Fe³+ aquoso fictício.
  Essa fase é um produto **irreversível** no escopo atual, sem dissolução ou
  equilíbrio de Fe(III); extrapolar para dissolução ácida/redutora é inadequado.
- A taxa tipo Stumm–Lee usa Fe²+ livre, O2 e OH- com unidades explícitas.
  O coeficiente refere-se à taxa de perda de Fe(II); a extensão que consome O2
  é um quarto dessa taxa. O solve implícito e os balanços/positividade
  substituem o corte de taxa por AT.
- Não se inclui a redução goethita/magnetita: sua reação, equilíbrio e dados
  cinéticos não estavam especificados. A expressão incorreta foi removida,
  não “calibrada” com números inventados.

Os sólidos são imobilizados. No volume adjacente ao metal, seu volume molar por
área define espessura efetiva do filme; no restante ficam como depósito local.
Há uma hipótese explícita de aderência completa no primeiro volume, sem
sedimentação, desprendimento ou transporte de partículas. A quantidade de ferro
é contabilizada uma vez; a espessura não recebe uma fonte independente.
Os volumes molares indicativos e essa hipótese morfológica precisam de calibração;
a execução é rejeitada se sólidos excederem 5% do volume local.

## Escoamento

Oseen transiente com velocidade advectante do último estado aceito é uma
linearização temporal de Navier–Stokes. Tensão=2 mu sym(grad u)-p I, velocidade
nula no metal e tração zero nos lados abertos. Não há âncora de pressão adicional
nessa configuração aberta. Malha FEniCS e malha de volumes usam o mesmo grading.
A densidade é rho0 + sum(b_i*(C_i-C_i_ref)); os b_i têm unidade kg/mol e precisam
ser fornecidos. Sem b_i, assume-se banho quiescente, em vez de inventar uma
relação de densidade a partir da massa molar. A aproximação Boussinesq é rejeitada
acima de 1% e Re é registrado. A corrente versão é adequada a escoamento laminar;
Re alto ou recirculação com forte backflow exigem estudo de estabilidade/contorno,
não certificação automática pelo fato de o solver terminar.

## Controle numérico e validação

Concentrações positivas usam variáveis logarítmicas, sem substituir valores
negativos por pisos nos coeficientes. Os limites de representação logarítmica
são explícitos; resíduo e inventários devem passar mesmo se um limite for ativo.
Minerais têm limites inferiores zero. Resíduos de cada componente são escalados
por sua própria unidade/escala; neutralidade, I e controle elétrico têm escalas
separadas. O solve só é aceito se resíduos e balanços integrais passarem.

Química e escoamento são iterados; cada tentativa hidrodinâmica começa do mesmo
estado aceito. Nenhuma tentativa rejeitada é incorporada. O erro de um passo
é estimado pela diferença entre um passo e dois meios passos implícitos; o
estado mais fino é aceito. O controle inclui componentes, minerais, potencial e
velocidade. dt termina exatamente no fim de cada fase. A corrente registrada é
a média dos fluxos dos dois meios passos aceitos, e não um recálculo após
modificar pH/filme. I, pH e Nernst pertencem ao mesmo estado acoplado.

Os testes demonstram propriedades e soluções conhecidas, não validade universal
das curvas de corrosão. O estudo de convergência incluído verifica ordem espacial
no caso difusivo neutro e ordem temporal para o operador discreto; não demonstra
que uma malha escolhida resolve todas as reações da experiência completa.
Comparar OCP, corrente, h e pH em malhas/passos sucessivos e com dados experimentais
continua necessário. Não foi executada uma simulação convergida de 3600 s.
