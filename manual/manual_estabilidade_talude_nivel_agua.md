# Manual do Usuário — Análise de Estabilidade de Talude com Nível d'Água

*Método de Bishop Simplificado em tensões efetivas — teoria, código Python e passo a passo*

- **Notebook de referência:** `analise_estabilidade_talude_bishop_nivel_agua.ipynb`
- **Extensão de:** `analise_estabilidade_talude_bishop.ipynb` (talude seco)
- **Material de apoio didático — Geotecnia**

![resultado do notebook: superfície de ruptura crítica com nível d'água (FS = 1,357) e sem água (FS = 1,832).](manual_estabilidade_talude_nivel_agua_imagens/cell19.png)

*Figura 1 — resultado do notebook: superfície de ruptura crítica com nível d'água (FS = 1,357) e sem água (FS = 1,832).*

## Sumário

- [1. Como usar este manual](#1-como-usar-este-manual)
  - [1.1 O que o notebook faz](#11-o-que-o-notebook-faz)
  - [1.2 Como este manual está organizado](#12-como-este-manual-está-organizado)
  - [1.3 Convenções](#13-convenções)
- [2. Início rápido](#2-início-rápido)
- [3. Conceitos fundamentais](#3-conceitos-fundamentais)
  - [3.1 Fator de Segurança e equilíbrio limite](#31-fator-de-segurança-e-equilíbrio-limite)
  - [3.2 O método das fatias](#32-o-método-das-fatias)
  - [3.3 Nível d'água, poropressão e tensão efetiva](#33-nível-dágua-poropressão-e-tensão-efetiva)
  - [3.4 A fórmula de Bishop Simplificado em tensões efetivas](#34-a-fórmula-de-bishop-simplificado-em-tensões-efetivas)
  - [3.5 Por que o cálculo é iterativo](#35-por-que-o-cálculo-é-iterativo)
- [4. Preparando os dados de entrada](#4-preparando-os-dados-de-entrada)
  - [4.1 Geometria do talude (Seção 2)](#41-geometria-do-talude-seção-2)
  - [4.2 A orientação do talude importa? (crista à esquerda, pé à direita)](#42-a-orientação-do-talude-importa-crista-à-esquerda-pé-à-direita)
  - [4.3 Nível d'água (Seção 3)](#43-nível-dágua-seção-3)
  - [4.4 Solo e água (Seção 4)](#44-solo-e-água-seção-4)
- [5. O notebook seção por seção](#5-o-notebook-seção-por-seção)
- [6. As funções explicadas](#6-as-funções-explicadas)
  - [6.1 surface_y(x)](#61-surface_yx)
  - [6.2 water_y(x)  [nova]](#62-water_yx--nova)
  - [6.3 circle_y(x, xc, yc, R)](#63-circle_yx-xc-yc-r)
  - [6.4 find_entry_exit(xc, yc, R)](#64-find_entry_exitxc-yc-r)
  - [6.5 build_slices(xc, yc, R, n_slices, com_agua)  [alterada]](#65-build_slicesxc-yc-r-n_slices-com_agua--alterada)
  - [6.6 bishop_fs(slices, phi_rad, c)  [alterada]](#66-bishop_fsslices-phi_rad-c--alterada)
  - [6.7 fs_for_circle(xc, yc, R, n_slices, com_agua)  [alterada]](#67-fs_for_circlexc-yc-r-n_slices-com_agua--alterada)
  - [6.8 buscar_circulo_critico(com_agua)  [nova]](#68-buscar_circulo_criticocom_agua--nova)
- [7. Exemplo resolvido: uma fatia passo a passo](#7-exemplo-resolvido-uma-fatia-passo-a-passo)
- [8. Lendo os resultados](#8-lendo-os-resultados)
  - [8.1 Superfície de ruptura crítica (Seção 8.1)](#81-superfície-de-ruptura-crítica-seção-81)
  - [8.2 Poropressão e peso das fatias (Seção 8.2)](#82-poropressão-e-peso-das-fatias-seção-82)
  - [8.3 Mapa de FS por centro testado (Seção 8.3)](#83-mapa-de-fs-por-centro-testado-seção-83)
  - [8.4 Avaliar um círculo específico (Seção 9)](#84-avaliar-um-círculo-específico-seção-9)
  - [8.5 Resumo final (Seção 10)](#85-resumo-final-seção-10)
- [9. Quanto a posição do nível d'água importa?](#9-quanto-a-posição-do-nível-dágua-importa)
- [10. Solução de problemas](#10-solução-de-problemas)
- [11. Exercícios propostos](#11-exercícios-propostos)
- [12. Simplificações e limitações](#12-simplificações-e-limitações)
- [13. Considerações finais](#13-considerações-finais)
  - [Glossário de variáveis do código](#glossário-de-variáveis-do-código)

---

## 1. Como usar este manual

Este manual acompanha o notebook `analise_estabilidade_talude_bishop_nivel_agua.ipynb` e foi escrito para quem quer **usar** o código com seus próprios dados e **entender** o que ele calcula, sem precisar ler o código-fonte linha a linha. Ele substitui o documento de apoio da versão anterior (talude seco), incorporando o nível d'água.

### 1.1 O que o notebook faz

O notebook calcula o **Fator de Segurança (FS)** de um talude pelo método de equilíbrio limite de **Bishop Simplificado**, com superfícies de ruptura **circulares**, considerando um **nível d'água (linha freática)** dentro do maciço. Em resumo:

- você informa o **perfil do terreno**, a **linha do nível d'água** e as **propriedades do solo e da água**;
- o código testa milhares de círculos e encontra o **círculo crítico** — aquele com o menor FS;
- o cálculo é feito em **tensões efetivas**: abaixo do N.A. o solo é mais pesado (γ_sat) e a água empurra a base de cada fatia (poropressão u);
- para comparação, o mesmo cálculo é repetido **sem água** (talude seco), mostrando quanto a água reduz a segurança.

### 1.2 Como este manual está organizado

| Se você quer… | Vá para |
|---|---|
| rodar o notebook o quanto antes | Capítulo 2 — Início rápido |
| entender a teoria (FS, fatias, poropressão, Bishop) | Capítulo 3 — Conceitos fundamentais |
| preparar seus dados — **inclusive a orientação do talude** | Capítulo 4 — Dados de entrada (seção 4.2) |
| saber o que cada seção e cada função faz | Capítulos 5 e 6 |
| ver uma conta feita à mão | Capítulo 7 — Exemplo resolvido |
| interpretar gráficos e números | Capítulos 8 e 9 |
| resolver um erro | Capítulo 10 — Solução de problemas |
| praticar | Capítulo 11 — Exercícios |

### 1.3 Convenções

- **Unidades:** comprimentos em metros (m), pesos específicos em kN/m³, coesão e poropressão em kPa, ângulos em graus. Forças por metro de talude (kN/m), pois a análise é bidimensional.
- **Sinais:** γ é o peso específico natural (acima do N.A.), γ_sat o saturado (abaixo do N.A.), γ_w o da água; c' e φ' são a coesão e o ângulo de atrito **efetivos**.
- Trechos de código aparecem em `fonte monoespaçada`. Os números do exemplo vêm da execução do próprio notebook.
- Caixas coloridas destacam **dicas** (azul), **atenções** (laranja) — pontos em que é fácil errar — e **exemplos** (verde).

## 2. Início rápido

Se você só quer ver o notebook funcionando com o exemplo, siga os passos abaixo. Cada passo é detalhado nos capítulos seguintes.

1. **Abra o notebook** no Jupyter Lab, VS Code (extensão Jupyter) ou Google Colab.
2. **Confira as bibliotecas**: numpy, matplotlib e scipy. Se faltar alguma: `pip install numpy matplotlib scipy`.
3. **Seção 2 — geometria:** edite a lista `superficie` com os pontos (x, y) do perfil, em x crescente, com o talude **descendo da esquerda para a direita** (veja a seção 4.2).
4. **Seção 3 — nível d'água:** edite a lista `nivel_agua` (x crescente, abaixo do terreno) e confira o gráfico gerado.
5. **Seção 4 — solo e água:** ajuste `gamma`, `gamma_sat`, `gamma_w`, `c` e `phi_deg`.
6. **Execute tudo** (Run All). A busca da Seção 7 roda duas vezes (com e sem água) e pode levar de segundos a alguns minutos.
7. **Leia o resultado:** FS com N.A., FS seco e a redução percentual aparecem na Seção 7 e no resumo da Seção 10; os gráficos estão na Seção 8.

> **📝 Exemplo**
>
> Com os dados que já vêm no notebook (talude de 10 m de altura, inclinação 2H:1V, γ = 18, γ_sat = 20, c' = 10 kPa, φ' = 28°), o resultado esperado é:
>
> **FS com nível d'água = 1,357** · **FS seco = 1,832** · **redução de 26,0 %**. Se você obtiver esses valores, sua instalação está funcionando.

## 3. Conceitos fundamentais

### 3.1 Fator de Segurança e equilíbrio limite

A análise por **equilíbrio limite** compara, ao longo de uma superfície de ruptura imaginada, a resistência ao cisalhamento que o solo consegue oferecer com o esforço que o peso do solo provoca. Essa comparação é o Fator de Segurança:

$$
FS = \frac{\text{resistência disponível}}{\text{solicitação atuante}}
$$

Para superfícies circulares, as duas grandezas são **momentos** em torno do centro do círculo. A leitura é direta:

| Faixa de FS | Interpretação |
|---|---|
| **FS < 1,0** | Instável — segundo o modelo, a ruptura já teria ocorrido. |
| **1,0 ≤ FS < 1,5** | Estável, mas abaixo da margem de segurança usual de projeto. |
| **FS ≥ 1,5** | Estável, dentro da margem usual (a referência exata depende da norma e do tipo de obra). |

Como não se sabe de antemão onde a ruptura ocorreria, testam-se muitos círculos. **O FS do talude é o menor FS encontrado** — o do chamado círculo crítico.

### 3.2 O método das fatias

A massa de solo acima do círculo é dividida em **fatias verticais**. Para cada fatia *i* conhecemos a largura bᵢ, a altura hᵢ, o peso Wᵢ e o ângulo αᵢ que a base faz com a horizontal. Somando as contribuições de todas as fatias obtêm-se os momentos resistente e atuante. O notebook usa 20 fatias durante a busca (mais rápido) e 50 no refinamento e nos resultados finais (mais preciso).

### 3.3 Nível d'água, poropressão e tensão efetiva

Com água no maciço, cada fatia passa a ter duas partes — e a base recebe a pressão da água:

![anatomia de uma fatia atravessada pelo nível d'água.](manual_estabilidade_talude_nivel_agua_imagens/fatia.png)

*Figura 2 — anatomia de uma fatia atravessada pelo nível d'água.*

| Grandeza | Fórmula no código | Significado |
|---|---|---|
| altura saturada | `h_sat = máx(0, mín(y_NA, y_topo) − y_base)` | trecho da fatia abaixo do N.A. |
| altura não saturada | `h_nat = h − h_sat` | trecho acima do N.A. |
| peso da fatia | `W = b·(γ·h_nat + γ_sat·h_sat)` | abaixo do N.A. o solo pesa mais |
| poropressão na base | `u = γ_w·máx(0, y_NA − y_base)` | pressão hidrostática da água |

O ponto central é o **princípio da tensão efetiva** de Terzaghi: só a parcela da tensão transmitida de grão para grão — a tensão efetiva — mobiliza atrito.

$$
\sigma' = \sigma - u \qquad\qquad \tau = c' + \sigma'\,\tan\varphi'
$$

A água nos poros “sustenta” parte do peso da fatia. Na base, a força que efetivamente aperta um grão contra o outro deixa de ser W e passa a ser **W − u·b**. Menos força normal efetiva significa menos atrito, e portanto menos resistência.

> **💡 Dica**
>
> Analogia: empurrar uma caixa pesada que está apoiada sobre uma boia. Quanto mais a água sustenta a caixa, menos ela “agarra” no chão — e mais fácil ela desliza.

#### Os dois efeitos da água

- **O solo fica mais pesado** (γ_sat > γ): aumenta o momento atuante. Efeito desfavorável, mas pequeno.
- **A poropressão reduz a força normal efetiva** (W − u·b): diminui o atrito. É o efeito dominante.

No círculo crítico do exemplo, o peso total das 50 fatias é ΣW = 1.824 kN/m e o empuxo da água na base é Σu·b = 580 kN/m: cerca de **32 % do peso deixa de gerar atrito**.

### 3.4 A fórmula de Bishop Simplificado em tensões efetivas

Bishop (1955) admite que as forças **cisalhantes** entre fatias vizinhas se anulam e impõe o equilíbrio vertical de cada fatia e o equilíbrio de momentos em torno do centro do círculo. O resultado, implementado no notebook, é:

$$
FS = \frac{\displaystyle\sum_i \frac{c'\,b_i + (W_i - u_i\,b_i)\tan\varphi'}{m_{\alpha,i}}}{\displaystyle\sum_i W_i \sin\alpha_i}
$$

$$
m_{\alpha,i} = \cos\alpha_i + \frac{\sin\alpha_i\,\tan\varphi'}{FS}
$$

| Termo | O que representa |
|---|---|
| `c'·bᵢ` | parcela de coesão ao longo da base da fatia |
| `(Wᵢ − uᵢ·bᵢ)·tan φ'` | **parcela de atrito com o peso efetivo** — é aqui que a água entra (termo novo em relação ao talude seco) |
| `mα,i` | fator geométrico que depende do ângulo da base e do próprio FS |
| `Σ Wᵢ·sin αᵢ` | momento motor do peso (dividido por R); W já inclui γ_sat abaixo do N.A. |

> **💡 Dica**
>
> Se u = 0 e γ_sat = γ, a fórmula volta exatamente à do talude seco. É assim que o notebook calcula o caso de comparação (parâmetro `com_agua=False`).

### 3.5 Por que o cálculo é iterativo

O FS aparece dos **dois lados** da equação (no resultado e dentro de mα,i). A solução é por **iteração de ponto fixo**: chuta-se FS = 1,0, calcula-se mα,i, obtém-se um novo FS, e repete-se até que a diferença entre duas iterações seja menor que 10⁻⁵. No círculo crítico do exemplo:

| Iteração | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| FS | 1,0000 | 1,2955 | 1,3481 | 1,3557 | 1,3567 | 1,3569 | 1,3569 | 1,3569 |

Bastaram 7 iterações. Na prática, a convergência é sempre rápida; o notebook limita a 100 iterações por segurança.

#### Instabilidade numérica de mα

Em fatias com base muito inclinada (α próximo de ±90°) e FS baixo, mα pode se aproximar de zero ou ficar negativo, gerando divisões absurdas. O notebook **rejeita** qualquer círculo em que algum mα fique abaixo de 0,2. Isso é um cuidado clássico em implementações do método de Bishop.

## 4. Preparando os dados de entrada

Todos os dados são informados em três células: Seção 2 (geometria), Seção 3 (nível d'água) e Seção 4 (solo e água). As regras abaixo evitam a maior parte dos erros.

### 4.1 Geometria do talude (Seção 2)

*Entrada da Seção 2 (exemplo do notebook)*

```python
superficie = [
    (-15.0, 10.0),   # início do platô superior (crista)
    (  0.0, 10.0),   # crista do talude
    ( 20.0,  0.0),   # pé do talude (inclinação 2H:1V)
    ( 40.0,  0.0),   # trecho horizontal inferior
]
```

- Pontos (x, y) em metros, com **x estritamente crescente** — o código verifica isso e para com erro se não for.
- O terreno entre pontos é interpolado em linha reta: use quantos pontos forem necessários para descrever o perfil.
- Inclua **platôs antes da crista e depois do pé** (no exemplo, 15 m e 20 m). Sem esse espaço, os círculos não têm onde entrar e sair do terreno.
- O talude deve **descer da esquerda para a direita**. Isso não é detalhe — veja a próxima seção.

### 4.2 A orientação do talude importa? (crista à esquerda, pé à direita)

> **✅ Resposta curta**
>
> **Sim, neste notebook a orientação importa.** O talude precisa ser informado **descendo da esquerda para a direita**: crista nos valores menores de x e pé nos valores maiores.
>
> A física não tem lado preferido — um talude e a sua imagem no espelho têm exatamente o mesmo FS. A restrição vem do **código**, que foi escrito assumindo essa orientação. Se o seu talude desce para a esquerda, a solução é simples: **espelhar as coordenadas** antes de rodar (receita na seção 4.2.3).

![A) orientação aceita pelo código; B) talude que desce para a esquerda, que não funciona diretamente; C) o mesmo talude B após espelhar as coordenadas, que volta a funcionar.](manual_estabilidade_talude_nivel_agua_imagens/orientacao.png)

*Figura 3 — A) orientação aceita pelo código; B) talude que desce para a esquerda, que não funciona diretamente; C) o mesmo talude B após espelhar as coordenadas, que volta a funcionar.*

#### 4.2.1 Por que o código exige essa orientação

Há **três pontos independentes** do notebook que assumem crista à esquerda e pé à direita. Basta um deles falhar para a análise não funcionar:

| Onde | O que o código assume | O que acontece com o talude invertido |
|---|---|---|
| **Seção 3.1** — detecção de crista e pé | Crista = **último** ponto de maior cota; pé = ponto mais baixo **depois** (à direita) da crista: `i_toe = i_crest + argmin(surf_y[i_crest:])`. | A crista é o último ponto do perfil; não há pontos depois dela, então **pé = crista**. Altura e largura da face dão 0 e os intervalos de busca da Seção 7 ficam degenerados. |
| **find_entry_exit** — filtro de plausibilidade | A entrada do círculo (menor x) fica perto da crista e a saída (maior x) perto do pé: rejeita se `x_entry > x_crest + margin` ou `x_exit < x_toe − margin`. | Com crista e pé trocados de lado, círculos fisicamente corretos são descartados. |
| **build_slices / bishop_fs** — sinal de α | Convenção `sin α = (xc − x)/R`, que torna o momento motor ΣW·sin α **positivo** quando a massa escorrega no sentido de +x (para a direita). | A massa escorrega para −x, ΣW·sin α fica **negativo** e `bishop_fs` devolve NaN (“momento motor ≤ 0”). |

#### 4.2.2 Comprovação: o exemplo espelhado

Para confirmar, o notebook foi executado com o **mesmo talude do exemplo espelhado** — crista em x = 0 a 15 m, pé em x = −20 m, descendo para a esquerda, e o N.A. espelhado da mesma forma. O resultado:

*Saída do notebook com o talude espelhado (sem correção)*

```text
Crista: x=15.00 m, y=10.00 m
Pe do talude: x=15.00 m, y=10.00 m
Altura do talude: 0.00 m | Largura da face: 0.00 m
...
RuntimeError: Nenhum círculo válido foi encontrado nos intervalos de busca.
```

Mesmo corrigindo à mão a crista e o pé, o círculo crítico espelhado (centro em xc = −14,47 m) dá **ΣW·sin α = −652 kN/m**, exatamente o oposto dos +652 kN/m do caso original — e o FS sai NaN. Ou seja: os três pontos da tabela acima realmente falham.

#### 4.2.3 Receita: como analisar um talude que desce para a esquerda

Espelhe as coordenadas trocando x por −x e **invertendo a ordem da lista** (para que x continue crescente). Acrescente uma linha logo depois de cada lista:

*Linhas a acrescentar nas Seções 2 e 3*

```python
# Seção 2 — logo depois da lista "superficie":
# talude desce para a esquerda? espelhe x -> -x e inverta a ordem
superficie = [(-x, y) for (x, y) in reversed(superficie)]

# Seção 3 — logo depois da lista "nivel_agua" (espelhe também o N.A.!):
nivel_agua = [(-x, y) for (x, y) in reversed(nivel_agua)]
```

> **⚠️ Atenção**
>
> Espelhe **as duas listas** — terreno e nível d'água. Espelhar só o terreno deixa o N.A. do lado errado e o resultado não tem significado.
>
> As linhas devem ficar **antes** de `surf_x = …` e `wt_x = …`, que estão na mesma célula logo abaixo de cada lista.

**Como ler os resultados depois de espelhar:**

- FS, raio R e cota do centro yc **não mudam** — use-os diretamente.
- Todas as coordenadas x mudam de sinal. Para voltar ao seu sistema original: **xc_real = −xc_crit** (e o mesmo para entrada, saída e qualquer x lido nos gráficos).
- Os gráficos aparecem espelhados (o talude descendo para a direita). Isso é só a forma de desenho; a análise é a mesma.

> **📝 Exemplo**
>
> Um perfil medido como `[(0, 0), (20, 0), (40, 10), (55, 10)]` (pé à esquerda, crista à direita) vira, após o espelhamento, `[(-55, 10), (-40, 10), (-20, 0), (0, 0)]`. Se o notebook encontrar xc = −24,5 m, no sistema original o centro está em x = +24,5 m.

#### 4.2.4 E se o terreno tiver duas faces (um aterro ou um morro)?

O código analisa **apenas uma face**: a que desce para a direita a partir do último ponto mais alto. Para um aterro com as duas faces, faça duas análises: uma com o perfil como está (face direita) e outra com o perfil espelhado (face esquerda). O FS do aterro é o menor dos dois.

> **💡 Dica**
>
> Uma melhoria possível do notebook é detectar automaticamente o sentido de descida (comparando a cota do primeiro e do último ponto) e espelhar sozinho os dados, devolvendo os resultados já no sistema original. Até lá, a receita acima resolve.

### 4.3 Nível d'água (Seção 3)

*Entrada da Seção 3 (exemplo do notebook)*

```python
nivel_agua = [
    (-15.0,  8.0),
    (  0.0,  7.5),
    ( 20.0, -0.5),
    ( 40.0, -0.5),
]
```

- Pontos (x, y) em metros, **x estritamente crescente**, mesma orientação do terreno.
- O N.A. deve ficar **abaixo do terreno ou coincidente** com ele. Se algum trecho ficar acima, o código **avisa** e limita o N.A. à cota do terreno (não há lâmina d'água livre sobre o talude).
- Fora do intervalo informado, o N.A. é **estendido na horizontal** (mesma cota do primeiro e do último ponto).
- Use o **gráfico de conferência** gerado pela célula para checar se a linha azul ficou onde você queria.

![gráfico de conferência da Seção 3: perfil do talude, nível d'água e região saturada.](manual_estabilidade_talude_nivel_agua_imagens/cell6.png)

*Figura 4 — gráfico de conferência da Seção 3: perfil do talude, nível d'água e região saturada.*

> **💡 Dica**
>
> Para simular o talude seco basta colocar o N.A. bem abaixo do pé (por exemplo, 10 m abaixo de toda a superfície). Mas não é preciso: o notebook já calcula o caso seco automaticamente.

### 4.4 Solo e água (Seção 4)

| Variável | Símbolo | Unidade | Exemplo | Observação |
|---|---|---|---|---|
| `gamma` | γ | kN/m³ | 18,0 | peso específico natural, acima do N.A. |
| `gamma_sat` | γ_sat | kN/m³ | 20,0 | peso específico saturado, abaixo do N.A.; normalmente maior que γ |
| `gamma_w` | γ_w | kN/m³ | 9,81 | peso específico da água (há quem use 10,0) |
| `c` | c' | kPa | 10,0 | coesão **efetiva** |
| `phi_deg` | φ' | graus | 28,0 | ângulo de atrito **efetivo** |

> **⚠️ Atenção**
>
> Como a análise é em tensões efetivas, c' e φ' devem vir de ensaios **drenados** (ou não drenados com medição de poropressão). Usar parâmetros totais/não drenados junto com poropressão conta o efeito da água duas vezes.

## 5. O notebook seção por seção

O notebook tem 10 seções e deve ser executado **em ordem**: cada célula usa variáveis criadas nas anteriores. A tabela indica o que mudou em relação ao notebook original do talude seco.

| Seção | Conteúdo | Você edita? | Em relação ao original |
|---|---|---|---|
| 1. Bibliotecas | numpy, matplotlib e `brentq`/`minimize` do scipy | não | igual |
| 2. Geometria | lista `superficie` e verificação de x crescente | **sim** | igual |
| 3. Nível d'água | lista `nivel_agua`, aviso de N.A. acima do terreno, gráfico | **sim** | **nova** |
| 3.1 Crista e pé | detecção automática para dimensionar a busca | raramente | igual |
| 4. Solo e água | γ, γ_sat, γ_w, c', φ' | **sim** | alterada |
| 5. Funções de geometria | `surface_y`, `water_y`, `circle_y`, `find_entry_exit` | não | + `water_y` |
| 6. Fatias e FS | `build_slices`, `bishop_fs`, `fs_for_circle` em tensões efetivas | não | alterada |
| 7. Busca do círculo | grid + Nelder-Mead, com e sem N.A. | às vezes (intervalos) | alterada |
| 8. Gráficos | ruptura crítica, poropressão e pesos, mapa de FS | não | alterada |
| 9. Círculo manual | FS de um círculo escolhido, com e sem N.A. | opcional | alterada |
| 10. Resumo | resultados, comparação e veredito | não | alterada |

> **💡 Dica**
>
> Depois de mudar qualquer dado de entrada, rode novamente **todas** as células a partir dela (ou simplesmente Run All). Rodar só a última célula mostra resultados antigos.

## 6. As funções explicadas

O cálculo inteiro está em oito funções pequenas, cada uma com uma única tarefa. O código abaixo é o do notebook, sem alterações. As marcas **[nova]** e **[alterada]** indicam o que mudou em relação ao notebook do talude seco.

### 6.1 surface_y(x)

Devolve a cota do terreno em qualquer x, por interpolação linear entre os pontos de `superficie`. Aceita um número ou um array do NumPy.

```python
def surface_y(x):
    """Altura da superfície do talude em x (interpolação linear)."""
    return np.interp(x, surf_x, surf_y)
```

### 6.2 water_y(x)  [nova]

Devolve a cota do nível d'água em x, interpolando `nivel_agua`. O `np.minimum` com `surface_y(x)` garante que o N.A. **nunca fica acima do terreno**. Fora do intervalo informado, `np.interp` repete a cota da extremidade (extensão horizontal).

```python
def water_y(x):
    """Cota do nivel d'agua em x (interpolacao linear), limitada a cota do terreno."""
    return np.minimum(np.interp(x, wt_x, wt_y), surface_y(x))
```

### 6.3 circle_y(x, xc, yc, R)

Devolve a cota do **arco inferior** do círculo de centro (xc, yc) e raio R, isto é, y = yc − √(R² − (x − xc)²). Fora do domínio do círculo (|x − xc| > R) devolve NaN, sinalizando “não definido aqui” em vez de gerar erro.

```python
def circle_y(x, xc, yc, R):
    """Altura do arco inferior do círculo (centro xc,yc raio R) em x. NaN fora do domínio do círculo."""
    val = R**2 - (x - xc)**2
    with np.errstate(invalid="ignore"):
        return np.where(val >= 0, yc - np.sqrt(np.maximum(val, 0)), np.nan)
```

### 6.4 find_entry_exit(xc, yc, R)

Encontra onde o círculo corta o terreno: a **entrada** (perto da crista) e a **saída** (perto do pé).

- amostra 400 pontos no trecho comum ao perfil e ao círculo e calcula `diff = terreno − círculo` (positivo onde há solo acima do círculo);
- procura **mudanças de sinal** de `diff` e refina cada cruzamento com o método de Brent (`brentq`);
- a menor raiz é a entrada e a maior, a saída;
- aplica o **filtro de plausibilidade**: o círculo deve atravessar praticamente toda a face, da crista ao pé (com folga `margin`). É um dos pontos que exigem a orientação esquerda → direita (seção 4.2);
- devolve `None` se o círculo não gerar uma superfície de ruptura válida.

```python
def find_entry_exit(xc, yc, R, n_samples=400):
    """
    Encontra as coordenadas x de entrada e saída do círculo de ruptura na
    superfície do talude (pontos onde superficie(x) == circulo(x)).
    Retorna (x_entrada, x_saida) ou None se o círculo não gerar uma
    superfície de ruptura válida.
    """
    dom_min = max(surf_x[0], xc - R)
    dom_max = min(surf_x[-1], xc + R)
    if dom_max - dom_min < 1e-6:
        return None

    xs = np.linspace(dom_min, dom_max, n_samples)
    diff = surface_y(xs) - circle_y(xs, xc, yc, R)

    # precisa haver trecho onde a superficie esta acima do circulo (massa de solo)
    if not np.any(diff > 0):
        return None

    # localizar mudancas de sinal e refinar com brentq
    roots = []
    for i in range(len(xs) - 1):
        d0, d1 = diff[i], diff[i + 1]
        if np.isnan(d0) or np.isnan(d1):
            continue
        if d0 == 0:
            roots.append(xs[i])
        elif d0 * d1 < 0:
            f = lambda x: surface_y(x) - circle_y(x, xc, yc, R)
            try:
                r = brentq(f, xs[i], xs[i + 1])
                roots.append(r)
            except ValueError:
                continue

    if len(roots) < 2:
        return None

    x_entry, x_exit = min(roots), max(roots)
    if x_exit - x_entry < 1e-3:
        return None

    # verificar que ha material entre entrada e saida (circulo abaixo da superficie)
    x_mid = 0.5 * (x_entry + x_exit)
    if surface_y(x_mid) - circle_y(x_mid, xc, yc, R) <= 0:
        return None

    # a superficie de ruptura deve atravessar essencialmente toda a face do
    # talude (da crista ao pe, com alguma folga), evitando "fatias" rasas e
    # sem sentido fisico que tocam apenas o platô superior ou inferior
    if x_entry > x_crest + margin or x_exit < x_toe - margin:
        return None

    return x_entry, x_exit
```

### 6.5 build_slices(xc, yc, R, n_slices, com_agua)  [alterada]

Divide a massa entre entrada e saída em fatias de mesma largura e calcula, no ponto médio de cada uma:

- `b`, `h` (= y_topo − y_base) e `alpha`, com a convenção `sin α = (xc − x)/R` (seção 4.2);
- com `com_agua=True`: `h_sat`, `h_nat`, o peso `W = b·(γ·h_nat + γ_sat·h_sat)` e a poropressão `u = γ_w·(y_NA − y_base)` (zero se a base estiver acima do N.A.);
- com `com_agua=False`: `h_sat = 0`, `u = 0` e `W = γ·b·h` — o talude seco;
- devolve um dicionário com todos os arrays, ou `None` se alguma fatia tiver altura ≤ 0.

```python
def build_slices(xc, yc, R, n_slices=30, com_agua=True):
    """
    Constrói as fatias entre a entrada e a saída do círculo de ruptura.
    Se com_agua=True, considera o nível d'água (peso saturado e poropressão).
    Retorna dict com arrays b, h, h_sat, alpha, W, u, x_mid — ou None se inválido.
    """
    entry_exit = find_entry_exit(xc, yc, R)
    if entry_exit is None:
        return None
    x_entry, x_exit = entry_exit

    edges = np.linspace(x_entry, x_exit, n_slices + 1)
    x_left, x_right = edges[:-1], edges[1:]
    x_mid = 0.5 * (x_left + x_right)
    b = x_right - x_left

    y_top = surface_y(x_mid)
    y_bot = circle_y(x_mid, xc, yc, R)

    if np.any(np.isnan(y_bot)):
        return None

    h = y_top - y_bot
    if np.any(h <= 0):
        return None

    # convencao: sin(alpha) = (xc - x)/R, positivo do lado da saida (pe) do talude
    arg = np.clip((xc - x_mid) / R, -1.0, 1.0)
    alpha = np.arcsin(arg)  # rad

    if com_agua:
        y_wt = water_y(x_mid)
        h_sat = np.clip(np.minimum(y_wt, y_top) - y_bot, 0.0, None)  # parcela saturada
        h_nat = h - h_sat                                             # parcela acima do N.A.
        W = b * (gamma * h_nat + gamma_sat * h_sat)
        u = gamma_w * np.clip(y_wt - y_bot, 0.0, None)               # poropressao na base
    else:
        h_sat = np.zeros_like(h)
        W = gamma * b * h
        u = np.zeros_like(h)

    return {"b": b, "h": h, "h_sat": h_sat, "alpha": alpha, "W": W, "u": u,
            "x_mid": x_mid, "x_entry": x_entry, "x_exit": x_exit}
```

### 6.6 bishop_fs(slices, phi_rad, c)  [alterada]

Aplica a iteração de Bishop (seção 3.5). A única linha nova em relação ao talude seco é o peso efetivo:

- `N_ef = np.clip(W − u·b, 0, None)` substitui W no termo de atrito. O `clip` impede um “atrito negativo” caso a poropressão seja muito alta;
- rejeita o círculo se o momento motor Σ W·sin α for ≤ 0 ou se algum mα ficar abaixo de 0,2;
- devolve o par (FS, convergiu).

```python
def bishop_fs(slices, phi_rad, c, fs_init=1.0, tol=1e-5, max_iter=100):
    """
    Calcula o FS pelo método de Bishop Simplificado (tensões efetivas, iterativo)
    para um conjunto de fatias. Retorna (FS, convergiu:bool).
    """
    b, alpha, W, u = slices["b"], slices["alpha"], slices["W"], slices["u"]
    tan_phi = np.tan(phi_rad)

    denom = np.sum(W * np.sin(alpha))
    if denom <= 0:
        return np.nan, False

    # parcela de atrito com o peso efetivo (W - u*b); nao pode ser negativa
    N_ef = np.clip(W - u * b, 0.0, None)

    fs = fs_init
    for _ in range(max_iter):
        m_alpha = np.cos(alpha) + np.sin(alpha) * tan_phi / fs
        if np.any(m_alpha < 0.2):  # instabilidade numerica tipica do metodo
            return np.nan, False
        numer = np.sum((c * b + N_ef * tan_phi) / m_alpha)
        fs_new = numer / denom
        if abs(fs_new - fs) < tol:
            return fs_new, True
        fs = fs_new
    return fs, False
```

### 6.7 fs_for_circle(xc, yc, R, n_slices, com_agua)  [alterada]

Atalho que encadeia `build_slices` e `bishop_fs`: recebe só o círculo e devolve o FS, ou NaN se o círculo for inválido ou não convergir. É chamada milhares de vezes durante a busca.

```python
def fs_for_circle(xc, yc, R, n_slices=30, com_agua=True):
    """Calcula o FS (Bishop) para um círculo (xc, yc, R). NaN se inválido."""
    slices = build_slices(xc, yc, R, n_slices=n_slices, com_agua=com_agua)
    if slices is None:
        return np.nan
    fs, converged = bishop_fs(slices, phi_rad, c)
    if not converged:
        return np.nan
    return fs
```

### 6.8 buscar_circulo_critico(com_agua)  [nova]

Reúne a busca do círculo crítico numa função, para poder executá-la duas vezes (com e sem N.A.). Funciona em duas etapas:

- **grid search**: testa todas as combinações de `xc_vals`, `yc_vals` e `R_vals` (22 × 16 × 22 = 7.744 círculos, com 20 fatias) e guarda o de menor FS. No exemplo com N.A., 211 círculos são válidos e o melhor tem FS = 1,391;
- **refinamento Nelder-Mead** (`scipy.optimize.minimize`): parte do melhor círculo do grid e ajusta (xc, yc, R) continuamente, com 50 fatias. Círculos inválidos recebem penalidade 1000. Resultado do exemplo: FS = 1,357;
- devolve xc, yc, R, FS, a lista de resultados do grid e o melhor FS do grid.

```python
def buscar_circulo_critico(com_agua=True):
    """Grid search + refinamento Nelder-Mead. Retorna (xc, yc, R, FS, resultados_grid)."""
    melhor_fs, melhor_circ, resultados = np.inf, None, []
    for xc in xc_vals:
        for yc in yc_vals:
            for R in R_vals:
                fs = fs_for_circle(xc, yc, R, n_slices=n_slices_grid, com_agua=com_agua)
                if np.isfinite(fs) and fs > 0:
                    resultados.append((xc, yc, R, fs))
                    if fs < melhor_fs:
                        melhor_fs, melhor_circ = fs, (xc, yc, R)

    if melhor_circ is None:
        raise RuntimeError("Nenhum círculo válido foi encontrado nos intervalos de busca. "
                           "Ajuste xc_vals / yc_vals / R_vals.")

    def objetivo(params):
        fs = fs_for_circle(*params, n_slices=50, com_agua=com_agua)
        return fs if np.isfinite(fs) else 1e3  # penaliza círculos inválidos

    res = minimize(objetivo, x0=np.array(melhor_circ), method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1e-5, "maxiter": 500})
    xc_c, yc_c, R_c = res.x
    return xc_c, yc_c, R_c, objetivo(res.x), resultados, melhor_fs
```

#### Intervalos de busca do grid

*Início da Seção 7: os intervalos são calculados a partir da crista, do pé e da diagonal do talude*

```python
# --- Intervalos de busca (ajuste se necessário) ---
n_xc, n_yc, n_R = 22, 16, 22   # resolução do grid

xc_vals = np.linspace(x_crest - 0.5 * slope_width, x_toe + 1.5 * slope_width, n_xc)
yc_vals = np.linspace(y_crest + 0.1 * diag, y_crest + 2.5 * diag, n_yc)
R_vals = np.linspace(0.5 * diag, 3.0 * diag, n_R)

n_slices_grid = 20  # menos fatias no grid (mais rápido); refina depois
```

Se o círculo crítico ficar na **borda** do mapa de FS (seção 8.3), amplie esses intervalos. Para uma busca mais fina, aumente `n_xc`, `n_yc` e `n_R` — o tempo de execução cresce na mesma proporção do número de círculos.

## 7. Exemplo resolvido: uma fatia passo a passo

Para ver os números por dentro, tomemos a **fatia 26 de 50** do círculo crítico com N.A. (xc = 14,47 m, yc = 18,64 m, R = 19,51 m).

| Dado | Valor | Origem |
|---|---|---|
| x (meio da fatia) | 8,84 m | divisão entre entrada e saída |
| b | 0,465 m | largura de cada uma das 50 fatias |
| y_topo | 5,58 m | `surface_y(8,84)` |
| y_NA | 3,97 m | `water_y(8,84)` |
| y_base | −0,04 m | `circle_y(8,84, xc, yc, R)` |
| α | 16,8° | `arcsin((xc − x)/R)` |

| Passo | Conta | Resultado |
|---|---|---|
| 1. Altura total | h = 5,58 − (−0,04) | **5,62 m** |
| 2. Altura saturada | h_sat = 3,97 − (−0,04) | **4,00 m** |
| 3. Altura não saturada | h_nat = 5,62 − 4,00 | **1,62 m** |
| 4. Peso | W = 0,465 · (18 · 1,62 + 20 · 4,00) | **50,8 kN/m** |
| 5. Poropressão | u = 9,81 · 4,00 | **39,3 kPa** |
| 6. Empuxo na base | u · b = 39,3 · 0,465 | **18,3 kN/m** |
| 7. Peso efetivo | W − u·b = 50,8 − 18,3 | **32,5 kN/m** |
| 8. Capacidade de atrito | (W − u·b) · tan 28° = 32,5 · 0,532 | **17,3 kN/m** |

> **📝 Exemplo**
>
> **Comparação com a mesma fatia seca:** W = 18 · 5,62 · 0,465 = 47,1 kN/m e a capacidade de atrito seria 47,1 · 0,532 = **25,0 kN/m**.
>
> Com água, a fatia ficou **mais pesada** (50,8 contra 47,1 kN/m), mas sua capacidade de atrito **caiu 31 %** (17,3 contra 25,0 kN/m). É exatamente o mecanismo descrito na seção 3.3.

## 8. Lendo os resultados

### 8.1 Superfície de ruptura crítica (Seção 8.1)

![talude com N.A., região saturada e círculos críticos com e sem água.](manual_estabilidade_talude_nivel_agua_imagens/cell19.png)

*Figura 5 — talude com N.A., região saturada e círculos críticos com e sem água.*

- **Faixa azul e linha tracejada:** região saturada e nível d'água.
- **Arco vermelho e “+”:** círculo crítico com N.A. (FS = 1,357) e seu centro.
- **Pontilhado preto e “×”:** círculo crítico do talude seco (FS = 1,832) e seu centro.
- **Linhas cinza:** as 50 fatias do círculo crítico.

Observe que, com água, o círculo crítico fica **menor** (R = 19,5 m contra 23,7 m) e com o centro mais baixo: a ruptura “procura” a zona saturada, onde a resistência é menor. Confira sempre se a ruptura parte de perto da crista e termina perto do pé — se não, revise os dados ou os intervalos de busca.

### 8.2 Poropressão e peso das fatias (Seção 8.2)

![em cima, poropressão na base de cada fatia; embaixo, peso total W e empuxo u·b.](manual_estabilidade_talude_nivel_agua_imagens/cell21.png)

*Figura 6 — em cima, poropressão na base de cada fatia; embaixo, peso total W e empuxo u·b.*

- No gráfico de cima, u acompanha a forma do arco: é maior onde a base da fatia está mais funda abaixo do N.A. Máximo no exemplo: **40,1 kPa**.
- No de baixo, a parte azul de cada barra é quanto do peso a água “tira” do atrito. Nas fatias centrais é cerca de **um terço** do peso.
- No exemplo, **45 das 50 fatias** têm a base abaixo do N.A.; as poucas fatias sem poropressão ficam junto à crista.

### 8.3 Mapa de FS por centro testado (Seção 8.3)

![cada ponto é um centro (xc, yc) do grid, colorido pelo menor FS entre os raios testados.](manual_estabilidade_talude_nivel_agua_imagens/cell23.png)

*Figura 7 — cada ponto é um centro (xc, yc) do grid, colorido pelo menor FS entre os raios testados.*

- Do vermelho (FS baixo) ao verde (FS alto). A **estrela azul** é o centro refinado pelo Nelder-Mead.
- **Bom sinal:** uma região vermelha concentrada — a busca encontrou um mínimo bem definido.
- **Sinal de alerta:** o mínimo encostado na borda do mapa — amplie `xc_vals`, `yc_vals` ou `R_vals` na Seção 7.

### 8.4 Avaliar um círculo específico (Seção 9)

Defina `xc_manual`, `yc_manual` e `R_manual` para calcular o FS de um círculo qualquer — por exemplo, uma superfície de ruptura observada em campo — com e sem água. Por padrão, a célula usa o círculo crítico com N.A., o que permite isolar o efeito da água **na mesma superfície**:

| Caso | Círculo (xc; yc; R) em m | FS |
|---|---|---|
| Seco — círculo crítico do caso seco | 17,39; 23,53; 23,67 | 1,832 |
| Seco — círculo crítico do caso com N.A. | 14,47; 18,64; 19,51 | 1,905 |
| Com N.A. — círculo crítico | 14,47; 18,64; 19,51 | **1,357** |

Na mesma superfície, a água sozinha reduz o FS em 29 % (1,905 → 1,357). E a água também **muda qual é o círculo crítico** — por isso a busca precisa ser refeita para cada cenário, e não apenas reavaliada no círculo do caso seco.

### 8.5 Resumo final (Seção 10)

*Saída da Seção 10 com os dados do exemplo*

```text
=======================================================
RESUMO DA ANALISE DE ESTABILIDADE DE TALUDE COM N.A.
=======================================================
Metodo: Bishop Simplificado (superficie circular, tensoes efetivas)
Solo: gamma = 18.0 | gamma_sat = 20.0 | gamma_w = 9.81 kN/m3
      c' = 10.0 kPa | phi' = 28.0 graus
-------------------------------------------------------
Circulo critico COM nivel d'agua:
  xc = 14.470 m | yc = 18.640 m | R = 19.510 m
  FS minimo = 1.357
Circulo critico SEM nivel d'agua (referencia):
  xc = 17.391 m | yc = 23.527 m | R = 23.671 m
  FS minimo = 1.832
Reducao do FS pelo N.A.: 26.0 %
-------------------------------------------------------
=> Talude com FS abaixo do minimo usual de projeto (1.5)
=======================================================
```

> **⚠️ Atenção**
>
> O limiar de 1,5 é uma referência comum, mas o FS mínimo aceitável depende da norma, do tipo de obra, da qualidade da investigação geotécnica e das consequências de uma ruptura. O notebook não substitui o julgamento de um engenheiro geotécnico responsável.

## 9. Quanto a posição do nível d'água importa?

Para mostrar a sensibilidade do resultado, o código do notebook foi executado deslocando todo o nível d'água do exemplo para cima ou para baixo (somando Δ às cotas de `nivel_agua`, sempre limitado ao terreno):

| Deslocamento do N.A. | FS crítico | Leitura |
|---|---|---|
| −6 m | 1,832 | N.A. abaixo da superfície de ruptura: igual ao talude seco |
| −4 m | 1,832 | ainda abaixo da ruptura: igual ao talude seco |
| −2 m | 1,676 | volta a ficar acima de 1,5 |
| **0 m (exemplo)** | **1,357** | abaixo do mínimo usual |
| +1 m | 1,172 | N.A. aflorando no pé |
| +3 m | 1,024 | maciço praticamente todo saturado: à beira da ruptura |

- Enquanto o N.A. está **abaixo da superfície de ruptura**, ele não afeta o FS.
- A partir daí, **cada metro conta**: entre −2 m e 0 m o FS cai de 1,68 para 1,36.
- Por isso **drenagem** (rebaixar o N.A.) é uma das medidas de estabilização mais eficientes: neste exemplo, rebaixar 2 m traz o FS de volta para cima de 1,5.

> **💡 Dica**
>
> Esses valores não são impressos pelo notebook — foram obtidos rodando o mesmo código com o N.A. deslocado. Você pode reproduzi-los no exercício 1 do capítulo 11.

## 10. Solução de problemas

| Sintoma | Causa provável | O que fazer |
|---|---|---|
| `AssertionError: As coordenadas devem ter x estritamente crescente` | pontos fora de ordem ou com x repetido (inclusive degrau vertical) | ordene os pontos por x; para um degrau, use dois x muito próximos (ex.: 10,00 e 10,01) |
| `RuntimeError: Nenhum círculo válido foi encontrado…` | talude descendo para a esquerda; platôs curtos; intervalos de busca inadequados | confira a orientação (seção 4.2); alongue os platôs; amplie `xc_vals`, `yc_vals`, `R_vals` |
| Seção 3.1 mostra **altura 0** ou crista = pé | talude invertido (crista à direita) | espelhe terreno **e** N.A. (receita da seção 4.2.3) |
| Aviso: “o N.A. informado fica acima do terreno…” | algum ponto do N.A. acima da superfície | corrija os pontos; se for intencional (N.A. aflorante), o código limita o N.A. ao terreno |
| FS seco e com N.A. **iguais** | N.A. abaixo da superfície de ruptura crítica | resultado correto: a água não atinge a massa que escorrega |
| Estrela azul na **borda** do mapa de FS | mínimo fora da região pesquisada | amplie os intervalos de busca na Seção 7 |
| Ruptura não passa perto da crista ou do pé | geometria complexa (bermas, vários trechos) | ajuste manualmente `x_crest`, `x_toe` ou os intervalos de busca |
| `ImportError` ao importar o scipy | instalação do scipy danificada no ambiente Python | reinstale (`pip install --force-reinstall scipy`) ou use outro ambiente (ex.: Colab) |
| Execução muito lenta | grid muito fino | reduza `n_xc`, `n_yc`, `n_R` para testar; aumente só na análise final |

## 11. Exercícios propostos

Mude **uma coisa de cada vez** e anote, para cada caso, o FS com N.A., o FS seco e o círculo crítico.

1. **Rebaixe o N.A. em 2 m** (subtraia 2 de cada y de `nivel_agua`). O FS passa de 1,5? *Resposta esperada: FS ≈ 1,68.*
2. **Sature o talude:** coloque o N.A. coincidindo com o terreno. Quão perto de FS = 1 você chega? O que acontece com o círculo crítico?
3. **Areia sem coesão:** faça `c = 0`. Como mudam o FS e a profundidade do círculo crítico? Por quê?
4. **Talude mais suave:** mude o pé para (30; 0), ou seja, inclinação 3H:1V (ajuste também o N.A.). Quanto o FS com N.A. melhora?
5. **Orientação:** monte o exemplo descendo para a esquerda, confirme que o notebook falha, aplique a receita da seção 4.2.3 e verifique que obtém FS = 1,357 e que, convertido para o sistema original, o centro fica em xc = −14,47 m.
6. **Efeito de γ_sat:** repita a análise com `gamma_sat = gamma = 18`. A diferença no FS é grande ou pequena? O que isso diz sobre qual dos dois efeitos da água domina?

## 12. Simplificações e limitações

- **Poropressão hidrostática:** u é medida na vertical a partir do N.A. Isso despreza a curvatura das equipotenciais de um fluxo real; para percolação importante, use uma rede de fluxo ou análise de percolação.
- **Sem sucção:** acima do N.A., u = 0. O ganho de resistência do solo não saturado é ignorado (a favor da segurança).
- **Solo homogêneo:** um único conjunto de γ, γ_sat, c' e φ'. Não há camadas.
- **Ruptura apenas circular:** superfícies planares ou compostas (ex.: uma camada fraca) exigem métodos como Spencer ou Morgenstern-Price.
- **Peso pela altura no ponto médio:** cada fatia é tratada como um retângulo b × h; com 50 fatias o erro é desprezível.
- **Sem cargas externas, sismo ou reforços** (sobrecargas, tirantes, geossintéticos, estacas).
- **Orientação fixa:** o talude deve descer da esquerda para a direita (seção 4.2), e só uma face é analisada por vez.
- **Geometria simples para a detecção de crista e pé:** um único trecho descendente entre dois platôs.
- **Busca heurística:** grid + Nelder-Mead dão ótima aproximação, sem garantia matemática de mínimo global.

## 13. Considerações finais

Este manual mostrou como o notebook incorpora o nível d'água a uma análise de Bishop Simplificado: o solo abaixo do N.A. fica mais pesado, mas sobretudo a poropressão reduz a tensão efetiva e, com ela, o atrito na base das fatias. No exemplo, isso basta para levar o FS de 1,83 para 1,36 e tirar o talude da faixa considerada segura.

Ao usar o notebook com dados reais: confira a **orientação** do talude, use **parâmetros efetivos**, verifique os **gráficos de conferência** e o **mapa de FS**, e submeta os resultados à avaliação de um profissional de geotecnia antes de qualquer uso em projeto.

### Glossário de variáveis do código

| Variável | Significado | Unidade |
|---|---|---|
| `superficie`, `surf_x`, `surf_y` | pontos do perfil do terreno | m |
| `nivel_agua`, `wt_x`, `wt_y` | pontos do nível d'água | m |
| `x_crest`, `y_crest`, `x_toe`, `y_toe` | crista e pé detectados (Seção 3.1) | m |
| `slope_width`, `slope_height`, `diag` | largura, altura e diagonal da face | m |
| `margin` | folga para entrada/saída do círculo (25 % da largura) | m |
| `gamma`, `gamma_sat`, `gamma_w` | pesos específicos natural, saturado e da água | kN/m³ |
| `c`, `phi_deg`, `phi_rad` | coesão efetiva e ângulo de atrito efetivo | kPa, °, rad |
| `xc`, `yc`, `R` | centro e raio do círculo de ruptura | m |
| `b`, `h`, `h_sat`, `alpha` | largura, altura, altura saturada e ângulo da base da fatia | m, rad |
| `W`, `u`, `N_ef` | peso da fatia, poropressão na base, peso efetivo | kN/m, kPa, kN/m |
| `m_alpha` | fator mα do método de Bishop | — |
| `com_agua` | liga (True) ou desliga (False) o nível d'água | — |
| `xc_crit`, `yc_crit`, `R_crit`, `fs_crit` | círculo crítico e FS com N.A. | m, — |
| `xc_seco`, `yc_seco`, `R_seco`, `fs_seco` | círculo crítico e FS do talude seco | m, — |
