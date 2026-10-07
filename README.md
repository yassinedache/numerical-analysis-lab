# numlab : laboratoire d'analyse numérique

![tests](https://github.com/yassinedache/numerical-analysis-lab/actions/workflows/tests.yml/badge.svg)

*English version: [README_EN.md](README_EN.md)*

J'ai implémenté en Python les méthodes classiques de l'analyse numérique, puis **vérifié expérimentalement que leurs vitesses de convergence correspondent à la théorie**. Chaque résultat ci-dessous est mesuré par le code du dépôt, et des tests automatiques échouent si un ordre de convergence s'écarte de sa valeur théorique.

| Domaine | Méthodes | Ce que l'on vérifie |
|---|---|---|
| Équations non linéaires | dichotomie, point fixe, Newton, sécante | ordres 1, 1, 2 et $\varphi \approx 1{,}618$ |
| Intégration numérique | point milieu, trapèzes, Simpson, Gauss-Legendre | ordres 2, 2, 4 et convergence exponentielle |
| Équations différentielles | Euler, Heun (RK2), Runge-Kutta 4 | ordres 1, 2, 4 et conservation d'un invariant |
| Interpolation | Lagrange (barycentrique), nœuds équidistants / de Tchebychev | phénomène de Runge |
| Probabilités | Monte-Carlo | erreur en $1/\sqrt{N}$ (théorème central limite) |

## Mesurer un ordre de convergence

Si une méthode de pas $h$ a une erreur $e(h) \approx C\,h^p$, alors $\log e = \log C + p \log h$ : l'ordre $p$ est la **pente** de la droite des moindres carrés dans le plan $(\log h, \log e)$.

Pour une méthode itérative, si $e_{n+1} \approx C\,e_n^{\,q}$, on estime l'ordre par

$$q \approx \frac{\log(e_{n+1}/e_n)}{\log(e_n/e_{n-1})}.$$

Ces deux mesures sont dans [`numlab/convergence.py`](numlab/convergence.py).

## 1. Résolution de $\cos x = x$

On cherche la racine $r \approx 0{,}739085$ de $f(x) = \cos x - x$.

- **Dichotomie** : l'intervalle est divisé par 2 à chaque étape, donc $|x_n - r| \le (b-a)/2^{n+1}$.
- **Point fixe** $x_{n+1} = \cos x_n$ : convergence linéaire de taux $|g'(r)| = \sin r \approx 0{,}674$.
- **Newton** $x_{n+1} = x_n - f(x_n)/f'(x_n)$ : convergence quadratique près d'une racine simple.
- **Sécante** : Newton avec la dérivée remplacée par un taux d'accroissement, d'ordre $\frac{1+\sqrt5}{2}$.

![Erreur en fonction des itérations](figures/1_roots.png)

| Méthode | Itérations (précision $10^{-12}$) | Mesuré | Théorie |
|---|---|---|---|
| Dichotomie | 40 | ordre 1, taux 0,502 | ordre 1, taux 1/2 |
| Point fixe | 69 | ordre 1, taux 0,674 | ordre 1, taux $\sin r = 0{,}674$ |
| Sécante | 6 | 1,60 | 1,618 |
| Newton | 5 | 2,00 | 2 |

Newton atteint la précision machine en 5 itérations, contre 69 pour le point fixe : à chaque itération, le nombre de chiffres exacts double.

## 2. Intégration numérique de $\int_0^\pi \sin x\,dx = 2$

Avec $n$ sous-intervalles de pas $h = \pi/n$ :

$$T_h = h\Big(\tfrac{f(x_0)}{2} + \sum_{i=1}^{n-1} f(x_i) + \tfrac{f(x_n)}{2}\Big), \qquad S_h = \tfrac{h}{3}\Big(f(x_0) + 4\!\!\sum_{i \text{ impair}}\!\! f(x_i) + 2\!\!\sum_{i \text{ pair}}\!\! f(x_i) + f(x_n)\Big).$$

La quadrature de **Gauss-Legendre** à $n$ points choisit les nœuds (racines du polynôme de Legendre $P_n$) et les poids pour être exacte sur les polynômes de degré $\le 2n-1$.

![Erreur des méthodes d'intégration](figures/2_quadrature.png)

| Méthode | Ordre mesuré | Théorie |
|---|---|---|
| Point milieu | 2,02 | 2 |
| Trapèzes | 2,01 | 2 |
| Simpson | 4,08 | 4 |

Gauss-Legendre atteint une erreur de $1{,}1 \times 10^{-7}$ avec seulement 5 points, et la précision machine avec 10 : sa convergence est plus rapide que toute puissance de $h$, car $\sin$ est analytique.

## 3. Équations différentielles

**Ordres de convergence.** Pour $y' = -2ty$, $y(0) = 1$, de solution exacte $y(t) = e^{-t^2}$, on mesure l'erreur en $t = 2$ :

![Ordres des méthodes de Runge-Kutta](figures/3_ode_orders.png)

| Méthode | Ordre mesuré | Théorie |
|---|---|---|
| Euler explicite | 1,01 | 1 |
| Heun (RK2) | 2,12 | 2 |
| Runge-Kutta 4 | 4,11 | 4 |

**Un modèle proie-prédateur.** Le système de Lotka-Volterra

$$x' = \alpha x - \beta xy, \qquad y' = \delta xy - \gamma y$$

conserve la quantité $V(x,y) = \delta x - \gamma \ln x + \beta y - \alpha \ln y$ : ses trajectoires exactes sont des courbes fermées. Avec le même pas $h = 0{,}05$ sur $[0, 60]$ :

![Lotka-Volterra : Euler contre RK4](figures/4_lotka_volterra.png)

| Méthode | Dérive relative maximale de $V$ |
|---|---|
| Euler | $4{,}8 \times 10^{-1}$ |
| Runge-Kutta 4 | $7{,}9 \times 10^{-9}$ |

Euler fait spiraler la solution vers l'extérieur : une erreur d'ordre 1 qui s'accumule change le comportement qualitatif du modèle. RK4 garde l'orbite fermée.

## 4. Interpolation et phénomène de Runge

On interpole $f(x) = \dfrac{1}{1+25x^2}$ sur $[-1, 1]$ par un polynôme de degré $n$, avec la formule barycentrique de Lagrange (stable numériquement). L'erreur d'interpolation vaut

$$f(x) - p_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!}\prod_{i=0}^{n}(x - x_i).$$

Les **nœuds de Tchebychev** $x_k = \cos\frac{(2k+1)\pi}{2(n+1)}$ minimisent $\max|\prod (x - x_i)|$ sur $[-1,1]$.

![Phénomène de Runge](figures/5_runge.png)

| Degré | Erreur max, nœuds équidistants | Erreur max, nœuds de Tchebychev |
|---|---|---|
| 10 | $1{,}9$ | $1{,}1 \times 10^{-1}$ |
| 20 | $6{,}0 \times 10^{1}$ | $1{,}5 \times 10^{-2}$ |
| 30 | $2{,}4 \times 10^{3}$ | $2{,}1 \times 10^{-3}$ |
| 40 | $1{,}1 \times 10^{5}$ | $2{,}9 \times 10^{-4}$ |

Avec des nœuds équidistants, augmenter le degré **aggrave** l'erreur près des bords. Avec Tchebychev, elle décroît de façon géométrique.

## 5. Méthode de Monte-Carlo

On estime $\pi$ en tirant $N$ points uniformes dans $[0,1]^2$ : la proportion de points dans le quart de disque tend vers $\pi/4$ (loi des grands nombres), avec une erreur de l'ordre de $\sigma/\sqrt N$ (théorème central limite).

![Monte-Carlo](figures/6_monte_carlo.png)

Pente mesurée sur 20 tirages : **−0,50**, comme prévu. Pour $\int_0^\pi \sin x\,dx$ avec $N = 100\,000$, on obtient $2{,}0019$ avec l'intervalle de confiance à 95 % $[1{,}9959 ; 2{,}0079]$, qui contient bien la valeur exacte 2. Cette méthode est lente en dimension 1, mais son erreur ne dépend pas de la dimension, ce qui la rend précieuse en grande dimension.

## Structure du projet

```
numlab/
  roots.py          dichotomie, point fixe, Newton, sécante
  quadrature.py     point milieu, trapèzes, Simpson, Gauss-Legendre
  ode.py            Euler, Heun, Runge-Kutta 4 (équations et systèmes)
  interpolation.py  Lagrange barycentrique, nœuds équidistants et de Tchebychev
  monte_carlo.py    estimation de pi, intégrales avec intervalle de confiance
  convergence.py    mesure des ordres et des taux de convergence
experiments/
  run_all.py        lance toutes les expériences, génère figures/ et results.md
tests/              26 tests, dont la vérification de chaque ordre théorique
```

## Utilisation

```bash
pip install -r requirements.txt
python experiments/run_all.py   # régénère les figures et results.md
pytest                          # lance les tests
```

```python
from numlab import roots, quadrature
import math, numpy as np

roots.newton(lambda x: x**2 - 2, lambda x: 2*x, 1.0).root   # 1.414213562373095 (racine de 2)
quadrature.simpson(np.exp, 0, 1, 10)                          # 1.7182828 (e - 1 = 1.7182818...)
```

## Pistes pour aller plus loin

- Méthodes à pas adaptatif (Runge-Kutta-Fehlberg) et intégrateurs symplectiques, qui conservent exactement la structure de Lotka-Volterra.
- Résolution de systèmes linéaires : Gauss, LU, méthodes itératives de Jacobi et Gauss-Seidel.
- Équation de la chaleur par différences finies et étude de la condition de stabilité CFL.

## Auteur

**Yassine Abderrezak Dache**, licence d'informatique (USTHB). Projet réalisé pour approfondir l'analyse numérique et les mathématiques appliquées.
GitHub : [@yassinedache](https://github.com/yassinedache)

Licence MIT.
