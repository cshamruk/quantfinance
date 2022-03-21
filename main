---
title: Quantitative Finance
---

## AR - MA conversion

### Scalar AR(1) process

$$y_t = \rho y_{t-1} + \epsilon_t, \quad \forall t = \dots, -2, -1, 0, 1, 2, \dots $$
​
 We can substitute $y_{t-1}$ for $\rho y_{t-2} + \epsilon_{t-1}$ in the first expression to get:
​
$$y_t = \rho^2 y_{t-2} + \rho \epsilon_{t-1}  + \epsilon_t,$$

then repeat the trick with $y_{t-2},y_{t-3},\dots y_{t-k}$ to get:

\begin{align}     
  y_t =& \rho^k y_{t-k} + \rho^{k-1} \epsilon_{t-k} + \dots \rho \epsilon_{t-1}  + \epsilon_t \\     
  =&\rho^k y_{t-k} + \sum_{j=1}^{k} \rho^{j-1} \epsilon_{t-j} 
\end{align}

If we keep going, extending the number of used historical observations to infinity, $k \to \infty$, we can get the following expression for $y_t$:

$$y_t = \lim_{k\to\infty} \rho^k y_{t-k} + \sum_{j=1}^{\infty} \rho^{j-1} \epsilon_{t-j} $$

 if $|\rho| < 1$ and $y_{t-k}$ is bounded, $\lim_{k\to\infty} \rho^k y_{t-k} = 0$ and we are left with

$$y_t = \sum_{j=1}^{\infty} \rho^{j-1} \epsilon_{t-j}, $$
 
 i.e. the current value of our AR(1) process is equal to the weighted sum of the infinite number of old disturbances $\epsilon_{t-j}$, MA($\infty$).

#### Operator notation
Operator notation is very handy when reasoning about the sequences:
​
 Denote $L$ the **lag operator**, such that $$L y_t \triangleq y_{t-1}.$$
 - Integer power means consecutive application of the operator several times, e.g.
 $$L^3 y_t = y_{t-3}$$,
 - Polynomial of the lag operator means that we are taking a weighted moving average of the series, e.g.  
 $$(1 + 4 L + 5L^2 -3 L^4) y_t = y_t + 4 y_{t-1} + 5 y_{t-2} - 3 y_{t-4}.$$
​
During some intermediate linear algebra manipulations we can treat the operator as any constant, for example we can write:
​
 $$ y_t = \rho L y_t + \epsilon_t, $$
 $$ (1- \rho L) y_t = \epsilon_t,  $$
 $$ y_t = (1 - \rho L)^{-1} \epsilon_t,$$
 and using the McLaurin series expansion $\frac{1}{1-x} = 1 + x + x^2 + x^3 +\dots$, we get
​
 $$ y_t = (1 + \rho L + \rho^2 L^2 +  \rho^2 L^3 \dots ) \epsilon_t = \epsilon_t + \rho \epsilon_{t-1} + \rho^2 \epsilon_{t-2} + \dots$$

### Vector AR(1) process

 For vector processes we now deal with vector-valued observations, $z_t = \begin{bmatrix} y^1_t \\ \dots \\ y^n_t \end{bmatrix}$ and the disturbances are also a vector $u_t = \begin{bmatrix} \epsilon^1_t \\ \dots \\ \epsilon^n_t \end{bmatrix}$.

 An vector AR(1) process will be represented by a system of $n$ equations, one per component of observed vector, which will have the following matrix form :

$$z_t = \underbrace{\begin{bmatrix} \rho^{11} & \dots & \rho^{1n}  \\ \dots & \dots & \dots \\ \rho^{n1} & \dots & \rho^{n1} \end{bmatrix}}_{A} z_{t-1} + u_t.$$

 The vector AR-MA inversion takes additional steps compared to scalar inversion :
 1. **Transorm observation space:** Convert $n$-vector AR into $n$ separate scalar $AR$ processes, by applying linear transformation of the observation vector (eigendecomposition)
 2. Invert each AR(1) separately, and write the MA($\infty$) in vector form,
3. **Reverse observation space transformation:** Apply the reverse linear transformation to the $MA(\infty)$ to get back to the original observation vector.

 #### Eigenvalue decomposition and trasnformation of observation vector.

 Matrix $A$ is square, of dimension $n \times n$. It can be decomposed into the matrix product $A = Q \Lambda Q^{-1},$ where $\Lambda = \begin{bmatrix} \lambda_1 & \dots & 0 \\ 0 & \dots &0 \\ 0 &\dots & \lambda_n \end{bmatrix}$ and $Q$ has corresponding eigenvectors as columns. Using this decomposition in our system:
 
> **Refresher on matrix multiplication**
> Row vector (matrix of dimensions $1 \times n$) can be multiplied by a column vector (matrix of dimensions $n \times 1$):
> $$\begin{bmatrix}a & b \end{bmatrix} \cdot \begin{bmatrix} x \\ y\end{bmatrix} = a x + b y$$
>Matrix of dimensions $m \times n$ can be multiplied by a vector of dimensions $n \times 1$, row by row:
>$$
  \begin{bmatrix}a & b \\ c & d\end{bmatrix} \cdot \begin{bmatrix} x \\ y\end{bmatrix} = \begin{bmatrix} a x + b y \\ c x + d y\end{bmatrix}
 $$
> Diagonal matrix multiplication yeilds a vector which has entries depending only on the corresponding entries of the original vector, i.e. represents a system where each equation can be solved separately: 
>$$ \begin{bmatrix}a & 0 \\ 0 & d\end{bmatrix} \cdot \begin{bmatrix} x \\ y\end{bmatrix} = \begin{bmatrix} a x \\ d y\end{bmatrix}
$$
>
 
\begin{align} Q^{-1} \times \vert\quad && z_t &= Q \Lambda Q^{-1} z_{t-1} + u_t   \\ && \underbrace{Q^{-1} z_t}_{\tilde z _t} &= \Lambda \underbrace{Q^{-1} z_t}_{\tilde z_{t-1}} + \underbrace{Q^{-1} u_t}_{\tilde u_{t}},\end{align}


 the latter in non-matrix format is just a system of AR(1) equations we can invert separately:

 $$\begin{cases} \tilde z^1_t &= \lambda_1 \tilde z^1_{t-1} + \tilde u^1_{t}\\ &\dots\\ \tilde z^n_t &= \lambda_n \tilde z^n_{t-1} + \tilde u^n_{t} \end{cases}$$

 Note that just as we transformed the original observation vector, $\tilde z_t = Q^{-1} z_t$, we can reverse the transformation to get the original observation vector, $z_t = Q \tilde z_t$.

 #### MA($\infty$) in the transformed observation space

 After inversion we obtain the following:

 $$ \tilde z^n_t = \sum_{j=0}^{\infty} \lambda_n^j \tilde u^n_{t-j}$$, or in vector form

 $$ \tilde z_t = \sum_{j=1}^{\infty} \Lambda^j \tilde u_{t-j}.$$
 
 Power of a diagonal matrix: 
 $$
  \begin{bmatrix}a & 0 \\ 0 & d\end{bmatrix}^n  = \begin{bmatrix}a^n & 0 \\ 0 & d^n\end{bmatrix} 
 $$

 #### Reversing the transformation

 To reverse the transformation, we matrix multiply both sides of the vector MA by $Q$: 

\begin{align} Q \times \vert\quad && \tilde z_t &= \sum_{j=0}^\infty \Lambda^j \tilde u_{t-j} \\ && z_t &= \sum_{j=0}^\infty Q \Lambda^j Q^{-1} u_{t-j}, \end{align}

\begin{align}
\begin{bmatrix} y_t^1 \\ \dots \\ y_t^n \end{bmatrix} &= \sum_j \begin{bmatrix} q_{11} & \dots & q_{1n}\\ \dots &&\dots\\ q_{n1} & \dots & q_{nn} \end{bmatrix} \cdot \begin{bmatrix} \lambda^j_1 & \dots & 0 \\ \dots & & \dots \\
0 &\dots  & \lambda^j_n \end{bmatrix} \cdot \begin{bmatrix} s_{11} & \dots & s_{1n}\\ \dots &&\dots\\ s_{n1} & \dots & s_{nn} \end{bmatrix} \cdot \begin{bmatrix} u^1_{t-j}\\ \dots \\ u^n_{t-j} \end{bmatrix}  \\
&= \sum_{j=0}^\infty\begin{bmatrix} \lambda_1^j q_11 & \dots & \lambda_1^j q_{1n} \\ \dots \\\lambda_1^j q_{n1} & \dots & \lambda_n^j q_{nn}\end{bmatrix} \cdot \begin{bmatrix} s_{11} & \dots & s_{1n}\\ \dots &&\dots\\ s_{n1} & \dots & s_{nn} \end{bmatrix} \cdot \begin{bmatrix} u^1_{t-j}\\ \dots \\ u^n_{t-j} \end{bmatrix}\\
&= \sum_{j=0}^\infty \begin{bmatrix} \lambda_1^j q_{11}s_{11} + \dots + \lambda_1^j q_{1n}s_{n1} & \dots &\dots  \\ \dots \\\dots &\dots \dots \end{bmatrix} \cdot \begin{bmatrix} u^1_{t-j}\\ \dots \\ u^n_{t-j} \end{bmatrix}
\end{align}

$$y^1_t  = \sum_{j=0}^\infty \lambda_1^j \underbrace{\sum_k \cdots}_{\alpha_k} u^k_{t-j}$$




### Inverting scalar AR( p)

> Any scalar AR( p)  can be represented as a vector AR(1)


#### Example :
Take an AR(2):

$$y_t = \alpha y_ {t-1} + \beta y_{t-2} + \epsilon_t,$$

> $$ [\alpha, \beta] \cdot \begin{bmatrix}y_{t-1} \\y_{t-2}\end{bmatrix} = \alpha y_ {t-1} + \beta y_{t-2}$$

Let $z_t = [y_t, y_{t-1}]^T$, then the AR-p equation can be rewritten in vector form as (add identity $y_{t-1} = y_{t-1}$ as the second equation to the system):

$$\begin{bmatrix}y_t \\y_{t-1}\end{bmatrix} =
\begin{bmatrix}\alpha & \beta \\ 1 & 0 \end{bmatrix} \cdot \begin{bmatrix}y_{t-1} \\y_{t-2}\end{bmatrix} + \begin{bmatrix}\epsilon_t \\ 0 \end{bmatrix}
$$

Eigenvalues of matrix $\begin{bmatrix}\alpha & \beta \\ 1 & 0 \end{bmatrix}$ must be less than one in absolute value, for this vector AR(1) to be invertible. 

Eigenvalues are the roots of the equation 

$$\det \begin{bmatrix}\alpha - \lambda & \beta \\ 1 & -\lambda \end{bmatrix} = 0$$
$$ (\alpha - \lambda)\cdot (-\lambda) - 1 \cdot \beta = 0$$
$$\lambda^2 - \alpha \lambda - \beta = 0$$

Equivalent common notation is using the inverse of eigenvalues $z = \frac{1}{\lambda},$

$$1 - \alpha z - \beta z^2 = 0,$$

and the invertibility condition requires that $|z| > 1$.

> #### Shortcut:
>
>Formula for the AR-p process : $y_t = \alpha y_ {t-1} + \beta y_{t-2} + \epsilon_t,$ $\Rightarrow$ characteristic polynomial takes form 
$$1 - \alpha z - \beta z^2$$
> 
> In the general case 
> $$y_t = (\rho_1 L + \rho_2 L^2 + \dots  + \rho_p L^p) y_t + \epsilon_t,$$
>$$ (1 - \rho_1 L - \rho_2 L^2 - \dots  - \rho_p L^p) y_t = \epsilon_t$$
>**The characteristic polynomial coincides with the lag-polynomial.**
>

#### Complex characteristic roots $z$

If a root is complex, its absolute value is determined as follows:

$$ z = a+ b \cdot i, \quad |z| = \sqrt{a^2 + b^2}$$

---

### Ergodicity in the mean:

If $y_t$ is called **ergodic in the mean** if

$$\lim_{T\to \infty}\frac{1}{T} \sum_{t=0}^T y_t  = \mathbb{E}(y_t)$$

**Sufficient condition** for the ergodicity in mean:

$$\lim_{T\to \infty} \frac{1}{T}\sum_{\tau = 1}^{T}  \gamma(\tau)  = 0$$

$\gamma(\tau) = Cov(y_t, y_{t-\tau})$ the same for any $t$ for a stationary process $y_t$.

-----
### Autocorrelation and partial autocorrelation

#### MA( p) processes

$y_t = \sum_{j=0}^p \alpha_j \epsilon_{t-j},\quad \alpha_0 = 1$

\begin{align}
Cov(y_t, y_{t-2}) &= Cov(\sum_{j=0}^p \alpha_j \epsilon_{t-j},\sum_{j=0}^p \alpha_j \epsilon_{t-2-j} )\\
&= \sum \alpha_j Cov(\epsilon_{t-j}, \sum_{j=0}^p \alpha_j \epsilon_{t-2-j}) \\
&= \sum_j \alpha_j \sum_{k=0}^p \alpha_k Cov(\epsilon_{t-j},  \epsilon_{t-2-k})
\end{align}

||$\epsilon_t$|$\epsilon_{t-1}$|$\epsilon_{t-2}$|$\epsilon_{t-3}$|
|---|---|---|---|---|
|$y_t$|$1$|$\alpha_1$|$\alpha_2$||
|$y_{t-1}$|$0$|$1$| $\alpha_1$|$\alpha_2$|
|$y_{t-2}$|$0$|$0$|$1$| $\alpha_1$|$\alpha_2$|
|$y_{t-3}$|$0$|$0$|$0$|$1$| $\alpha_1$|$\alpha_2$|


$\epsilon$ are iid, withe noise disturbances, so that $Cov(\epsilon_{t-j},  \epsilon_{t-2-k}) = 0$ for all $k \neq j-2$

$$ \gamma(2) = Cov(y_t, y_{t-2}) = \sum_{j=0}^{p} \alpha_j \alpha_{j-2}$$

$$\gamma_{p+1} = \sum_{j=0}^{p} \alpha_j \alpha_{j-p-1} = \alpha_0 \alpha_{-p-1} + \dots + \alpha_p \alpha_{-1} = 0 + \dots + 0 = 0$$

**Key takeaway :** Autocovariance of order $k>p$ is zero for MA of order $p$.

> We used the bilinearity of the covariance, i.e. the property that:
>
> \begin{align}
Cov(Ax + B y, C z + D w) &= A\cdot Cov(x, C z + D w) + B\cdot Cov( y, C z + D w)=\\ 
&= AC\cdot Cov(x,z) + BC\cdot Cov(y,z) + AC\cdot Cov(x,z) 
\end{align}


#### Example : AR( p)

$$y_t = (\rho_1 L + \rho_2 L^2 + \dots  + \rho_p L^p) y_t + \epsilon_t \quad \text{[ AR(p) ]}$$

An invertible AR( p) process is efectively an MA($\infty$), i.e. there exists $\{\beta_i\}_{i=0}^\infty$ such that 

$$y_t = \sum_i \beta_i \epsilon_{t-i},$$

where $\beta_i = \sum_{k=1}^p \alpha_k \lambda_k^i$

This implies that the autocovariance of lag $\tau$ takes form:

$$\gamma(\tau) = \sum_{i=\tau}^\infty \beta_i \beta_{i-\tau}$$

As $\beta_i$ is proportional to $\lambda^i$, autocovariance $\gamma(\tau)$ is proportional to $\lambda^i \lambda^{i-\tau} = \lambda^{2 i - \tau}$, that is $\gamma(\tau)$ is decreasing exponentially. 

![](https://i.imgur.com/nvKmaBi.png)

### Yule walker equations

Autocovariance at lag $\tau$ is given by the formula:

\begin{align}
\gamma_\tau &= Cov(y_t, y_{t-\tau})\\
& = \mathbb{E}(y_t y_{t-\tau}) -\mathbb{E}(y_t)\mathbb{E}( y_{t-\tau}) \\
\end{align}

Note that the subtracted term is just zero as the expectetion of $y_t$ is zero:
$$\mathbb{E}(y_t) = (\rho_1  + \rho_2  + \dots  + \rho_p) \mathbb{E}(y_t) + 0 \Rightarrow \mathbb{E}(y_t) = 0$$

Thus, we have : 

$$\gamma_\tau = \mathbb{E}(y_t y_{t-\tau})$$

Multiplying the equation $\text([ AR(p) ])$ by $y_{t}, y_{t-1}, .. y_{t-p}$ and taking expectations we get:

\begin{align}
\tau = 0 : &\quad \mathbb{E}(y^2_t) = \rho_1 \mathbb{E}( y_{t-1} y_{t}) + \rho_2 \mathbb{E}(y_{t-2} y_{t}) + \dots + \rho_p\mathbb{E}(y_{t-p} y_{t}) + \underbrace{\mathbb{E}(\epsilon_t y_{t})}_{= \mathbb{E}\epsilon^2_t = \sigma^2} \\
& \gamma_0 = \rho_1 \gamma_1 + \rho_2 \gamma_2 + \dots \rho_p \gamma_{p}\\ \\
\tau = 1 : &\quad \mathbb{E}(y_t y_{t-1}) = \rho_1 \mathbb{E}( y^2_{t-1}) + \rho_2 \mathbb{E}(y_{t-2} y_{t-1}) + \dots + \rho_p\mathbb{E}(y_{t-p} y_{t-1}) + \underbrace{\mathbb{E}(\epsilon_t y_{t-1})}_{=0}\\
& \gamma_1 = \rho_1 \gamma_0 + \rho_2 \gamma_1 + \dots \rho_p \gamma_{p-1}\\ \\
\tau = 2 : &\quad \gamma_2 = \mathbb{E}(y_t y_{t-2})= \rho_1 \mathbb{E}(y_{t-1} y_{t-2}) + \rho_2 \mathbb{E}( y^2_{t-2}) + \dots + \rho_p \mathbb{E}(y_{t-p} y_{t-2})\\\\
& \gamma_2 = \rho_1\gamma_1 + \rho_2 \gamma_0 + \dots + \rho_p \gamma_{p-2}\\
&\dots\\
&\gamma_p = \rho_1 \gamma_{p-1} + \rho_2 \gamma_{p-2} + \dots + \rho_p \gamma_0
\end{align}

$$
\begin{cases}
\gamma_1 = \rho_1 \gamma_0 + \rho_2 \gamma_1 + \dots \rho_p \gamma_{p-1}\\
\gamma_2 = \rho_1\gamma_1 + \rho_2 \gamma_0 + \dots + \rho_p \gamma_{p-2}\\
\dots\\
\gamma_p = \rho_1 \gamma_{p-1} + \rho_2 \gamma_{p-2} + \dots + \rho_p \gamma_0
\end{cases}
$$

The system of Yule-walker has $p$ equations  with $p$ unknowns, given sample covarinces $\hat \gamma$ we can solve it with respect to $\rho$ to figure out the coefficients of our AR process.

System of $k$ Yule walker equations yeilds a solution $[\phi_{k,1}, \dots, \phi_{k,k}]$.
By increasing the size of the system and keeping each time only $\phi_{k,k}$ component of the solution, we get the sequence PACF $\{\phi_{k,k}\}_{k=1}^\infty$.

![](https://i.imgur.com/FUK7LyS.png)


### Dickey-Fuller test for Unit root

take AR(1) process:

$$y_t = y_{t-1} + \epsilon_t$$

$$(1-L) y_t  = \epsilon_t$$

The $\lambda$ of this polynomial is exactly 1.

This means that this process has an exploding variance and autocovariance.

If we invert it, we get MA($\infty$) with coefficients $\beta_i$ proportional to $1$ - the effect of very old errors does not vanish as for the invertible AR process.

Autocovariance is an infinite sum of weighted betas, and as betas are constant, do not decrease exponentially, the infinite sum tends to +/- infinity. This means that we get no ergodicity, and no empirical estimate is correct or reliable. 

If there is a unit root, we cannot trust our estimation procedures to take care of it, we must take the first difference BY HAND.

$$y_t = \rho y_{t-1} + \epsilon_t $$

$$y_t - y_{t-1} = (\rho - 1) y_{t-1} + \epsilon_t $$

**H0** : 
- $\epsilon_t \sim iid N(0,\sigma^2)$ 
- AND $y_t - y_{t-1} =  (\rho - 1) y_{t-1} + \epsilon_t$ 
- AND $\rho-1 = 0$.

**Ha** : 
- any of these bullet points fail.
 
Under H0 take t-statistic of $\hat{\rho}-1$:

$$\frac{\hat \rho - 1}{\hat{Var(\hat \rho)}}$$

![](https://i.imgur.com/P9UAnsH.png)

If H0 is not rejected - we do have a unit root.

**Rejection of H0 does not imply absence of unit root**, because H0 can be rejected for two other reasons - errors are not normal, or the model equation is incorrect (for example, MA term is specified incorrectly, or intercept is really zero, or we forgot to include a determinstic trend.)

For this reason, Dickey Fuller test is usually performed in several different variations before moving on to the conclusion. Variations are usually different combinations of the following specifications:

- DF wit hintercept : $y_t - y_{t-1} = \rho_0 + (\rho - 1) y_{t-1} + \epsilon_t$
- DF with constant and a trend : $y_t - y_{t-1} = \rho_0 + \delta t + (\rho - 1) y_{t-1} + \epsilon_t$
- Augmented Dickey-Fuller : adding AR (and MA) terms.  


#### Augmented Dickey-Fuller

In general a model with unit root, ARIMA($\cdot, 1, \cdot$), has form

$$(1 - \alpha_1 L - \alpha_2 L^2 +...) (1-L) y_t = (1 + \beta_1 L + \beta_2 L^2 +\dots )\epsilon_t$$

Augmented Dickey Fuller takes care of such specifications:

$$(1 - \alpha_1 L - \alpha_2 L^2 +...) \Delta y_t = (1 + \beta_1 L + \beta_2 L^2 +\dots )\epsilon_t,$$ or equivalently 
$$\Delta y_t = \alpha_1 \Delta y_{t-1} + \alpha_2 \Delta y_{t-2} +... +\epsilon_t + \beta_1 \epsilon_{t-1} + \beta_2 \epsilon_{t-2} +\dots$$

**Augmented DF** test adds $y_{t-1}$ to the right side of the equation:

$$\Delta y_t = (\rho-1)y_{t-1} + \underbrace{\alpha_1 \Delta y_{t-1} + \alpha_2 \Delta y_{t-2}}_\text{AR terms} +... +\epsilon_t + \underbrace{\beta_1 \epsilon_{t-1} + \beta_2 \epsilon_{t-2} +\dots}_\text{MA terms}.$$


## Forecasting

### 1. ARIMA fitting outcomes

- outcomes : coefficient estimates $\{\hat \alpha_k\}_{k=0}^{p}, \{\hat \beta_k\}_{k=1}^{q}$ for ARMA(p,q), $\{\hat \epsilon_t\}_{t=0}^T$.
- in-sample prediction : $\{\hat \epsilon_t\}_{t=0}^T$ is estimated using data $\{y_t\}_{t=0}^T$

### 3. Loss functions
- optimality = minimization of expected loss
Having a model of data, some idea of probability distribution of the actual data $y_t$, we seek to minimize the **expected** loss of prediction. 
$$L(y_{t+h},\hat y_{t+h|t})$$
$L()$ is higher if $|y_{t+h} - \hat y_{t+h|t}|$ grows.
- types of loss functions
    - symmetric
        MSE $(y_{t+h}-\hat y_{t+h|t})^2$, MAE $|y_{t+h} - \hat y_{t+h|t}|$
    - asymmetric
        
    - if $y_t$ can tale only two values "good", "bad", then you can use *entropy loss*.
    
### 2. Forecasts
- MSE-loss forecast
MSE-loss based optimal forecast is a lways just the expectation of $y_{t+h}$ based on information available at $t$ an our estimated model.

How the forecasts are done?:
**Input:** $\{y_t\}_{t=0}^T$
1. Estimate ARIMA model.
    - Check for unit roots using ADF
        - with trend, with intercept, and for a range of different ARMA structures, i.e. determine $d$.
    - Use ACF and PACF to narrow down the options for $p,q$, orders of AR and MA.
    - Use any of the estimation methods to get the coiefficient and error estimates of ARIMA(p,d,q)
2. $\hat y_{t+h|t} = \mathbb{E}\left(y_{t+h}\Bigg|\{y_t\}_{t=0}^T, \{\hat \alpha_k\}_{k=0}^{p}, \{\hat \beta_k\}_{k=1}^{q},\{\hat \epsilon_t\}_{t=0}^T\right)$

## Examples

### 1.a

$y_t = 0.6 y_{t-1} - 0.3 y_{t-2} + \epsilon_t$

$y_t - 0.6 y_{t-1} + 0.3 y_{t-2} = \epsilon_t$

$(1 - 0.6 L +0.3 L^2 )y_t = \epsilon_t$

Charactestic polynomial is $1 - 0.6 L +0.3 L^2$, its roots are complex $z_{1,2} = 1 \pm i\sqrt{\frac{0.84}{0.36}}$. 
For invertibility $z_{1,2}$ must "lie outside the unit disk", i.e. $||z|| > 1$
-$||z_1|| = \sqrt{1^2 + \frac{0.84}{0.36}} =\sqrt{1.2/.36} > 1$
-$||z_2|| = \sqrt{1^2 + \frac{0.84}{0.36}} =\sqrt{1.2/.36} > 1$
Both roots are outside of unit disk and the AR is invertible. Therefore stationary, coavariance stationary, ergodic etc.

### 1.b

$1 - .3 L - .3 L^2 = 0$

Roots : $z_{1,2} = \frac{1 \pm \sqrt 1.29}{-.6}$
- $z_1  = - 2.14/.6$ - greater than 1 in absolute value
- $z_2 = -(1 - 1.14)/.6 = -.14/.6$ - is less than 1 in absolute value, therefore it is inside the unit dist, and AR is not invertible.


### 1.c

$1 + 1.2 L + L^2 = 0$

The roots: $z_{1,2} = \frac{-1.2 \pm i\sqrt{2.56}}{2}$
- $|z_1| = 2.56/4 + .6^2 = .64 + .36 = 1$ - unit root
- If your polynomial has a complex root $a + ib$, there always exists the second complex root $a - ib$, and their absolute values are the same.
$x^2 + 1 = 0$ has roots $x_{1,2} = \pm i$.

$d=2$ for this process.

### 1.d
finite order MA is always stationary.

### 2.

$y_t = .5 y_{t-1} + \epsilon_t$

$\hat y_{t+1} = \mathbb{E}_ty_{t+1} =\mathbb{E}_t(.5 y_t + \epsilon_{t+1})$

at $t=100$ we know:
- $y_t =-.2$
- $\epsilon_{101} \sim WN(0,1)$, i.e. $\mathbb{E}_t \epsilon_{t+1} = 0$

$\hat y_{t+1} = \mathbb{E}_t(.5 y_t + \epsilon_{t+1}) = .5 \cdot y_t + 0$
$\hat y_{t+2} = \mathbb{E}_t(.5 y_{t+1} + \epsilon_{t+2}) = \mathbb{E}_t(.5 y_{t+1}) = \mathbb{E}_t(.5 (.5 y_t + \epsilon_{t+1})) = .25 y_t + 0$
$\hat y_{t+3} = \mathbb{E}_t(.5 y_{t+2} + \epsilon_{t+3}) =.5\mathbb{E}_t y_{t+2} + \mathbb{E}_t\epsilon_{t+3} = .5 \hat y_{t+2} =0.125 y_t$

$y_{103|100} =0.125*(-.2) =-.025$

#### b. j-step ahead forecast

$\hat y_{t+j} = \mathbb{E}_t{y_{t+j}} = \mathbb{E}_t(.5y_{t+j-1}+ \epsilon_{t+j}) = .5 \hat y_{t+j-1}$

By iterated substitutions (or by math induction) you can see that : $\hat y_{t+j} = .5^{j} y_t$

To get unconditional expectation $\mathbb{E}y_t$, I take expectations of both sides of our model:

$$\mathbb{E}y_t = .5\mathbb{E}y_{t-1} + \mathbb{E}\epsilon_t$$
Unconditional expectation of $y_t$ does not depend on $t$ for a stationary process:
$$\mathbb{E}y_t = .5\mathbb{E}y_{t}$$
$$\mathbb{E}y_t = 0 \quad vs \quad \lim_{j \to \infty} \hat y_{t+j}=\lim_{j \to \infty} .5^j y_t = 0$$

#### c.

$$e_{t+j} = y_{t+j} - \hat y_{t+j} = y_{t+j} - .5^j y_t$$

$$e_{t+j} = .5 y_{t+j-1} + \epsilon_{t+j} - .5^j y_t = .5^{j-1} y_t + \sum_{k=0}^{j-1}.5^k\epsilon_{t+j-k} - .5^j y_t = $$

$\mathbb{E}_t e_{t+j} = 0$

$$(1-.5L)y_t = \epsilon_t\quad \Bigg| \div (1-.5L)$$
$$y_t = \frac{1}{1 - .5 L} \epsilon_t  = (1 + .5 L + ^2.5 L^2 + \dots )\epsilon_t$$

\begin{align}
y_{t+j} &= (1 + .5 L + ^2.5 L^2 + \dots )\epsilon_{t+j} \\
&=(1 + .5L + \dots + .5^{j-1} L^{j-1})\epsilon_{t+j} + .5^jL^j(1 + .5L + \dots )\epsilon_{t+j}\\
& = (1 + .5L + \dots + .5^{j-1} L^{j-1})\epsilon_{t+j} + .5^j\underbrace{(1 + .5L + \dots )\epsilon_t}_{=y_t}\\
&= (1 + .5L + \dots + .5^{j-1} L^{j-1})\epsilon_{t+j} + .5^j y_t
\end{align}

$$e_{t+j} = y_{t+j} - \hat y_{t+j} = \sum_{k=0}^{j-1}.5^k\epsilon_{t+j-k}$$

Variance of forecast error is obtained by direct application of the variance formula to the last expression

> $$Var(aX + bY) = a^2 Var(X) + b^2 Var(Y)$$
> Partial sum of a gemetric progression: $$\sum_{t=0}^{T} \rho^t= \frac{1 -\rho^{T+1}}{1-\rho}\to_{T\to \infty} 1$$

$$Var_t(e_{t+j}) = \sum_{k=0}^{j-1} .5^{2k} \underbrace{Var_t(\epsilon_{t+j-k})}_{=1} = \sum_{k=0}^{j-1} .25^{k} = \frac{1 - .25^j}{1-.25} = \frac{4}{3}(1 - .25^j) \to_{j\to \infty} \frac{4}{3}$$

### 3. Forecasting pure MA processes

$$y_t = 1 + \epsilon_t - .4 \epsilon_{t-1}$$

$$y_{t+h|t} = \mathbb{E}_t(y_{t+h}) = \mathbb{E}(y_{t+h} |y_t, y_{t-1}, \dots ,y_{t-n}, \dots)$$

Ma process is invertible under some good conditions:

Just as AR(p ) can be represented as MA($\infty$),, MA(q ) can be represented as AR($\infty$), if coeffients allow.

$$\text{MA : } y_t = \alpha_0 + \epsilon_t + \alpha_1 \epsilon_{t-1} + \dots + \alpha_q \epsilon_{t-q}$$
$$y_t = \alpha_0 + (1 + \alpha_1 L + \alpha_2 L^2 +\dots +\alpha_q L^q)\epsilon_t$$
$$\epsilon_t = \frac{1}{1 + \alpha_1 L + \alpha_2 L^2 +\dots +\alpha_q L^q}(y_t - \alpha_0)$$

> Main theorem of algebra says that any polynomial can be represented as a product of monomials:
> 
> $$1 + \alpha_1 L + \alpha_2 L^2 +\dots +\alpha_q L^q = (1 - z_1L)\cdot(1 - z_2 L)\cdot\dots\cdot(1 - z_q L)$$

$$\epsilon_t = \frac{1}{(1 - z_1L)\cdot(1 - z_2 L)\cdot\dots\cdot(1 - z_q L)}(y_t - \alpha_0)$$

We can invert the right hand side in q steps, monomial by monomial:

> $$y_t = \beta_0 + \rho y_{t-1} + \epsilon_t \Rightarrow y_t = \frac{\beta_0}{1 - \rho} + \frac{1}{1 - \rho L}\epsilon_t = \frac{\beta_0}{1 - \rho} + \sum_{j=0}^{\infty} \rho^j \epsilon_{t-j} $$

$$\epsilon_t = \frac{1}{(1 - z_1L)\cdot(1 - z_2 L)\cdot\dots\cdot}\left(\frac{1}{(1 - z_q L)}(y_t - \alpha_0)\right) = \frac{1}{(1 - z_1L)\cdot(1 - z_2 L)\cdot\dots\cdot} \sum_{j=0}^\infty z_q^j (y_{t-j}-\alpha_0) $$

$$\epsilon_t = \frac{\alpha_0}{1 - \alpha_1 - \dots \alpha_q} + \sum_{j=0}^{\infty} \phi_j y_{t-j}$$

If the MA lag polynomial has unit roots inside the unit circle, i.e. if it is invertible, $\epsilon_t$ can be approximated via a weighted sum of $y_t,\dots,y_{t-K}$, K historical observations of our variable of interest, and the approximation error (due to truncation of history) will be limited (small enough).

$$y_{t+h|t} = \mathbb{E}_t(y_{t+h}) = \mathbb{E}(y_{t+h} |y_t, y_{t-1}, \dots ,y_{t-n}, \dots)$$

I.e. observation of some history of $y_t$ is equivalent to observing the current value of perturbation $\epsilon_t,\epsilon_{t-1},\dots $.

When calculating a forecast $t+h|t$, you should treat $\epsilon_t, \epsilon_{t-1},\dots$ as known value, and $\epsilon_{t+1},\epsilon_{t+2},\dots$:

\begin{align}
y_{t+1|t} &= E(1 + \epsilon_{t+1} - .4 \epsilon_{t}| y_t, y_{t-1}\dots) \\
&= E(1 + \epsilon_{t+1} - .4 \epsilon_{t}| \epsilon_t, \epsilon_{t-1}\dots)\\& = 1 + \underbrace{E(\epsilon_{t+1}| \epsilon_t, \epsilon_{t-1}\dots)}_{=0} - .4 \epsilon_t \\
&= 1 - .4 \epsilon_{t} \\
&\\
y_{t+2|t} & = \mathbb{E}_t(1 + \epsilon_{t+2} - .4 \epsilon_{t+1}) = 1 \\
y_{t+3|t} & = \mathbb{E}_t(1 + \epsilon_{t+3} - .4 \epsilon_{t+2}) = 1
\end{align}

**Conclusion :** Best forecast of MA(q) on a horizon $\geq q$ is the constant.

$$y_{101|100}= 1 -.4*.2 = .92 \quad y_{102|100} = y_{103|100} = 1$$

### 4.b

$$y_{t+1|t} = \mathbb{E}_t(-1.3 - .6 y_{t} + \epsilon_{t+1}  + .4 \epsilon_{t}) = -1.3 - .6 y_{t} + .4 \epsilon_{t}$$

Forecast error:

$$e_{t+1|t} = y_{t+1} - y_{t+1|t} = $$

\begin{align}
e_{t+1|t} & = y_{t+1} - y_{t+1|t} \\
& = -1.3 - .6 y_{t} + \epsilon_{t+1}  + .4 \epsilon_{t} - (-1.3 - .6 y_{t} + .4 \epsilon_{t}) \\
&= \epsilon_{t+1}\\
\mathbb{V}_t(e_{t+1|t}) =& \mathbb{V}_t(\epsilon_{t+1}) = \mathbb{V}(\epsilon_{t}) = \sigma^2 
\end{align}

### 6.

Conjecture : $p,d,q = 4, \dots, 1$

Order of actions:
1. Determine $d$
    1. If $\beta(1) = 0$, then $d\geq 1$, find $\beta^{(1)}(L) = \frac{\beta(L)}{1-L}$;
        otherwise $d=0$
    2. If $\beta^{(1)}(1)=0$, then $d\geq 2$ find $\beta^{(2)}(L) = \frac{\beta^{(1)}(L)}{1-L}$;
        otherwise $d=1$
    3. Continue until $\beta^{(k)}(1) \neq 0$, when we are sure that $d=k$.
2. Take $d$th order difference of the series:
$$\beta^{(d)}(L)\Delta^d y_t = \alpha(L)\epsilon_t$$
3. $p$ is the order of the last polynomial $\beta^{(d)}(L)$

> Example: 
> $$y_t = y_{t-1} + \epsilon_t$$
> This is an ARIMA(0,1,0) process.


$$G_t = 4 + 3.5 G_{t-1} -4.5 G_{t-2} +2.5 G_{t-3} - .5 G_{t-4} + \epsilon_t -0.7 \epsilon_{t-1}$$

Characteristic polynomial:

$$\beta(L) = 1 - 3.5 L + 4.5 L^2 - 2.5 L^3 + .5 L^4 $$

To check, if there is a unit root, calculate $\beta(1)$:


$$\beta(1) = 1 - 3.5 + 4.5 -2.5 +.5 = 0,$$

$L=1$ is one of the roots (the "unit root") of the polynomial, $d \geq 1$, and therefore we have to take the difference of our series:

$$\beta(L) = (1 - L)(1-2.5L + 2L^2 -.5L^3)$$

|$1 - 3.5 L + 4.5 L^2 - 2.5 L^3 + .5 L^4$|$1 - L$|
|---|---|
|$1 - 3.5 L + 4.5 L^2 - 2L^3$|$-.5L^3$|
|$1 - 3.5L + 2.5 L^2$|$2L^2 -.5L^3$|
|$1 - L$|$-2.5L + 2L^2 -.5L^3$|
|$0$|$1-2.5L + 2L^2 -.5L^3$|

There is still a unit root, as $1-2.5\cdot 1 + 2\cdot 1^2 -.5\cdot 1^3 = 0$

|||
|---|---|
|$1 - 2.5 L + 1.5 L^2$|$.5L^2$|
||$-1.5L+ .5L^2$|

$$\beta(L) = (1-L)^3(1-.5L)$$

After decomposing the lag polynomial into monomials we see, that the series $G$ is an instance of ARIMA(1,3,1):

$$\Delta^3 G_t=\Delta^2 (G_{t} - G_{t-1}) = \dots$$

$$(1 - .5L) \Delta^3 G_t = \epsilon_t -0.7 \epsilon_{t-1}$$
$$\Delta^3 G_t = .5 \Delta^3 G_{t-1} + \epsilon_t -0.7 \epsilon_{t-1}$$

$\Delta^3 G$ is ARMA(1,1).



## Volatility models

![](https://i.imgur.com/t45kmmE.png)

$$y_t = \mu_t + \epsilon_t,$$

where $\mu_t = \beta(L)y_t$, and $\epsilon_t = \sigma_t \nu_t$, with $\nu_t \sim N(0,1)$, and 

$$\sigma^2_t = \alpha_0 + \alpha_1 \epsilon_{t-1}^2 + \dots + \alpha_k \epsilon_{t-k}^2.$$

Conditional on history by date $t$, $\epsilon \sim N(0,\sigma_t^2),$ $Var(\epsilon_t) = \mathbb{E}(\epsilon^2_t) = \mathbb{E}(\mathbb{E_t}(\epsilon^2_t))=\mathbb{E} (\mathbb{E}_t\sigma_t^2) = \mathbb{E} (\sigma_t^2)$.

Unconditional variance of $\sigma^2_t$:

$$\mathbb{E}\sigma^2_t = \alpha_0 + \alpha_1 \mathbb{E}\epsilon_{t-1}^2 + \dots + \alpha_k \mathbb{E}\epsilon_{t-k}^2 $$

$$\mathbb{E}\sigma^2_t = \alpha_0 + \alpha_1 \mathbb{E}\sigma^2_t + \dots + \alpha_k \mathbb{E}\sigma^2_t $$

$$
\mathbb{E}\sigma^2_t = \frac{\alpha_0}{1 - \alpha_1 - \dots - \alpha_k}
$$

Requirements for the ARCH model to make sense:

1. $\alpha_0 > 0$, so that a series of $\epsilon_{t-k}=\dots = \epsilon_{t-1}=0$ does not turn $\epsilon_t$ into a non-random variable
2. $\alpha_1, \dots, \alpha_k \geq 0$ so that conditional variance stays positive
3. $\sum_{j_1}^{k} \alpha_j < 1$ - so that unconditional variance stays positive

If estimation of the model results in coefficients $\alpha$ that violate any of these conditions - you have find another specification, this one cannot adequately describe stochastic data.

ARCH estimation:

0. Determine ARIMA structure of the main series
    - ADF to check for trend and unit roots
    - PACF - to check for $p$, the order of autoregression
    - ACF on the residuals to check for $q$ MA order

1. Based on calculated regression residuals $\hat \epsilon_t = e_t$, construct series $\{\hat \epsilon_t^2\}$ 
    - PACF to determine the AR order $k$   
    - and run regression 
$$\hat \epsilon_t^2 = \alpha_0 + \alpha_1\hat \epsilon^2_t + \dots + \alpha_k\hat \epsilon_{t-k}^2$$

### GARCH model

GARCH model allows to concisely represent higher order ($\infty$) ARCH models,

$$\sigma^2_t = \alpha_0 + \alpha_1 \epsilon_{t-1}^2 + \gamma_1 \sigma_{t-1}^2$$

$$\text{GARCH(1,1) :} \quad (1 - \gamma_1 L) \sigma_t^2 = \alpha_0 + (\alpha_1 L ) \epsilon_t^2$$
$$\sigma_t^2 = \frac{\alpha_0}{(1 - \gamma_1 L)} + \frac{\alpha_1 L }{(1 - \gamma_1 L)}\epsilon_t^2 = \frac{\alpha_0}{1- \gamma_1} + \alpha_1 \sum_{j=0}^{\infty}\gamma_1^j \epsilon_{t-1-j}^2 $$


Invertibility requirements remain active: $0 <\gamma_1 < 1$ for the infinite sum to converge, and for both variances, unconditional and conditional, to be positive. 

### Leverage effect

In vanilla (G)ARCH volatility $\sigma^2_t$ only depends on the magnitude of $\epsilon_t$ (absolute value of the price move).

In reality it is often observed that volatility increases significantly after big downward moves in the prices and insignificantly after equivalent upward moves.

Simplest way to introduce leverage effect: introduce a dummy variable - Threshold ARCH

![](https://i.imgur.com/xMes4yk.png)

Another way: let the $\log(\sigma^2_t)$ folow the ARCH structure, then $\sigma^2_t$ will have the desired asymmetry:

$$\ln\sigma_t^2 = \alpha_0 + \sum_{i=1}^p {\alpha_i \left(\left|\epsilon_{t-i}\right|+\gamma_i\epsilon_{t-i}\right )}+\sum_{j=1}^q{\beta_j \ln\sigma_{t-j}^2}$$

- two sources of asymmetry:
    1. logarithm
    2. adding absolute value of $\epsilon_i$
    
TGARCH does nest the GARCH model - if $\alpha_i^- = 0$
EGARCH does not due to nonlinearity introduced by the logarithm.

### Restricted vs unrestricted models

TGARCH :

$$\sigma_t^2 = \alpha_0 + \alpha_1 \epsilon^2_{t-1}+\color{red}{\lambda \delta _t \epsilon^2_{t-1}}+\beta_1 \sigma_{t-1}^2,$$

TGARCH adds  a term $\lambda \delta_t \epsilon^2_{t-1}$ into the expression for GARCH(1,1), meaning that TGARCH is an unrestricted (richer) model and GARCH is the restricted one.

Estimation is done by maximizing the (log)likelihood.

Unrestricted model will have a higher likelihood.

We pick the simplest model (the restricted one) if the loss of likelihood is statistically insignificant.

Loss of likelihood is measured in relative terms by the fraction $\frac{L_u}{L_r}$ (conceptually). In practice we go further, we take the logarithm of this ratio:

$$\text{LR-statistic} = -2 \ln(\frac{L_r}{L_u}) \sim \chi^2(m),$$

where $m$ is the number of restrictions (1 if we test for TGARCH against GARCH(1,1)).

Where does 2 come from?

$$l(\theta_R) = l(\theta_U) + l^\prime(\theta_U)(\theta_R-\theta_U) + \frac{1}{2} l^{\prime\prime}(\theta_R-\theta_U)^2 + \dots$$

Let $l(\theta)$ be the log-likelihood function (dependence on the sample $ \{y_1, \dots, y_T\}$ is implicit).

Maximum of log-likelihood for the unrectricted model is achieved for $\theta_U$ such that:

$$l^\prime(\theta_U) = 0,$$

therefore our second order expansion is 

$$l(\theta_R) \approx l(\theta_U) + \frac{1}{2} l^{\prime\prime}(\theta_U)(\theta_R-\theta_U)^2 $$

$$ - l^{\prime\prime}(\theta_U)(\theta_R-\theta_U)^2 = -2(l(\theta_R) - l(\theta_U)) = -2 \ln (\frac{L_R}{L_U})$$

### Example:

$$z_t = u_t + \alpha_1 u_{t-1} + \alpha_2 u_{t-2}$$

Autocorrelation:

$$\rho(z_t,z_{t-j}) = \frac{Cov(z_t,z_{t-j})}{\sqrt{Var(z_t)Var(z_{t-j})}} = \frac{Cov(z_t,z_{t-j})}{Var(z_t)}$$

$$
Cov(z_t,z_{t-j}) = \begin{cases}
\alpha_1 \sigma^2 + \alpha_1 \alpha_2 \sigma^2, \text{ if $j=1$},
\end{cases}
$$


### Exam problem

We have two processes 
$$X_t = u_t + \delta u_{t-1} \quad Y_t = X_t + \nu_t$$

We can represent $Y_t$ as an MA process:

$$Y_t = \epsilon_t + \theta \epsilon_{t-1}$$

$\theta$ and $\epsilon$ are such taht 

$$u_t + \delta u_{t-1} + \nu_t = \epsilon_t + \theta \epsilon_{t-1} \tag{identity}$$

It is obvious that there gonna be coefficients $a,b$ such that $\epsilon_t = a u_t + b \nu_t + [c u_{t-1} + d \nu_{t-1}]$, we have only have to find these coefficients:

$$u_t + \delta u_{t-1} + \nu_t + 0 \cdot \nu_{t-1} = a u_t + b \nu_t + \theta a u_{t-1} +  \theta b \nu_{t-1}$$

The above equation must hold for all values of $u_t, u_{t-1}, \nu_t$

$$
\begin{cases}
a &= 1\\
\delta &= \theta a\\
1 &= b\\
0 &= \theta b 
\end{cases}
$$

The remaining question is whether there exists a $\theta$ such that this system has a solution w.r.t. $a,b$? The first and the third equations impose that $a=b=1$ and the remaining two are inconsistent, as they imply $\theta = \delta = 0$.


Take $\epsilon_t = \alpha(L)u_t + \beta(L)\nu_t$. Substituting this into the equality (identity) we get the following equation:

$$
u_t + \delta u_{t-1} + \nu_t + 0 \cdot \nu_{t-1} = \alpha(L) u_t + \beta(L) \nu_t + \theta L \alpha(L) u_t +  \theta L \beta(L)\nu_{t}\\
(1 + \delta L) u_t + (1 + 0\cdot L)\nu_t=(1+\theta L)\alpha(L)u_t + \beta(L)(1 + \theta L )\nu_t\\
$$

The coefficients of the lag polynomials  $\alpha,\beta$ must satify the system

$$
\begin{cases}
\alpha(L) &= \frac{1 + \delta L}{1+\theta L}\\
\beta(L) &= \frac{1}{1+\theta L}
\end{cases}
$$

Knowing that $\alpha(L) = \alpha_0 + \alpha_1 L + \alpha_2 L^2 + \dots$

$$
\begin{cases}
\alpha_0 + \alpha_1 L + \alpha_2 L^2 + \dots &= (1 + \delta L)\sum_{i=0}^\infty (-\theta L)^i \\
\beta_0 + \beta_1 L + \beta_2 L^2 + \dots &= \frac{1}{1+\theta L} = 1 - \theta L + \theta^2 L^2 -\theta^3 L^3 +\dots
\end{cases}
$$

We use Taylor expansions:

> **Reminder :** Taylor expansion around $x=0$ $$\frac{1}{1 + bx} = f(0) + f^\prime(0) (x-0) + \frac{1}{2!}f^{\prime\prime}(0)(x-0)^2 + \dots = \\= \frac{1}{1 + b\cdot 0} + \frac{-b}{(1 + b\cdot 0 )^2} x + \frac{1}{2!} \frac{2 b^2}{(1 + b\cdot 0)^3} x^2 + \dots =\\= 1 - bx + b^2 x^2 - b^3 x^3 + \dots =\\= \sum_{i=0}^\infty (-bx)^i$$
> Sum of geometric progression is a special case of Taylor expansion formula $$\frac{1}{1 - \rho L} = 1 + \rho L + \rho^2 L^2 + \dots = \sum_{i=0}^\infty (\rho L)^i $$


$$
\sum_{i=0}^\infty (-\theta L)^i + \delta L \sum_{i=0}^\infty (-\theta L)^i = \sum_{i=0}^\infty (-\theta L)^i + \sum_{i=0}^\infty (-\theta)^i \delta L^{i+1} = \\
= 1 - \theta L + \theta^2 L + \dots + \\
\quad + \delta L - \delta \theta L^2 + \delta \theta^2 L^3 \\
= 1 + (\delta - \theta)L + (\theta^2 - \delta \theta)L^2 + (\delta \theta^2 - \theta^3)L^3 + \dots$$

$$
\begin{cases}
\alpha_0 &= 1 \\
\alpha_1 &= \delta - \theta\\
\alpha_2 &= \delta \theta - \theta^2 \\
\alpha_i &= (\delta-\theta)(-\theta)^{i-1}
\end{cases}
$$

Quastion: $X,Y$ are ARMA processes with perturbations sequences $u_t, \nu_t$,  And $Z = aX + bY$. Show that $Z$ is ARMA as well and find its representation (coefficients and perturbation $\epsilon_t$ as a function of ) 
1. Say that there will be some lag polynomials $\alpha(L), \beta(L)$
2. The perturbations of the $Z$ process must be a combination of perturbations $u,\nu$:
$$\epsilon_t = \alpha(L)u_t + \beta(L) \nu_t$$
3. Substitute the above expression into the formula for $Z$, this gives you an identity that must hold for all possible $u_t,\nu_t$
4. Formulate a system of equations defining the coefficients of the lag polynomials.

We can pin down the value of $\theta$ by matching the autocovariances of $Y_t = \epsilon_t + \theta \epsilon_{t-1}$ and $Y_t = X_t + \nu_t$.

$$\text{$\gamma_s$ : }\quad Cov(\epsilon_t + \theta \epsilon_{t-1}, \epsilon_{t-s} + \theta \epsilon_{t-1-s}) = Cov(X_t + \nu_t, X_{t-s}+ \nu_{t-s}) $$

$$\text{$\gamma_0 : $}\quad \sigma_\epsilon^2 + \theta^2 \sigma_\epsilon^2 = \sigma^2_\nu + \delta^2 \sigma^2_u + \sigma^2_u$$
$$\text{$\gamma_1$ : } \quad \theta \sigma_\epsilon^2 = \delta \sigma_u^2$$

We have two equations and two unknown coefficients $\theta, \sigma^2_u$.

$\theta$ solves the equation:

$$(\frac{1}{\theta} + \theta) = \frac{1}{\delta}\frac{\sigma_\nu^2}{\sigma_u^2} + \frac{1}{\delta} + \delta.$$

Note that if $\frac{\sigma_\nu}{\sigma_u} \to 0$ then $\theta \to \delta$. Otherwise measurment error $\nu$ makes us observe a biased MA process (with MA(1) coefficient slightly smaller in absolute value than the original $\delta$.)

---

## Homework 3

### Problem 1

#### a)

\begin{align}
Var_t(Y_{t+1}) &= Var_t(b_0 + b_1 Y_t + u_{t+1}) =\\
&=Var_t(u_{t+1}) \quad \text{ because $Y_t$ is known}\\
&= Var_t(\sigma_{t+1} \nu_{t+1}) \\
& = Var_t(\sqrt{\alpha_0 + \alpha_1 u_t^2 + \alpha_2 u_{t-1}^2} \nu_{t+1})\\
& = (\alpha_0 + \alpha_1 u_t^2 + \alpha_2 u_{t-1}^2) Var_t (\nu_{t+1}) \\
&= \alpha_0 + \alpha_1 u_t^2 + \alpha_2 u_{t-1}^2
\end{align}

#### b)

> $Z(x,y)$ is our variable of interest , it depends on the realizations of two random variables $x,y \sim f(x,y)$.
> $\mathbb{E} Z = \int \int  Z(x,y) f(x,y) \mathrm{d}y \mathrm{d} x$
> Bayes law says that $P(A|B) = \frac{P(A\cap B)}{P(B)}$ or equivalently $P(A\cap B) = P(A|B) \cdot P(B)$ which in density terms means $f(x,y) = f(y|x) \cdot f(x)$.
> The expectation therefore can be rewritten as:
> $\mathbb{E} Z = \int \int  Z(x,y) f(y|x)f(x) \mathrm{d}y \mathrm{d} x = \int \underbrace{\int  Z(x,y) f(y|x) \mathrm{d}y}_{\mathbb{E}(Z(y,x)|x)} f(x) \mathrm{d} x = \mathbb{E}(\mathbb{E}(Z|x))$
> $$\mathbb{E} Z = \mathbb{E}(\mathbb{E}(Z|x))$$
> $$Var Z = \mathbb{E} (Z-\mathbb{E} Z)^2 = \mathbb{E} \left(\mathbb{E}((Z-\mathbb{E} Z)^2|x)\right) = \mathbb{E} Var(Z|x)$$

\begin{align}
Var(Y_t) &= Var(b_0 + b_1 Y_{t-1} + u_t) \\
&= b_1^2 Var(Y_{t-1}) + Var(u_{t}) \quad \text{$Y_t$ assumed stationary, and $Var(Y_t) = Var(Y_{t-1})$}\\
&= \frac{1}{1 - b_1^2} Var(u_t)\\
&= \frac{1}{1 - b_1^2} Var(\sigma_t \nu_t)\\
&= \frac{1}{1 - b_1^2} Var(\sqrt{\alpha_0 + \alpha_1 u_{t-1}^2 + \alpha_2 u_{t-2}^2} \nu_t)\\
& = \frac{1}{1 - b_1^2} \mathbb{E} \left(Var_{t-1} (\sigma_t \nu_t)\right) \\
& = \frac{1}{1 - b_1^2} \mathbb{E} \left( \alpha_0 + \alpha_1 u_{t-1}^2 + \alpha_2 u_{t-2}^2\right) \\
&= \frac{1}{1 - b_1^2}(\alpha_0 + \alpha_1 \mathbb{E} u_{t-1}^2 + \alpha_2 \mathbb{E} u_{t-2}^2)\\
& = \frac{1}{1 - b_1^2} (\alpha_0 + \alpha_1  \sigma_u^2 + \alpha_2 \sigma_u^2)
\end{align}

> $$Var(Z) = \mathbb{E} (Z - \mathbb{E} Z)^2 = \mathbb{E} (Z^2) - (\mathbb{E} Z)^2$$

To find the unconditional variance of $u$ we can use the ARCH formula:

$$\sigma_u^2 = \mathbb{E} \sigma_t^2 = \alpha_0 + \alpha_1 \mathbb{E} u_{t-1}^2 + \alpha_2 \mathbb{E} u_{t-2}^2 = \alpha_0 + \alpha_1 \sigma_u^2  + \alpha_2 \sigma_u^2,$$

i.e. from the equation of the ARCH component we get $$\sigma_u^2 = \frac{\alpha_0}{1 - \alpha_1 - \alpha_2}$$

and the unconditional variance of $Y$ is thus

$$Var(Y_t) =\frac{1}{1 - b_1^2} (\alpha_0 + \frac{(\alpha_1 + \alpha_2)\alpha_0}{1 - \alpha_1 - \alpha_2}) = \frac{1}{1 - b_1^2} \frac{\alpha_0}{1 - \alpha_1 - \alpha_2}$$

### Problem 2

#### a)

ARMA(p,q) takes form:

$$y_t = \alpha_1 y_{t-1} + \epsilon_t + \beta_1 \epsilon_{t-1}$$

Our formula looks like

$$\sigma_t^2 = \alpha_0 +\alpha_1 u_{t-1}^2 + \beta_1 \sigma_{t-1}^2 + \beta_2 \sigma^2_{t-2}$$

Let $\eta_t = u_t^2 - \sigma_t^2$. $\eta_t$ is white noise as $u_t$ has zero mean, and we can sustitute away $\sigma_t^2$:

$$u_t^2 - \eta_t = \alpha_0 +\alpha_1 u_{t-1}^2 + \beta_1 (u_{t-1}^2 + \eta_{t-1}) + \beta_2 (u^2_{t-2} + \eta_{t-2})$$

$$u_t^2 = \alpha_0 + (\alpha_1 + \beta_1) u^2_{t-1} + \beta_2 u_{t-2}^2 + \eta_t + \beta_1 \eta_{t-1} + \beta_2 \eta_{t-2},$$

thus $u_t^2$ is an ARMA(2,2) process.

#### b)

Necessary and sufficient condition for stationarity of $u_t^2$ is that the roots of the characteristic polynomial are outside the unit circle.

The characteristic polynomial looks like:

$$1 - (\alpha_1 + \beta_1)L - \beta_2 L^2$$

Its roots are $z_1,z_2 = \frac{\alpha_1 + \beta_1 \pm \sqrt{(\alpha_1 + \beta_1)^2 + 4 \beta_2}}{-2\beta_2}$ they must begreater than 1 in absolute value.

#### c)

$$Var_t(Y_{t+2}) = Var_t(\mu + u_{t+2}) = Var_t(u_{t+2})= Var_t(\sigma_{t+2}\nu_{t+2})$$

$$Var_t(\sigma_{t+2}\nu_{t+2}) = \mathbb{E} (\sigma^2_{t+2}\nu^2_{t+2})- \left(\mathbb{E} (\sigma_{t+2}\nu_{t+2})\right)^2$$

$$\sigma_{t+2}^2$$ is independent of $\nu_{t+2}^2$ therefore the first expectation is equal to the product of expectations:

$$\mathbb{E}_t \sigma^2_{t+2} = \mathbb{E}_t (\alpha_0 + \alpha_1 u^2_{t+1} + \beta_1 \sigma^2_{t+1} + \beta_2 \sigma^2_t)\\
=\alpha_0 + \alpha_1 \mathbb{E}u^2_{t+1} + \beta_1 \mathbb{E}_t \sigma^2_{t+1} + \beta_2 \sigma^2_t\\
=\alpha_0 + (\alpha_1  + \beta_1 )\sigma_{t+1}^2 + \beta_2 \sigma^2_{t},$$

as $\sigma_{t+1}$ depends only on the information we already have by period $t$.

---
title: interpolated Historical simulation VaR
---

## Interpolated Historical Simulation VaR

You are given a sample $\{r_1, \dots, r_T\}$.

1. Construct the empirical CDF
    a. sort $\{r_t\}$ in increasing order. they break the domain into $T+1$ interval. Sorted indexes are $t_1, t_2, t_3 \dots$ : 
    $$t_1 = \arg\min_t\{r_t\}, \quad t_2 = \arg\min_t\{r_t \text{ excluding } r_{t_1}\},$$
    
    b. to each interval (excluding the leftmost one) you assign pdf $\frac{1}{T}$
    c. CDF for each interval is obtained by summing up all pdfs of all the preceding intervals (including the current)
    
    $$\hat F_{HS} (r) = 
    \begin{cases}
     0 & \quad \text{ if $r < r_{t_1}$}\\
     (\hat F_1 =) \frac{1}{T} & \quad \text{ if $ r_{t_1} \leq r < r_{t_2}$}\\
      (\hat F_2 =) \frac{2}{T} & \quad \text{ if $ r_{t_2} \leq r < r_{t_3}$}\\
     \dots \\
     (\hat F_T =) 1 & \quad \text{ if $ r_{t_T} \leq r$}
    \end{cases}$$
    
    This empirical CDF has a "stairs" graph
    
2. Find the interval $r_{t_i}, r_{t_{i+1}}$ such that $\hat F_{HS}(r) > \alpha$ and the smallest of all such intevals. Note down $r_{t_i}, r_{t_{i+1}}, \hat F_{i}, \hat {F}_{i-1}$.

3. Interpolate $\hat F^{-1}(\alpha)$ by using the proportiality of the triangles:
    - the bigger one : with base $r_{t_{i+1}} - r_{t_i}$ and with height $\hat F_i - \hat F_{i-1}$ (which is equal to $\frac{1}{T}$ for historical simulation approach).
    - the smaller one: with base $-VaR - r_{t_i}$ and height $\alpha - \hat F_{i-1}$.
    
    $$\frac{\alpha - \hat F_{i-1}}{\hat F_i - \hat F_{i-1}} =\frac{\alpha - \frac{i-1}{T}}{\frac{1}{T}}= \frac{- VaR - r_{t_i}}{r_{t_{i+1}} - r_{t_i}} \quad \Rightarrow \\
    \quad \frac{1}{T}VaR = -r_{t_i} \frac{1}{T} - (\alpha - \frac{i-1}{T})(r_{t_{i+1}} - r_{t_i}).$$
    
    $$VaR = T\left( -r_{t_i} (\frac{1}{T} - \alpha + \frac{i-1}{T} ) - (\alpha - \frac{i-1}{T})r_{t_{i+1}} \right) \\= -r_{t_i} (i - \alpha) - (\alpha T - (i-1))r_{t_{i+1}}$$
    

### Weighted Historical Simulation

This approach is different in the approximation of CDF:

- HS : pdf was given by $\frac{1}{T}$ for all intervals except the leftmost one.
- WHS : pdf is given by $\frac{1}{T} \times \frac{\lambda^{t_i}}{\sum_{j=0, T} \lambda^j}$

Construction of empirical CDF:
a. Sort the observations -> $i = 1, \dots T$, and $t_i$ is the date of the $i$-th smallest observation, and $r_{t_i}$ - $i$-th smallest observation itself. 
b. when assigning pdf probabilities we diminish the importance of the interval by a factor increasing the in the historical order number $t_i$:

--- 
### h3
A single project net value cdf:

$$
F(x) = P(X \leq x) = 
\begin{cases}
 0 & \quad \text{ if $x <0$}\\
 .04 & \quad \text{ if $0 \leq x < 1000$}\\
 1  & \quad \text{ if $1000 \leq x$}
\end{cases}
$$

A portfolio CDF:
$$
F(x) = P(X \leq x) = 
\begin{cases}
 0 & \quad \text{ if $x <-2000$}\\
 0.04^2 = .0016 & \quad \text{ if $-2000 \leq x <-1000$}\\
 .0016 + 2 \times .04\times .96 = 0.0784 & \quad \text{ if $-1000 \leq x < 0$}\\
 1  & \quad \text{ if $0 \leq x$}
\end{cases}
$$

$$
P(P_1+P_2 = x) = 
\begin{cases}
    .000064 &\quad \text{ if $x = -2$}\\
    .0768 \times .04, &\quad \text{ if $x = -1.5$}\\
    .04 \times .04 = , &\quad \text{ if $x = -.5$}\\
\end{cases}
$$


---

Friedman lifecycle model

$$E\left(\underbrace{\beta\frac{u^\prime(c_{t+1})}{u^\prime(c_t)}}_{SDF}(1+r_{t+1}) \right) = 1$$

$$E(m_{t+1}R_{t+1}) = 1$$

If we have $n$ scenarios with probabilities $p_k$:

$$\sum_{k} p_k m_{t+1}^k R_{t+1}^k = 1$$


https://en.wikipedia.org/wiki/VIX

- interest rate
- exchange rate

---

Question : stability of supply chains in consumer electronics industry

- factors that affect stability (main risks)
    + market risk
        + interest rate
        + currency risk
        + commodity risk
    + counterparty/competitor risks
    + political risks
    https://www.yourarticlelibrary.com/politics/10-different-ways-in-which-political-risk-can-be-managed-investment/5772
    + rare event (epidemics)
    https://www.tradetalkspodcast.com/
    https://www.tradetalkspodcast.com/podcast/151-container-shipping-costs-are-through-the-roof-whos-paying/
    

$$
\text{Profit} = \text{Sales} - \text{Costs}
$$

Consumer electronics supply chain:
    - commodity:
        - aluminium (intensive and rare)
        - other metals
    - Production of intermediate inputs
        - production of inputs (semiconductors)
        - shipping
        - assembly
        - shipping retail
            - inventory of final goods
                - apple : centralized stock in US - more exposed to shipping risks at this stage

https://en.wikipedia.org/wiki/2020%E2%80%93present_global_chip_shortage


---

## Christophersen test

Assume you have a model that gives $VaR^\alpha_{t+\delta}$ based on a sample of observations $\{X_{t-T+1},\dots, x_t\}$.

Assume you have data $\{x_{l}, x_{u}\}$.

Given these data and the VaR model you can construct a series of VaR thresholds:

$$\{VaR_{l+T+\delta}, \dots, VaR_{u}\}$$

0. Construct the (synthetic) sequence of hit-miss measurements:

$$\{Hit_t\}_{t=l+T+\delta}^{u} = \left\{ \begin{cases} 0, \text{ if $-VaR_t < r_t$}\\ 1 , \text{ otherwise} \end{cases} \right\}_{t=l+T+\delta}^{u}.$$

1. If our VaR model is correct, Hit sequence is iid and 0's happen with $\alpha$ probability on average.

$$\text{H}_0 : Hit \sim Bernoulli(\alpha).$$

Let $n = u - (l+T+\delta)$ be the number of VaR's that we have computed.

- First approach : simple t-test (test of population proportion):
    a. Compute the sample average : $\bar{Hit} = \frac{1}{n}\sum_{t=l+T+\delta}^u Hit_t$
    b. If $Hit_t \sim Bernoulli(\alpha)$ then 
    $$E(\bar{Hit}) = \frac{1}{n}\sum E(Hit_t) = \frac{1}{n} \cdot n\cdot (0 \cdot (1-\alpha) + 1 \cdot \alpha) = \alpha$$
    $$Var(\bar{Hit}) = \frac{1}{n^2} \cdot Var(\sum Hit_t) = \frac{1}{n^2} \cdot \sum Var(Hit_t) \\
    = \frac{1}{n^2} \cdot n \cdot (E(Hit^2) - E(Hit)^2) = \frac{1}{n^2}\cdot n (0^2 \cdot (1-\alpha) + 1^2 \cdot \alpha - \alpha^2) = \frac{1}{n}(\alpha - \alpha^2) \\
    = \frac{\alpha(1-\alpha)}{n}$$
    
    Using CLT we can write that
    
    $$\text{t-statistic} = \frac{\bar{Hit} - \alpha}{\sqrt{\alpha(1-\alpha)/n}} \sim N(0,1).$$
    
    Confidence level for this test can generally be different from $\alpha$ (this $\alpha$ is the VaR percentage of worst scenarios).
    
- Second approach : test the same null hypothesis against the best maximum likelihood iid Bernoulli model:

    a.If the alternative hypothesis

    $$\text{H}_a : Hit \sim Bernoulli(\pi),$$

    is correct, the unknown $\pi$ can be estimated by maximizing the likelihood of the $\{Hit_t\}$ sample.

    Likelihood = Probability of observing the given sample under the assumed distribution of the underlying random variable.

    If $H_a$ is correct then

    $$P(Hit_t = 1) = \pi, P(Hit_t = 0) = 1-\pi, \\
\\
L(\{Hit_t\}; \pi) = P(H_{1} = 1, H_2 = 0, H_3 = 1, \dots) \\= P(H_1 = 1) \cdot P(H_2 = 0)\cdot P(H_3= 1) \cdot \dots \\
= \pi \cdot (1-\pi ) \cdot \pi \cdot \dots \\
=\pi^{\# \text{ of Hit = 1}} \cdot (1 - \pi)^{\# \text{of Hit = 0}}$$

    Maximize $L(\{Hit_t\}; \pi)$ with respect to $\pi$:

    $$\frac{\partial \log L}{\partial \pi} = 0 \quad \iff \frac{\# \text{ of Hit = 1}}{\pi} - \frac{\# \text{ of Hit = 0}}{1 - \pi} = 0 \\ \Rightarrow \hat \pi = \frac{\# \text{ of Hit = 1}}{n}$$

    Denote $n_1 = \# \text{ of Hit = 1}, n_0 = \# \text{ of Hit = 0}$. Note that $n = n_0 + n_1$. Then $\hat \pi = \frac{n_1}{n}$.

    The corresponding likelihood takes form

$$L_a = \frac{n_1^{n_1}}{n^{n_1}} \frac{n_0^{n_0}}{n^{n_0}} = \frac{n_1^{n_1} n_0^{n_0}}{n^n}.$$
    
Likelihood for the $H_0$ is given by:

$$L(\{Hit_t\}, \alpha) = \alpha ^{n_1} \cdot (1 - \alpha)^{n_0}$$

Plug the two likelihoods into the likelihood ratio test:

$$LR = -2 \ln (\frac{L(\dots, \alpha)}{L(\dots, \hat \pi)}) \sim \chi^2(1).$$

If the LR statistic is higher than the critical value, the difference between the models is significant and the estimated model is preferred.


> If the simulated LR is frequently higher that the data LR (high P-value) then we can say that the data LR is not too big, and H0 is not rejected.

If $Hit_t$ is first order Markov chain we have 

$$P(Hit_{t+1} = 1 | Hit_t = 0) = \pi_{01}, \quad P(Hit_{t+1} = 1 | Hit_t = 1) = \pi_{11}$$

Let our sample is $\{\dots ,0,1,1,0,1\}$:

> Bayes law $$P(A\cap B) = P(A|B)\cdot P(B)$$

$$L(\{Hit_t\}_{t=1}^T; \pi_{01}, \pi_{11}) \\
= P(Hit_T=1 | Hit_{T-1} = 0, Hit_{T-2} = 1, \dots )\cdot P(Hit_{T-1} = 0, Hit_{T-2} = 1, \dots)\\
= P(Hit_T=1 | Hit_{T-1} = 0)\cdot P(Hit_{T-1}=0 | Hit_{T-2} = 1)\cdot P(Hit_{T-2} = 1, \dots)\\
= \pi_{01}^{n_{01}}\cdot (1-\pi_{01})^{n_0 - n_{01}} \cdot \pi_{11}^{n_{11}}\cdot(1 - \pi_{11})^{n_1 - n_{11}}$$

Try to write down the likelihood of a sample for a second order markov chain:

$$P(Hit_{t+1} = 1 | Hit_t = 0,Hit_{t-1} = 0) = \pi_{001},\\P(Hit_{t+1} = 1 | Hit_t = 0,Hit_{t-1} = 1) = \pi_{101} \\
P(Hit_{t+1} = 1 | Hit_t = 1,Hit_{t-1} = 0) = \pi_{011}
\\P(Hit_{t+1} = 1 | Hit_t = 1,Hit_{t-1} = 1) = \pi_{111}$$.


## Midterm practice

### 

![](https://i.imgur.com/51ZgauX.png)

The characteristic polynomial of the AR part is 

$$Y_t = X + (X-1)Y_{t-1} + (.5X - .25 X^2)Y_{t-2} + \epsilon_t - .25 X^2 \epsilon_{t-1}$$

$$\alpha (L) = 1 - (X-1) L - (.5 X - .25X^2)L^2$$

Its roots are 

$$z_1 = \frac{X-1 - 1}{-2 (.5X - .25 X^2)} = \frac{X - 2}{.5 X (X-2)} = \frac{2}{X}\quad \text{if $X \neq 2$}.$$


$$z_2 = \frac{X-1 + 1}{-2 (.5X - .25 X^2)} = \frac{X}{.5 X (X-2)} = \frac{2}{X-2} \quad \text{if $X \neq 0$}$$

The AR is stationary if $z_1, z_2 $ are outside the unit circle.

$$
\begin{cases}
-2 < X < 2\\
-2 < X-2 < 2
\end{cases}
\iff 
\begin{cases}
-2 < X < 2\\
0 < X < 4 
\end{cases}
\iff
0< X < 2
$$

Invertibility is the property of the MA part. MA(1) is invertible if its coefficient is less then 1 in absolute value:

$$.25X^2 < 1 \iff  x\in (-2,2).$$

The process $y_t$ is both invertible and stationary for $X \in (0, 2)$.

#### b)

If $X \sim N(0, 1.0204^2),$ probability that $y_t$ is invertible and stationary is

$$\mathbb{P}(0 < X < 2) = \mathbb{P}(X < 2) - \mathbb{P}(X < 0) \\
= \mathbb{P}((X - 0)/1.0204 < 2/1.0204) - \mathbb{P}((X - 0)/1.0204 < 0) \\
= \Phi(2/1.0204) - \Phi(0)\\ 
= \Phi(1.96) - \frac{1}{2}\\
= .975 - .5 = .475$$

### c)

$$z_1, z_2 = \frac{X - 1 \pm 1}{.5 X (X - 2)}$$

if $X = 1$, $z_1,z_2 = \pm 2$, both are outside the unit circle, there are no unit roots, so $p=2$, $d=0$, $q=1$.

### d)

$$z_1, z_2 = \frac{X - 1 \pm 1}{.5 X (X - 2)}$$

As $x=2$ 

$$\alpha(L) = 1 - L,$$

only has 1 root equal to 1, therefore $p = 0$, $d = 1$, $q = 1$.

We have 
$$y_t = 2 + y_{t-1} + \epsilon_t - \epsilon_{t-1}\\
=2 + 2 + y_{t-2} + \epsilon_{t-1} - \epsilon_{t-2} + \epsilon_t - \epsilon_{t-1} \\
= 2\cdot 2 + y_{t-2} + \epsilon_t - \epsilon_{t-2} \\
= 2t + y_0 + \epsilon_{t} - \epsilon_{0}$$

$$\mathbb{E}(y_t) = \mathbb{E}(\mathbb{E}_0 y_t) = \mathbb{E}(2t + y_0 - \epsilon_0) =2t + \mathbb{E}(y_0 - \epsilon_0),$$

conditional expectation can be computed, unconditional expectation does not exist (as it has to depend on $t$).

$$Var(y_t) = Var(2t + y_0 + \epsilon_t - \epsilon_0) = Var(\epsilon_t) + Var(\epsilon_0) = 2 \sigma^2.$$

Alternatively, suppose that the unconditional expectation exists. Then, $E y_t = E y_{t-1}$ taking expectaion of both sides of the expression for $y_t$ we get:

$$E y_t = 2 + E y_{t} + 0 - 0 \iff 0 = 2.$$

We get a contradiction, implying that the premise is false.

### e)

$$\mathbb{E} y_t = X - (X-1) \mathbb{E}y_t +.25 X (2-X) \mathbb{E} y_t + 0 + 0 \Rightarrow \mathbb{E} y_t = \frac{X}{X-1 -.25 X (2-X)}.$$

## h2

![](https://i.imgur.com/3U8OVCN.png)

### a) is $Z_t$ homoscedastic?

Yes, as $\epsilon_t$ and $\nu_t$ are iid and cross-independent.

$$Z_t = u_t + \epsilon_t + \beta \epsilon_{t-1}$$

$$Var_t(Z_{t+1}) = Var_t(u_{t+1}) + Var_t(\epsilon_{t+1}) + \beta^2 Var_t(\epsilon_{t}) = 1 + \sigma^2$$

$$Var(Z_{t}) = Var(u_{t}) + Var_t(\epsilon_{t+1}) + \beta^2 Var_t(\epsilon_{t}) = 1 + \sigma^2$$


### b)

$$\gamma_0 = Var(Y_t) = \sigma_\nu^2 + \theta^2 \sigma_\nu^2$$
$$\gamma_1 = Cov(Y_t, Y_{t-1}) = \theta \sigma_\nu^2$$
$$\gamma_s = 0 \quad \text{for $s\geq 2$}$$

Autocorrelation is given by 

$$\rho_1 = \frac{\gamma_1}{\gamma_0} = \frac{\theta}{1 + \theta^2},$$

and $\rho_s = 0$ for $s > 1$.

### c)

$$Z_t = u_t + \epsilon_t + \beta \epsilon_{t-1}$$

$$Y_t = \nu_t + \theta \nu_{t-1}$$

If $Y_t$ and $Z_t$ are identical, then for practical purposes, their moments (expectation, variance, covariances, autocorrelations) must be the same.

$$Var(Z_t) = \sigma_u^2 + \sigma^2_\epsilon + \beta^2 \sigma_\epsilon^2$$
$$Cov(Z_t, Z_{t-1}) = \beta \sigma^2_\epsilon$$

Equality of moments suggests:

$$\mathbb{E} Z_t = \mathbb{E} Y_t = 0,$$
$$\mathbb{Var} Z_t = \mathbb{Var} Y_t \iff \sigma_u^2 + \sigma^2_\epsilon + \beta^2 \sigma_\epsilon^2 = \sigma_\nu^2 + \theta^2 \sigma_\nu^2$$
$$\mathbb{Cov} (Z_t, Z_{t-1}) = \mathbb{Cov} (Y_t, Y_{t-1}) \iff \beta \sigma_\epsilon^2 = \theta \sigma_\nu^2,$$

$$\begin{cases}
\sigma_\nu^2( 1 + \theta^2) &= \sigma_u^2 + \sigma^2_\epsilon + \beta^2 \sigma_\epsilon^2\\
\theta \sigma_\nu^2 &= \beta \sigma_\epsilon^2
\end{cases}$$

### d)

$\beta = 1$, $\epsilon_t=.5 u_t$.

$Z_t$ is still an MA process, which is weakly stationary.

### e)

If $\epsilon_t = u_t^2$.

use the definition of weak stationarity:

$$E Z_t = E(u_t + u_t^2 + \beta u_{t-1}^2) = 0 + Var(u_t) + \beta Var(u_t) = 1 + \beta$$

$$Var Z_t = Var(u_t + u_t^2 + \beta u_{t-1}^2) \\= Var(u_t) + Var(u^2_t) + \beta^2 Var(u^2_{t-1}) \\+ 2 Cov(u_t,u_t^2) + 2\beta Cov(u_t, u_{t-1}^2) + 2 \beta Cov(u^2_t, u_{t-1}^2)\\
= 1 + (E(u_t^4) - (E u_t^2)^2)(1 + \beta^2) + 2 (E(u_t^3) - Eu_t \cdot Eu_t^2) \\
= 1 + (1 + \beta^2)(E(u_t^4) - 1) + 2 Eu_t^3 
$$

$$Cov(Z_t,Z_{t-1}) = Cov(u_t + u_t^2 + \beta u_{t-1}^2, u_{t-1} + u_{t-1}^2 + \beta u_{t-2}^2) \\=  Cov(\beta u_{t-1}^2, u_{t-1} + u_{t-1}^2) \\
= \beta Cov(u_{t-1}^2,u_{t-1})+$$

$$Cov(Z_t,Z_{t-2}) = Cov(u_t + u_t^2 + \beta u_{t-1}^2, u_{t-2} + u_{t-2}^2 + \beta u_{t-3}^2) = 0$$


$Z_t$ remains weakly stationary.


## Q.4

1. Realized volatility measures:
- If we want to compute *daily* realized volatility, we need a higher frequency dataset (say hourly or 5-minutes):

$$\{p_{00:00}, p_{00:05}, p_{00:10}, \dots, p_{23:55}, p_{24:00}\}$$

Calculate the 5-minute returns:

$$r_{00:05} = \frac{p_{00:05} -p_{00:00}}{p_{00:00}}$$

$$RV = \sum_{\tau = 00:05}^{24:00} r^2_\tau$$

- OHLC - "open, high, low, close" - based estimators of realized volatility

    - Roger and Satchell volatility
    - Garman Klass volatility

---

## Value At Risk
  
$W_1, \dots W_p$ - amounts in USD of positions in assets $1,\dots p$

$r^1_t, \dots r^p_t$ - returns on assets $1, \dots p$

The whole portfolio will earn:

$$Profit_t = W_1 r^1_t + \dots W_p r^p_t \\ =W (w_1 r_t^1 + \dots + w_p r^p_t)\\
r_t = \frac{Profit_t}{W} = w_1 r_t^1 + \dots + w_p r^p_t$$

In the first semester we studied a lot of models or returns all of which implied normal distribution of future returns $r^1, \dots r^p$ given todays information.

$$r^k_t \sim N(\mu^k_t, \sigma_t^{k2}) $$

This means that $r_t$ given today's information is distributed as a normal as well.

$$r_t \sim N(w_1 \mu_t^1 + \dots + w_p \mu_t^p, w_1^2 \sigma^{12}_t +\dots + w_p^2 \sigma^{p2}_t + \text{Covarinaces})$$

![](https://i.imgur.com/XJhYQLN.png)

#### Definition

$VaR(5\%)$ - **the level of return on the portfolio, such that in 95% of cases the return is higher than this value.**

$$\mathbb{P}(r_t < VaR(5\%)) = 5\%$$

It is also possible to compute $VaR(1\%), Var(10\%)$ or for any critical loss probability really.

VaR is equivalent conceptually to a critical value of some statistical test.

If $r_t \sim N(\mu, \sigma^2)$, then the VaR can be specified in the following way:

$$\mathbb{P}(r_t < VaR_{5\%}) = \mathbb{P}(\frac{r_t - \mu_t}{\sigma_t} < \frac{VaR_{5\%} - \mu_t}{\sigma_t}) = \Phi(\frac{VaR_{5\%} - \mu_t}{\sigma_t}) = 5\%$$

where $\Phi(x) = \mathbb{P}_{X \sim N(0,1)}(X < x)$.

$$\frac{VaR_{5\%} - \mu_t}{\sigma_t} = \Phi^{-1}(5\%)\\
$$

$\Phi^{-1}(5\%) = - \Phi^{-1}(95\%)$ because the standard normal distribution is symmetric. $\Phi^{-1}(95\%)$ is simply the 95% critical value that is usually looked up in a table during an exam.

$$\Phi^{-1}(95\%) = 1.65\quad \Phi^{-1}(97.5\%) = 1.96 \quad \Phi^{-1}(99\%) \approx 2$$$

$$VaR_{5\%} = \mu_t + \sigma_t \times \Phi^{-1}(5\%)$$


![](https://i.imgur.com/h2CWm3Y.png)

### Revision exercise

![](https://i.imgur.com/eLY5Unk.png)

$$\gamma_0 = (1 + \theta_1^2 + \theta_2^2)\sigma^2$$

$$\gamma_1 = (\theta_1 + \theta_1 \theta_2)\sigma^2, \gamma_2= \theta_2 \sigma^2$$


#### 

![](https://i.imgur.com/L3GBFO9.png)

**Auxiliary question:**

GARCH is a model of the *conditional variance* of the returns. What is the *unconditional* (long-term) variance of the returns?

$$\sigma^2 = Var(r_t) = Var(\epsilon_t) = E(\epsilon_t^2) = E(E_t (\epsilon_{t+1}^2)) = E(\sigma^2_{t+1})$$

When we speak of the unconditional (long-term) moments of the distribution
$$E(\epsilon^2_t) = E(\epsilon^2_{t\pm k})$$
$$E(\sigma^2_t) = E(\sigma^2_{t+1})$$

$$\sigma_{t+1}^2 - \beta \sigma^2_t = \omega + \alpha \epsilon_t^2 $$

Take expectaiton of both sides 
$$(1 - \beta) E(\sigma^2_t) = \omega + \alpha E(\epsilon_t^2) \\ E(\sigma^2_t)(1 - \alpha - \beta) = \omega \\ Var(\epsilon_t) = \frac{\omega}{1 - \alpha - \beta}
$$

#### Main question:

We are asked to write down:

$$\sigma^2_{t+1} = f(\sigma^2, \sigma^2_t - \sigma^2, \epsilon^2_t - \sigma^2)$$

We have 

$$\sigma^2_{t+1} = \omega + \beta \sigma^2_t + \alpha\epsilon^2_t \to f(\sigma^2_t, \epsilon^2_t)\\
\sigma^2_{t+1} = \omega + \beta (\sigma^2_t - \sigma^2 + \sigma^2) + \alpha(\epsilon^2_t - \sigma^2 + \sigma^2)\\
\sigma^2_{t+1} = \omega + (\alpha + \beta)\sigma^2 + \beta (\sigma^2_t - \sigma^2) + \alpha(\epsilon^2_t - \sigma^2)
$$

b) 

$\sigma_t^2$ is the variance of $\epsilon_t$ as expected from period $t-1$, i.e. $\sigma^2_t = E_{t-1}{\epsilon_t^2}$. 

At moment $t$ we observe $\epsilon_t$, we remember $\sigma_t$ calculated in $t-1$ and we can calculate $\sigma^2_{t+1}$ for the next period. 
At moment $t+1$ we learn $\epsilon_{t+1}$ and can additionally compute $\sigma^2_{t+2}$.

Unconditional expectation (long-term, zero-information expectation, formaed before the history started) of $\sigma_t^2$ and $\epsilon_t^2$ and $\sigma_{t+1}^2$ are the same.

We need to compute the 2step ahead forecast of variance:

$$E_{t}(\sigma^2_{t+2}) = E_{t}(\omega + \beta \sigma^2_{t+1} + \alpha \epsilon^2_{t+1}) = \omega + \beta \sigma^2_{t+1} + \alpha E_t(\epsilon^2_{t+1}) = \omega + \beta \sigma^2_{t+1} + \alpha \sigma^2_{t+1} \\= \omega + (\alpha + \beta)(\omega + \alpha \epsilon^2_t + \beta \sigma^2_t) $$

c) 

$$E_t(\sigma^2_{t+3}) = \omega + \beta E_t(\sigma^2_{t+2}) + \alpha E_t(\epsilon_{t+2}^2)\\ = \omega + \beta E_t(\sigma^2_{t+2}) + \alpha E_t(\underbrace{E_{t+1}(\epsilon_{t+2}^2)}_{\sigma_{t+2}^2}) \\= \omega + (\alpha + \beta) (\omega + (\alpha + \beta)(\omega + \alpha \epsilon^2_t + \beta \sigma^2_t))$$

$$\sigma_{t+3,t}^2 = (1 + (\alpha + \beta) + (\alpha + \beta)^2)\omega + (\alpha + \beta)^2(\alpha \epsilon_t^2 + \beta \sigma^2_t)$$

$$\sigma_{t+h,t}^2 = \omega \sum_{k = 0}^{h} (\alpha + \beta)^k + (\alpha + \beta)^{h-1}(\alpha \epsilon_t^2 + \beta \sigma^2_t)=
\\=\frac{1 - (\alpha + \beta)^{h}}{1 - (\alpha + \beta)}\omega + (\alpha + \beta)^{h-1}(\alpha \epsilon_t^2 + \beta \sigma^2_t)$$

> $1 + \rho + \dots + \rho^t = \frac{1 - \rho^{t+1}}{1 - \rho}$

#### Example question: VaR

a) 

By definition of VaR:

$$P(R_{t+1} < VaR) = P(\sigma_{t+1} \nu_{t+1} < VaR) = P(\nu_{t+1} < \frac{VaR}{\sigma_{t+1}}) = \Phi(\frac{VaR}{\sigma_{t+1}}) = \alpha$$
$$VaR_{t+1} = \sigma_{t+1}\Phi^{-1}(\alpha)$$

b) Expected shortfall:

Mean return conditional on the sceanrio realized being one of the $(1-\alpha)$% most favourable: 

$$ES_{t+1}^\alpha = E(R_{t+1}| R_{t+1} < VaR)$$

Plug in the model for $R$:

$$ES = E(\sigma \nu |\frac{R}{\sigma} < \frac{VaR}{\sigma}) =\sigma E( \nu |\nu < \frac{VaR}{\sigma})\\
= \sigma(- \frac{\phi(\frac{VaR}{\sigma})}{\Phi(\frac{VaR}{\sigma})}) \\
= \sigma \frac{\phi(\frac{VaR}{\sigma})}{\alpha}$$

From part a) we know that $\frac{VaR}{\sigma} = \Phi^{-1}(\alpha)$

$ES_{t+1}^\alpha = \sigma_{t+1} \frac{\phi(\Phi^{-1}(\alpha))}{\alpha}$

> $\Phi(x) = P(X < x)$ for $X \sim N(0,1)$
> $\phi(x) = \Phi(x)^\prime_x$

![](https://i.imgur.com/hJNwQZd.png)


c) 

$$\lim_{\alpha \to 0} \frac{ES - VaR}{VaR} = \\
= \lim_{\alpha \to 0} \frac{\frac{\sigma}{\alpha} \phi(\Phi^{-1}(\alpha)) - \sigma \Phi^{-1}(\alpha)}{\sigma \Phi^{-1}(\alpha)} \quad \text{(let $z\triangleq \Phi^{-1}(\alpha)$) }\\
= \lim_{z \to -\infty} \frac{\phi(z)}{\alpha z} - 1$$

So far $VaR$ was the critical value of returns separating 5% worst scenarios from the 95% best scenarios (it was negative). Let's change the definition for VaR to be the absolute value of that.

Now $VaR > 0$ is defined implicitly by 
$$P(R < -VaR) = \alpha$$
and the expected shortfall $ES > 0$ is defined by 
$$ES = - E(R | R < -Var)$$

In this exercise we get $$VaR = -\sigma \Phi^{-1}(\alpha) \quad ES = \frac{\sigma}{\alpha} \phi(\Phi^{-1}(\alpha))$$

Substitutte these formulas into the limit and get:

$$\lim_{\alpha \to 0} \frac{ES - VaR}{VaR} = \lim_{z\to -\infty \quad \alpha \to 0 } \frac{\frac{\sigma}{\alpha}\phi(z)}{- \sigma z} - 1\\
= -\lim \frac{\partial_\alpha \phi(z)}{\partial_\alpha (\alpha z)} - 1 = $$

> Derivative of a standard normal pdf:
> $$\phi(z) = -\frac{1}{\sqrt {2 \pi}} \exp(- \frac{z^2}{2}))$$ 
> $$\ln \phi(z) = \dots - \frac{z^2}{2}$$
> Taking derivatives wrt $z$ on both sides
> $$\frac{1}{\phi(z)}\phi^\prime(z) = - z \quad \Rightarrow \phi^\prime(z) = - z \phi (z)$$

$$\Phi(z) = \alpha \quad \Rightarrow \phi(z) \partial_\alpha z = 1 \quad \Rightarrow \partial_\alpha z = \frac{1}{\phi(z)}$$

Going back to our limit fomula:

$$\lim \frac{ - z \phi \frac{1}{\phi}}{z + \alpha\frac{1}{\phi}} = \lim \frac{ - z}{z + \alpha\frac{1}{\phi(z)}}$$

$$\lim_{\alpha \to 0} \frac{ES - VaR}{VaR} = \lim \frac{z}{z + \alpha\frac{1}{\phi(z)}} - 1 = 0,$$
as $\alpha / \phi \to 0$.



#### Example VaR:

$$R_1 = \begin{cases}0 \text{ with prob=96%}\\ -1 \text{ with prob=4%} \end{cases}$$

Mean returns are not zero:

$$E(R) = .04 * (-1) + .96 * 0 = -.04$$

$E(P_1) = 1000*(-.04) = -40$

VaR is given by the formula:

$$
P(P_1 < E(P_1) - VaR) = 5\% \iff P(P_1 > E(P_1)-VaR) = 95%
$$


## UoL exam 2019A, question 2

### 1.
Deviation from unconditional variance:

$$\sigma^2_t - \sigma^2$$

$$\sigma^2_{t+1} - \sigma^2 = \omega - \sigma^2 + \beta \sigma^2_t + \alpha \epsilon_t^2 \\
= \omega - \sigma^2 + \beta (\sigma^2_t \pm \sigma^2) + \alpha (\epsilon_t^2 \pm \sigma^2)\\
= \omega - \sigma^2 + (\beta + \alpha)\sigma^2 + \beta (\sigma^2_t - \sigma^2) + \alpha (\epsilon_t^2 - \sigma^2)\\
= \omega + (\beta + \alpha - 1)\sigma^2 + \beta (\sigma^2_t - \sigma^2) + \alpha (\epsilon_t^2 - \sigma^2)$$

### 2.

We know that $\epsilon_{t+2} = \sigma_{t+2}\nu_{t+2},$ where $\nu_{t+2} \sim iid (0,1)$ in heteroskedastic models.

\begin{align}
    \sigma^2_{t+2,t} &= E_t(\epsilon^2_{t+2})\\ 
    &=E_t(\sigma^2_{t+2} \nu^2_{t+2}) \quad \text{ $\nu_{t+2}$ is independent of $\sigma_{t+2}$}\\
    & = E_t(\sigma^2_{t+2})E_t(\nu^2_{t+2})\\
    & = E_t(\omega + \beta\sigma^2_{t+1} + \alpha \epsilon_{t+1}^2)\\
    & = \omega + \beta E_t(\underbrace{\sigma^2_{t+1}}_\text{known ad date $t$}) + \alpha E_t(\underbrace{\epsilon_{t+1}^2}_\text{random at $t$, with variance $\sigma_{t+1}^2$})\\
    & = \omega + (\alpha + \beta)\sigma_{t+1}^2
\end{align}

### 3. 

\begin{align}
\sigma_{t,t+3}^2 &= E_t(\epsilon_{t+3}^2)\\
&= E_t(\sigma_{t+3}^2)\\
&= \omega + \beta E_t(\sigma_{t+2}^2) + \alpha E_t(\epsilon_{t+2}^2)\\
&= \omega + \beta(\omega + (\alpha +\beta) \sigma^2_{t+1}) + \alpha \sigma_{t+2,t}^2 \\
& = \omega + (\alpha + \beta)\sigma^2_{t+2,t}\\
&= \omega(1 + \alpha + \beta) + (\alpha + \beta)^2 \sigma_{t+1}^2
\end{align}

$$\sigma^2_{t+h,t} = \omega \sum_{j=0}^{h-2}(\alpha + \beta)^j + (\alpha + \beta)^{h-1}\sigma^2_{t+1}\\
 = \omega \frac{1 - (\alpha + \beta)^{h-1}}{1 - \alpha - \beta} + (\alpha + \beta)^{h-1}\sigma^2_{t+1}$$
 
> $$(1-x)(1 + x + x^2 + x^3 + \dots + x^{h-1}) \\
= 1 + x + x^2 + x^3 + \dots + x^{h-1} - x - x^2 - x^3 - \dots - x^{h}  \\
= 1 - x^h$$
> $$1 + x + x^2 + x^3 + \dots + x^{h-1} = \frac{1 - x^h}{1 - x}$$

$$\sigma^2_{t+h,t} - \frac{\omega}{1 - \alpha - \beta} = (\alpha + \beta)^{h-1}(\sigma^2_{t+1} - \frac{\omega}{1 - \alpha - \beta})$$

Note that $\sigma^2 = \frac{\omega}{1 - \alpha - \beta}$:

$$E(\sigma_{t+1}^2) = \omega +\beta E(\sigma^2_t) + \alpha E(\epsilon_{t}^2) \Rightarrow \sigma^2 = \omega  + \beta \sigma^2 + \alpha \sigma^2$$

$$\sigma^2_{t+h,t} - \sigma^2 = (\alpha + \beta)^{h-1}(\sigma^2_{t+1} - \sigma^2)$$

$\sigma^2$ is the unconditional expectation, i.e. forecast with zero (null) information. As the forecast horizon $h$ grows to infinity, the conditional forecast of variance approaches this zero-information forecast.

This model can be used for isntance to compute the $h$-step ahead value-at-risk.

For normal returns Value-at-risk(5%) was calculated as
$$VaR = -r_{t+h,t} - \sigma_{t+h,t} \times 1.65,$$

where 
- $r_{t+h,t}$ is the $h$-step ahead forecast of the returns,
- $\sigma_{t+h,t}$ is the $h$-step ahead forecast of volatility.

## Question 2

### 1.

A GARCH process 

$$\sigma^2_{t+1} = \omega + \sum_{j=0}^{q} \beta_j \sigma^2_{t-j} + \sum_{i=0}^p \epsilon^2_{t-i}$$

is covariance stationary if 


$$\sum_i\alpha_i + \sum_j \beta_j < 1.$$


### 2.

Take the uncondtional (zero-information) expectaiton of both sides:

$$E\sigma^2_{t+1} = \omega + \beta E(\sigma^2_t) + \alpha E\epsilon_{t}^2$$

$$\sigma^2 = \frac{\omega}{1 - \alpha - \beta}.$$

### 3. 

The $h$-step ahead volatility forecast is given by the formula: 

$$\sigma^2_{t+h,t} = \omega \sum_{i=1}^{h-2} (\alpha + \beta)^{i} + (\alpha + \beta)^{h-1} \sigma^2_{t+1} $$

> $$\sum_{i=0}^{k-1} \delta^i = \frac{1 - \delta^k}{1 - \delta}$$

The deviation of the $h$-step ahead forecast from the long-term average volatility is given by:

\begin{align}
      \sigma^2_{t+h,t} - \sigma^2 &= \omega \frac{1}{1 - \alpha - \beta} + (\alpha + \beta)^{h-1}(\sigma^2_{t+1} - \sigma^2)  - \sigma^2\\
      &= (\alpha + \beta)^{h-1}(\sigma^2_{t+1} - \sigma^2)
\end{align}

To compute the $h$-step ahead forecast it is sufficient to know 3 things:
- $\sigma^2$ - the long-term average volatility
- $\sigma^2_{t+1} - \sigma^2$ - the deviation of future volatility from the long-term average volatility
- $\alpha+\beta$

$$\sigma^2 = \frac{\omega }{1 - \alpha - \beta} = \frac{.000004}{1 - .93 - .06} = 4 \times 10^{-6} \times 10^2 = .0004$$

$$\sigma_{t+1} = .025, \quad \sigma^2_{t+1} - \sigma^2 = .000625 - .0004 = .000225$$

$$\sigma^2_{t+20,t} = .0004 + .99^{19} \times .000225 = 0.0006231775$$

$$\sigma^2_{t+40,t} = .0004 + .99^{39} \times .000225 = 0.0006236275 $$


---
## Expectations primer:

- **Conditional expectation:** $E_t(\cdot)$ - all variables with subscript $t, t-1, t-k$ under the expectation bracket are treated as constants (knowns), therefore they can be factorized.

- **Unconditional expectation:** $E(\cdot)$ - no variable with index $t, t-1, t\pm k$ can be treated as a known, non-random varaible, but $\epsilon_t$, $\epsilon_{t+1}, \epsilon_{t\pm k}$ **all have the same expectaion** (distribution).


For example:

Take GARCH(1,1) model:

$$\epsilon_t = \sigma_t \cdot \nu_t, \quad \nu_t\sim N(0,1)$$
where
$$\sigma^2_{t+1} = \omega + \beta \sigma^2_t + \alpha \epsilon_t^2.$$

And consider the distribution of $\epsilon_{t+2}$:

$$E(\epsilon_{t+2}) = E(\sigma_{t+2} \nu_{t+2}) = E(\sigma_{t+2})E(\nu_{t+2})= 0,$$

> $E(h(X) g(Y)) = E(h(X)) \cdot E(g(Y)),$ if $X,Y$ are independent.

$$E(\epsilon^2_{t+2}) = E(\sigma_{t+2}^2) E(\nu_{t+2}^2) = E(\sigma_{t+2}^2) = \sigma^2,$$

where $\sigma^2 = E(\sigma_{t}^2) = E(\sigma_{t+1}^2) = E(\sigma_{t+2}^2) = \dots$ is the long-term average volatility (unconditional volatility).

The unconditional volatility must satisfy the GARCH model:

$$\sigma^2 = \omega + (\alpha + \beta) \sigma^2.$$

$$E_t(\epsilon_{t+2}) = E_t(\sigma_{t+2}\nu_{t+2}) = E_t(\sigma_{t+2})E_t(\nu_{t+2})= 0$$
$$E_t(\epsilon_{t+2}) = E_t(\sigma_{t+2}^2)E_t(\nu^2_{t+2}) = E_t(\omega + \beta \sigma_{t+1} + \alpha \epsilon^2_{t+1})\\= \omega + \beta E_t(\sigma^2_{t+1}) + \alpha E_t(\epsilon_{t+1}^2) = \omega + \beta E_t(\omega + \beta\sigma^2_t + \alpha \epsilon_{t}^2) +  \alpha \sigma^2_{t+1}\\
= \omega + \beta(\omega + \beta\sigma^2_t + \alpha \epsilon_{t}^2) +  \alpha \sigma^2_{t+1}\\
=\omega + \beta\sigma^2_{t+1} +  \alpha \sigma^2_{t+1}
$$


---
title: HA 4
---

# Homework 4

## Problem 1


$t=15$
### a)

You are given a sample $\{r_k\}_{k=1}^{15} = \{r_1, \dots, r_{15}\}$.

To calculate $\hat F$:
1. Take several (three-four) smallest values from the sample, order them.
2. For each compute $j = t - k + 1$. let $a_k$ be the position of observation $k$ in the sorted list.
In our case (indices reflect $a_k$)
$$r_0 = -2.95\%, \quad r_1 = -1.76\%, \quad r_2 = -.79\%, \text{ etc.}$$

3. The historical simulation CDF:
$$\hat F_{HS}(r) = \begin{cases} 
0 & \text{ if $r <r_0$}\\ 
1/15 \approx .067 & \text{ if $r_0 \leq r < r_1$}\\
2/15 \approx .133 & \text{ if $r_1 \leq r < r_2$}\\
\text{etc.}
\end{cases}$$

4. Obviously $\alpha = 5\%$ falls into the bin with bounds $[-2.95\%, -1.76\%)$, $\Rightarrow$ 
\begin{align}
\text{Conservative VaR : }\quad &2.95\%,\\
\text{Formal VaR : }\quad &1.76\%,\\
\text{Interpolated VaR : }\quad & (.75 - 0)\times 1.76\% + (1 - .75)\times 2.95\% =2.06\%,
\end{align}

    $\alpha = 10\%$ falls into the bin $[-1.76\%,-.79\%)$ and 
\begin{align}
\text{Conservative VaR : }\quad &1.76\%,\\
\text{Formal VaR : }\quad &.79\%,\\
\text{Interpolated VaR : }\quad & (1.5 - 1)\times 1.76\% + (2 - 1.5)\times .79\% =1.28\%,
\end{align}

6. Expected shortfall for $\alpha = 10\%$ is

\begin{align}
\text{Conservative VaR : }\quad &\frac{1}{2} \times 1.76\% + \frac{1}{2}\times 2.95\% = ,\\
\text{Formal VaR : }\quad &\frac{1}{3}\times .79\% + \frac{1}{3} \times 1.76\% + \frac{1}{3}\times 2.95\% = ,\\
\text{Interpolated VaR : }\quad & ,
\end{align}

## b)

|$i$|${t_i}$|$r_{t_i}$|$\hat F$ from a)|$\hat F^w$ from b)|
|---|---|---|---|---|
|1|11|-2.95|$1/15 \approx .067$| $\lambda^{11}\frac{1 - \lambda}{1 - \lambda^{15}} = .95^{11} \times .05/(1-.95^{15}) \approx .053$|
|2|2|-1.76|$2/15 \approx .133$|$\dots +\lambda^{2}\frac{1 - \lambda}{1 - \lambda^{15}} = .053 + .95^2 \times .05/(1-.95^{15}) \approx .137$|
|3|10|-0.79|$3/15 = .200$|$\dots + \lambda^{10}\frac{1 - \lambda}{1 - \lambda^{15}} = .137 + .95^{10} \times .05/(1-.95^{15}) \approx .193$|


$$

$$\frac{\lambda^{t_i}}{\sum_{k=0}^{m-1} \lambda^k} = \frac{\lambda^{t_i}}{(1 - \lambda^m)/(1-\lambda)} = \frac{\lambda^{t_i} (1 - \lambda)}{1 - \lambda^m}$$

The weighted historical simulation CDF is given by

$$
\hat F^w (x) = \begin{cases}
0, &\quad \text{ if $x < -2.95$}\\
.053, &\quad \text{ if $-2.95\leq x < -1.76$}\\
.137, &\quad \text{ if $-1.76\leq x < -.79$}\\
\text{ etc. }
\end{cases}
$$

For $\alpha = 5\%$, we fall into the bin $[-2.95, -1.76)$:

\begin{align}
\text{Conservative VaR : }\quad &2.95\%,\\
\text{Formal VaR : }\quad &1.76\%,\\
\text{Interpolated VaR : }\quad & \frac{-VaR - (-2.95)}{-1.76 - (-2.95)} = \frac{.05 - .0}{.053 - .0} : VaR = 2.95 - 1.19\times \frac{50}{53} = ,
\end{align}


For $\alpha = 10\%$ we fall into the bin $[-1.76, -.079)$

\begin{align}
\text{Conservative VaR : }\quad &1.76\%,\\
\text{Formal VaR : }\quad &.79\%,\\
\text{Interpolated VaR : }\quad & \frac{-VaR - (-1.76)}{-.79 - (-1.76)} = \frac{.1 - .053}{.137 - .053} : VaR = 1.76 - 0.97\times \dots = ,
\end{align}

Expected shortfall for $\alpha= 10\%$:

We must take the average of the observations that fall beyond the -VaR, **weighted by the new WHS probabilities** ($\frac{\lambda^{t_i}}{\sum \lambda^k}$):

\begin{align}
\text{Conservative ES : }\quad& 2.95, \text{ (single observation beyond the VaR)}\\
\text{Formal ES : }\quad &(2.95 \times .053 + 1.76 \times (.137 - .053))/.137 \approx 2.22,\\
\text{Interpolated ES : }\quad & \frac{-VaR - (-1.76)}{-.79 - (-1.76)} = \frac{.1 - .053}{.137 - .053} : VaR = 1.76 - 0.97\times \dots = ,
\end{align}


---
## Problem 2

The returns are zero mean and normal with conditional variance given by a TGARCH model:

$$\sigma_{t+1}^2 = \begin{cases}
.1 +  .5 \sigma^2_t + .45 r_t^2 &\quad \text{ if $r_t <0 $}\\
.1 +  .5 \sigma^2_t + .25 r_t^2 & \quad\text{ if $r_t \geq0 $}
\end{cases}
$$

$$P(r_{t+1} < \bar r) = P(\frac{r_{t+1}-0}{\sigma_{t+1}} < \frac{\bar r}{\sigma_{t+1}}) = \Phi(\frac{\bar r }{\sigma_{t+1}})$$

If we set $\bar r = - VaR^\alpha$ : 

$$\Phi(\frac{-VaR - \mu_{t+1}}{\sigma_{t+1}}) = \alpha$$

$$VaR^\alpha = - \mu_{t+1} - \sigma_{t+1} \underbrace{\Phi^{-1}(\alpha)}_{\alpha \text{ critical value}}$$

For a normal distribution:

$\Phi^{-1}(.01) = -2.33$
$\Phi^{-1}(.05) = -1.65$


### a)

Conditional VaR for $t+1$ ($r_t = -1.74, \sigma_t^2 =4.85$):

$$\sigma_{t+1}^2 = .1 + .5 \times 4.85 + .45 \times 1.74^2 = 3.31 $$

A short reminder on how to obtain the VaR:

\begin{align}
P(r_{t+1} \leq -VaR^\alpha) &= \alpha
\\
P(\frac{r_{t+1}-\mu_{t+1}}{\sigma_{t+1}} \leq \frac{-VaR^\alpha-\mu_{t+1}}{\sigma_{t+1}}) &=\alpha \\
\Phi(\frac{-VaR^\alpha-\mu_{t+1}}{\sigma_{t+1}}) &= \alpha \\
\frac{-VaR^\alpha-\mu_{t+1}}{\sigma_{t+1}} &= \Phi^{-1}(\alpha)
\end{align}

We get the formula
$$\color{red}{VaR^\alpha = - \mu_{t+1} - \sigma_{t+1} \Phi^{-1}(\alpha)}$$

According to the problem statement $\mu_{t+1}=0$. Then

$$VaR^\alpha = -\mu_{t+1} - \sigma_{t+1} \Phi^{-1}(\alpha) = \begin{cases} 
-\sqrt {3.31} \times (- 2.33) = 4.24 & \quad \text{ for $\alpha = .01$}\\
-\sqrt{3.31} \times(-1.65) = 3.00 & \quad \text{ for $\alpha = .05$}
\end{cases}$$

> Truncated expectation of a normal:
> https://en.wikipedia.org/wiki/Truncated_normal_distribution
> $\beta = \frac{b - \mu}{\sigma}, \alpha_0 = \frac{a - \mu}{\sigma}$
> $$E(X\mid a < X < b) = \mu - \sigma\frac{\phi(\beta) - \phi(\alpha)}{\Phi(\beta) - \Phi(\alpha)}$$


Expected shortfall is given by

$$ES^\alpha = -E(r_{t+1} | r_{t+1} < -VaR^\alpha) = -E(r_{t+1} | -\infty < r_{t+1} < -VaR^\alpha)$$

Knowing that $r_{t+1} \sim N(0, \sigma^2_{t+1})$:

$$X : r_{t+1},\quad a : -\infty,\quad b : -VaR^\alpha $$
$$\alpha_0 \to -\infty,\quad \beta = \frac{-VaR^\alpha}{\sigma_{t+1}} = \frac{\mu_{t+1} +\sigma_{t+1} \Phi^{-1}(\alpha)}{\sigma_{t+1}} = \Phi^{-1}(\alpha)$$
$$\phi(\alpha_0) = 0,\quad  \phi(\beta) =  \phi(\Phi^{-1}(\alpha)) = \frac{1}{\sqrt{2\pi}}\exp(-\frac{1}{2}(\dots)^2)$$
$$\Phi(\alpha_0) =0,\quad \Phi(\beta) = \alpha$$

We get:

$$ES = - \left(\mu_{t+1} - \sigma_{t+1}\frac{\phi(\Phi^{-1}(\alpha)) -0}{\Phi(\Phi^{-1}(\alpha)) - 0}\right) = \sigma_{t+1}\frac{\phi(\Phi^{-1}(\alpha))}{\alpha}$$


---
$$
ES^\alpha = \frac{\sigma_{t+1}}{\alpha} \phi(\Phi^{-1}(\alpha)) = 
\begin{cases}
\frac{\sqrt {3.31}}{.01} \phi(- 2.33) = & \quad \text{ for $\alpha = .01$}\\
\frac{\sqrt {3.31}}{.05} \phi(- 1.65) = & \quad \text{ for $\alpha = .05$}
\end{cases}
$$

> Skewness & Kurtosis
> Mean $E(X)$
> Variance $E(X^2) - E(X)^2$
> Skewness $\propto E(X^3)$ - мера перекоса распределения. If pdf is asymmetric, skewness is different from zero. Skewness is always zero for a normal distribution. $\chi^2(k)$ has a positive skewness $\sqrt{\frac{8}{k}}$.
> Kurtosis $\propto E(X^4)$ - мера толщины хвостов. If the tails are thicker kurtosis is higher. $\chi^2(k)$ has kurtosis $\frac{12}{k}$. Kurtosis of a normal distribution is zero.

### b)

Unconditional VaR requires computing the unconditional variance $\sigma^2$. For the TGARCH model it is found from the equation:

$$\sigma^2 = E(r_{t}^2) = E(\sigma_{t+1}^2) = E(\sigma_t^2)$$

Applying the expectation operator to the TGARCH model we get:

$$\sigma^2 = .1 + .5 \sigma^2 + .25 \times E(r_t^2\mid r_t\geq 0) \cdot P(r_t\geq0) + .45 \times E(r_t^2\mid r_t<0) \cdot P(r_t<0)$$
Note that $r_t$ has a symmetric distribution, centered around zero, meaning that $P(r_t<0) =  P(r_t>0) = \frac{1}{2}$, and $r$ is ditstributed in the same way as $-r$, i.e. $E(r^2|r>0) = E(r^2|r<0)$.
$$E(r_t^2) = E(r_t^2\mid r_t<0) \cdot P(r_t<0) + E(r_t^2\mid r_t\geq 0) \cdot P(r_t\geq0) \\
\Rightarrow \sigma^2 = E(r^2|r<0) = E(r^2|r>0).$$

$$\sigma^2 = .1 + .5 \sigma^2 + .35 \sigma^2 \Rightarrow \sigma^2 = .1/(1 - .5 - .35) = .667$$

$$VaR^{1\%} = -.667 \times (-2.33) = \dots$$


### c)

Probability of breaching $VaR = -2.52$ is:

\begin{align}
P(r_{t+1} < -2.52) &= P(\frac{r_{t+1}}{\sigma_{t+1}} < \frac{-2.52}{\sigma_{t+1}})\\
&=\Phi(\frac{-2.52}{\sigma_{t+1}})\\
&=\Phi(\frac{-2.52}{\sqrt {3.31}})\\
&\approx .083
\end{align}

The historical simulation VaR corresponds to the 8.3% VaR based on the TGARCH model.

## d)

Replace the critical values $-2.33, -1.65$ by the t-Student(6) critical values $-3.14, -1.94$.

## e)

1. $VaR^\alpha_{HS}$ series were constructed based on the sample
2. 


## Problem 3

$$\text{Revenue}_i = 
\begin{cases}
-\$4M &\quad \text{ with prob=.02}\\
-\$1M &\quad \text{ with prob=.07}\\
+\$2M &\quad \text{ with prob=.91}
\end{cases}
$$

### a)

The PMF (probability density function for a discrete random variable) is given by:

$$P(\text{Revenue} = x) = \begin{cases}
.02 &\quad \text{  if $x = -\$4M$}\\
.07 &\quad \text{ if $x = -\$1M$}\\
.91 &\quad \text{ if $x = \$2M$}\\
0 & \quad \text{ otherwise}
\end{cases}
$$

The corresponding CDF of the revenue is given by

$$
F(x)=P(\text{Revenue} \leq x) = \begin{cases}
0 & \quad \text{ if $x < -4$}\\
0.02 & \quad \text{ if $-4 \leq x < -1$}\\
0.09 & \quad \text{ if $-1 \leq x < 2$}\\
1 & \quad \text{ if $2 \leq x $}
\end{cases}
$$

$$P(\text{Revenue} < -\$VaR^\alpha) = \alpha$$

We can provide conservative, formal and interpolated VaR estimates. $\alpha = 5\%$ falls in the bin $[-4,-1)$:

\begin{align}
\text{Conservative \$VaR : }\quad &\$4M,\\
\text{Formal \$VaR : }\quad &\$1M,\\
\text{Interpolated \$VaR : }\quad & \frac{-\$VaR - (-4)}{5-2} = \frac{-1 - (-4)}{9-2} \Rightarrow \$VaR = 4 - 9/7 = \$2.7M,
\end{align}

There is only one possible outcome in which the losses are lower than VaR= \$2.1 M , it is the 2%-chance event when we lose \$4M. This is the expected shortfall.

$$ES = - \frac{\sum_{R_i < -VaR} R_i \cdot P(R_i)}{\sum_{R_i < - VaR} P(R_i)} = - \frac{- 4 \cdot .02}{.02} = \$4M$$


### b)

Distribution of the portfolio with two equally weighted projects is given by the PMF:

>(constructed with the help of a table:)
>|Revenue : proba|-4|-1|2|
>|---|---|---|---|
>|-4|-4 : .02²| -2.5 : .02*.07 |-1 : .02*.91|
>|-1|-2.5 : .07*.02|-1 : .07²|.5 : .07*.91|
>|2|-1: .02*.91|.5 : .07*.91 |2 : .91²|

$$\text{Revenue}_i = 
\begin{cases}
-\$4M &\quad \text{ with prob=} .02^2 = .0004\\
-\$2.5M &\quad \text{with prob=} 2\times .02\times . 07 = .0028\\
-\$1M &\quad \text{ with prob=}.07^2 + 2\times .02\times.91 = .0413\\
+\$.5M &\quad \text{with prob=} 2\times .07\times .91 = .1274\\
+\$2M &\quad \text{ with prob=}.91^2 = .8281
\end{cases}
$$

$$P(\text{Revenue} = x) = \begin{cases}
.0004 &\quad \text{  if $x = -\$4M$}\\
.0028 &\quad \text{ if $x = -\$2.5M$}\\
.0413 &\quad \text{ if $x = -\$1M$}\\
.1274 &\quad \text{ if $x = \$.5M$}\\
.8281 &\quad \text{ if $x = \$2M$}\\
0 & \quad \text{ otherwise}
\end{cases}
$$

The corresponding CDF is

$$
F(x)=P(\text{Revenue} \leq x) = \begin{cases}
0 & \quad \text{ if $x < -4$}\\
.0004& \quad \text{ if $-4 \leq x < -2.5$}\\ 
0.0032 & \quad \text{ if $-2.5 \leq x < -1$}\\
0.0445 & \quad \text{ if $-1 \leq x < .5$}\\
0.1719 & \quad \text{ if $.5 \leq x < 2$}\\
1 & \quad \text{ if $2 \leq x $}
\end{cases}$$

$\alpha = 5\%$ falls into the fourth bin $[-1, .5)$ so that the interpolated VaR is calculated from:

$$\frac{-$VaR - .5}{.05 - .0445} = \frac{2 - .5}{.1274} \Rightarrow VaR = -.5 - 1.5 \times \frac{.0055}{.1274} = .935.$$

For a portfolio it is possible to observe 3 outcomes when losses are higher than VaR = .935, namely 

$$
\begin{cases}
-\$4M &\quad \text{ with prob=}.0004\\
-\$2.5M &\quad \text{with prob=}.0028\\
-\$1M &\quad \text{ with prob=}.0413\\
\end{cases}
$$

The ES is calculated as follows:

$$ES = \frac{4 \times .0004 + 2.5 \times .0028 + 1 \times.0413}{.0004 + .0028 + .0413} \approx \$1.12M$$

### c) 


---

### Quantile regression

Consider the conventional regression:

$$y_t = \alpha + \beta x_t + \epsilon_t,$$

it can be equivalently rewritten down as:

$$E(y) = \alpha + \beta x$$

Expectation of predicted variable $y$ is obtained when we minimize the mean squared error

$$MSE = \frac{1}{n}\sum_{t}(y_t - \alpha - \beta x_t)^2$$

If instead we minimized the mean absolute error:

$$MAE = \frac{1}{n} \sum |y_t - \alpha - \beta x_t|,$$

we construct a model for the *median* of $y$, (not the mean). Median is the 50% quantile.

Both mean and median are statistics of the entire distribution that do not describe it fully. Just as we model these statistics we can go further and try to model others, e.g. 5% quantile ($-VaR^{5\%}$).
