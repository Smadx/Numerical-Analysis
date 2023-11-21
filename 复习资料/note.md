# 计算方法复习笔记
## 0.矩阵和向量范数
### 0.1向量范数
向量$X(x_i)$的$L_p$范数定义为： 
$$||X||_p=(\sum_{i=1}^n|x_i|^p)^{ \frac{1}{p} }$$
其中$p\geq 1$，当$p=2$时，称为欧几里得范数，即： 
$$||X||_2=(\sum_{i=1}^n|x_i|^2)^{ \frac{1}{2} }$$
当$p=1$时，称为$L_1$范数，即： 
$$||X||_1=\sum_{i=1}^n|x_i|$$
当$p=\infty$时，称为$L_\infty$范数，即： 
$$||X||_\infty=max_{1\leq i\leq n}|x_i|$$
### 0.2矩阵范数
矩阵$A(a_{ij})$的$L_p$范数定义为： 
$$||A||_p=(\sum_{i=1}^m\sum_{j=1}^n|a_{ij}|^p)^{ \frac{1}{p} }$$
其中$p\geq 1$，当$p=2$时，称为Frobenius范数，即： 
$$||A||_F=(\sum_{i=1}^m\sum_{j=1}^n|a_{ij}|^2)^{ \frac{1}{2} }$$
当$p=1$时，称为$L_1$范数，即： 
$$||A||_1=max_{1\leq j\leq n}\sum_{i=1}^m|a_{ij}|$$
当$p=\infty$时，称为$L_\infty$范数，即： 
$$||A||_\infty=max_{1\leq i\leq m}\sum_{j=1}^n|a_{ij}|$$
当$p=2$时，称为谱范数，即： 
$$||A||_2=max_{1\leq i\leq m}|\lambda_i|$$
其中$\lambda_i$为$A^TA$的特征值。 
### 0.3矩阵的条件数和谱半径
矩阵$A$的条件数定义为：
$$cond(A)=||A||\cdot||A^{-1}||$$
其中$||\cdot||$为矩阵范数，$A^{-1}$为$A$的逆矩阵。 
矩阵$A$的谱半径定义为：
$$\rho(A)=max_{1\leq i\leq n}|\lambda_i|$$
其中$\lambda_i$为$A$的特征值。 
矩阵的条件数越大，矩阵越病态，矩阵的谱半径越大，矩阵越不稳定。 
## 1.插值
### 1.1Langrange插值
给定$n+1$个点$(x_i,y_i)$，其中$x_i$两两不同，$y_i$为函数$f(x)$在$x_i$处的函数值，$i=0,1,2,...,n$，则称函数$P_n(x)$为函数$f(x)$在区间$[x_0,x_n]$上的Langrange插值多项式，其表达式为：
$$P_n(x)=\sum_{i=0}^ny_i\prod_{j=0,j\neq i}^n\frac{x-x_j}{x_i-x_j}$$
误差为：
$$f(x)-P_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\prod_{i=0}^n(x-x_i)$$
其中$\xi\in[x_0,x_n]$。
### 1.2Newton插值
给定$n+1$个点$(x_i,y_i)$，其中$x_i$两两不同，$y_i$为函数$f(x)$在$x_i$处的函数值，$i=0,1,2,...,n$，则称函数$P_n(x)$为函数$f(x)$在区间$[x_0,x_n]$上的Newton插值多项式，其表达式为：
$$P_n(x)=\sum_{i=0}^na_i\prod_{j=0}^{i-1}(x-x_j)$$
其中$a_i$为差商，其表达式为：
$$a_i=f[x_0,x_1,...,x_i]=\frac{f[x_1,x_2,...,x_i]-f[x_0,x_1,...,x_{i-1}]}{x_i-x_0}$$
误差为：
$$f(x)-P_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\prod_{i=0}^n(x-x_i)$$
其中$\xi\in[x_0,x_n]$。
### 1.3Hermite插值
给定$2n+2$个点$(x_i,y_i)$，其中$x_i$两两不同，$y_i$为函数$f(x)$在$x_i$处的函数值，$i=0,1,2,...,n$，$y_{n+i}$为函数$f(x)$在$x_i$处的导数值，$i=0,1,2,...,n$，则称函数$P_{2n+1}(x)$为函数$f(x)$在区间$[x_0,x_n]$上的Hermite插值多项式，其表达式为：
$$P_{2n+1}(x)=\sum_{i=0}^n\left(y_i\prod_{j=0,j\neq i}^n\frac{x-x_j}{x_i-x_j}+y_{n+i}\prod_{j=0,j\neq i}^n\frac{x-x_j}{x_i-x_j}\int_{x_0}^x\prod_{k=0,k\neq i}^n\frac{t-x_k}{x_i-x_k}dt\right)$$
误差为：
$$f(x)-P_{2n+1}(x)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\prod_{i=0}^n(x-x_i)^2$$
其中$\xi\in[x_0,x_n]$。
### 1.4三次样条插值
给定$n+1$个点$(x_i,y_i)$，其中$x_i$两两不同，$y_i$为函数$f(x)$在$x_i$处的函数值，$i=0,1,2,...,n$，则称函数$S(x)$为函数$f(x)$在区间$[x_0,x_n]$上的三次样条插值函数，其表达式为：
$$S(x)=\left \{ \begin{array}{ll}
S_i(x)=a_i+b_i(x-x_i)+c_i(x-x_i)^2+d_i(x-x_i)^3,&x\in[x_i,x_{i+1}]\\
S(x_i)=y_i,&i=0,1,2,...,n
\end{array}\right.$$
其中$a_i=y_i$，$b_i=\frac{y_{i+1}-y_i}{h_i}-\frac{h_i}{3}(2c_i+c_{i+1})$，$c_i=\frac{1}{h_i}(m_i-m_{i-1})$，$d_i=\frac{1}{3h_i}(m_{i-1}+2m_i)$，$h_i=x_{i+1}-x_i$，$m_i=S''(x_i)$，$i=0,1,2,...,n$，$S''(x)$为$S(x)$的二阶导数。
其中$m_i$满足以下方程组：
$$\left \{ \begin{array}{ll}
\frac{h_i}{6}m_{i-1}+\frac{h_i+h_{i+1}}{3}m_i+\frac{h_{i+1}}{6}m_{i+1}=f[x_i,x_{i+1}],&i=1,2,...,n-1\\
m_0=0\\
m_n=0
\end{array}\right.$$
## 2.最小二乘拟合
### 2.1最小二乘法
给定$n+1$个点$(x_i,y_i)$，其中$x_i$两两不同，$y_i$为函数$f(x)$在$x_i$处的函数值，$i=0,1,2,...,n$，则称函数$P_n(x)$为函数$f(x)$在区间$[x_0,x_n]$上的最小二乘拟合多项式，其表达式为：
$$P_n(x)=\sum_{i=0}^na_ix^i$$
其中$a_i$满足以下方程组：
$$\left \{ \begin{array}{ll}
\sum_{i=0}^na_ix_i^j=\sum_{i=0}^ny_ix_i^j,&j=0,1,2,...,n\\
\end{array}\right.$$
### 2.3矛盾方程组
当线性方程组的系数矩阵不是满秩矩阵时，称其为矛盾方程组。使$||A\alpha-Y||_2$最小的$\alpha$称为矛盾方程组的最小二乘解。 
给定矩阵$A\in R^{m\times n}$，向量$Y\in R^m$，有:
$rank A^TA=rank A^T$ 
线性方程组$A^TA\alpha=A^TY$恒有解$\alpha$ 
$||A\alpha -Y||_2$达最小值,当且仅当$A^TA\alpha=A^TY$ 
使$||A\alpha-Y||_2$最小的$\alpha$是唯一的当且仅当$rank A=n$ 
## 3.非线性方程求解
### 3.1对分法
给定函数$f(x)$，若在区间$[a,b]$上$f(a)f(b)<0$，则在区间$[a,b]$上至少存在一个根$x^*$，迭代格式：
$$x_{n+1}=\frac{a_n+b_n}{2}$$
且有：
$$|x^*-x_n|\leq\frac{b-a}{2^n}$$
其中$x_n$为区间$[a,b]$上的第$n$次迭代值。
### 3.2不动点迭代
对给定的方程$f(x)=0$，将它转化为等价形式$x=\phi (x)$，其中$\phi (x)$为迭代函数，若在区间$[a,b]$上$\phi (x)\in [a,b]$，且$\phi '(x)$在$[a,b]$上存在且连续，且$|\phi '(x)|\le L<1$，则在区间$[a,b]$上至少存在一个根$x^*$，迭代格式:
$$x_{n+1}=\phi (x_n)$$
且有：
$$|x^*-x_n|\leq\frac{L^n}{1-L}|x_1-x_0|$$
其中$x_n$为区间$[a,b]$上的第$n$次迭代值。
### 3.3牛顿迭代
对给定的方程$f(x)=0$，若在区间$[a,b]$上$f'(x)$在$[a,b]$上存在且连续，且$f'(x)\neq 0$，则在区间$[a,b]$上至少存在一个根$x^*$，迭代格式:
$$x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)}$$
且有：
$$|x^*-x_n|\leq\frac{f(x_n)}{f'(x_n)}$$
其中$x_n$为区间$[a,b]$上的第$n$次迭代值。
收敛阶:
$$\lim_{n\rightarrow\infty}\frac{|x_{n+1}-x^*|}{|x_n-x^*|^2}=\frac{1}{2}$$
### 3.4弦截法
对给定的方程$f(x)=0$，若在区间$[a,b]$上$f'(x)$在$[a,b]$上存在且连续，且$f'(x)\neq 0$，则在区间$[a,b]$上至少存在一个根$x^*$，迭代格式:
$$x_{n+1}=x_n-\frac{f(x_n)(x_n-x_{n-1})}{f(x_n)-f(x_{n-1})}$$
且有：
$$|x^*-x_n|\leq\frac{f(x_n)(x_n-x_{n-1})}{f(x_n)-f(x_{n-1})}$$
### 3.5收敛阶分析
设$\alpha$为方程$f(x)=0$的根，若存在常数$C>0$，使得：
$$|x_{n+1}-\alpha|\leq C|x_n-\alpha|^p$$
则称迭代格式具有$p$阶收敛性，$p$为收敛阶。
### 3.6求解非线性方程组的牛顿迭代法
对给定的方程组$f(x)=0$，若在区间$[a,b]$上$f'(x)$在$[a,b]$上存在且连续，且$f'(x)\neq 0$，则在区间$[a,b]$上至少存在一个根$x^*$，迭代格式:
$$x_{n+1}=x_n-J^{-1}(x_n)f(x_n)$$
其中$J(x_n)$为$f(x)$在$x_n$处的雅可比矩阵。
## 4.求解线性方程组的直接法
### 4.1.1高斯顺序消元法
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，高斯顺序消元法的步骤如下：
1. 化方程组为上三角形方程组
2. 回代求解
### 4.1.2高斯列主元消元法
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，高斯列主元消元法的步骤如下：
1. 选出第$n$列中绝对值最大的元素所在的行，将该行与第$n$行交换
2. 化方程组为上三角形方程组
3. 回代求解
### 4.2.1Doolittle分解
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，Doolittle分解的步骤如下：
1. $L_{11}=1$，$U_{11}=a_{11}$
2. $U_{1j}=a_{1j}$，$L_{j1}=\frac{a_{j1}}{U_{11}}$
3. $U_{ij}=a_{ij}-\sum_{k=1}^{i-1}L_{ik}U_{kj}$，$L_{ji}=\frac{1}{U_{ii}}(a_{ji}-\sum_{k=1}^{i-1}L_{jk}U_{ki})$
### 4.2.2Crout分解
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，Crout分解的步骤如下：
1. $U_{11}=1$，$L_{11}=a_{11}$
2. $L_{i1}=a_{i1}$，$U_{1i}=\frac{a_{1i}}{L_{11}}$
3. $L_{ij}=a_{ij}-\sum_{k=1}^{j-1}L_{ik}U_{kj}$，$U_{ji}=\frac{1}{L_{jj}}(a_{ji}-\sum_{k=1}^{j-1}L_{jk}U_{ki})$
### 4.3.1追赶法求解三对角线性方程组
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，追赶法求解三对角线性方程组的步骤如下：
1. $c_1=\frac{a_{12}}{a_{11}}$，$d_1=\frac{b_1}{a_{11}}$
2. $c_i=\frac{a_{i,i+1}}{a_{ii}-a_{i-1,i}c_{i-1}}$，$d_i=\frac{b_i-a_{i-1,i}d_{i-1}}{a_{ii}-a_{i-1,i}c_{i-1}}$
3. $x_n=d_n$，$x_i=d_i-c_ix_{i+1}$
### 4.3.2$LDL^T$分解
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，$LDL^T$分解的步骤如下：
1. $L_{11}=1$，$D_{11}=a_{11}$
2. $L_{i1}=\frac{a_{i1}}{D_{11}}$，$D_{ii}=a_{ii}-\sum_{k=1}^{i-1}L_{ik}^2D_{kk}$，$L_{ji}=\frac{1}{D_{jj}}(a_{ji}-\sum_{k=1}^{j-1}L_{jk}L_{ik}D_{kk})$
## 5.求解线性方程组的迭代法
### 5.1Jacobi迭代法
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，Jacobi迭代法的步骤如下：
1. $x^{(0)}$为任意给定的初始向量
2. $x^{(k+1)}=D^{-1}(b-(A-D)x^{(k)})$
其中$D$为$A$的对角线元素组成的对角矩阵。 
收敛条件：$A$为严格对角占优矩阵
### 5.2Gauss-Seidel迭代法
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，Gauss-Seidel迭代法的步骤如下：
1. $x^{(0)}$为任意给定的初始向量
2. $x^{(k+1)}=-(D+L)^{-1}Ux^{(k)}+(D+L)^{-1}b$
其中$D$为$A$的对角线元素组成的对角矩阵，$L$为$A$的下三角矩阵，$U$为$A$的上三角矩阵。 
收敛条件：$A$为严格对角占优矩阵 
### 5.3超松弛迭代法
对于线性方程组$Ax=b$，其中$A\in R^{n\times n}$，$b\in R^n$，$x\in R^n$，若$A$非奇异，则存在唯一解$x=A^{-1}b$，超松弛迭代法的步骤如下：
1. $x^{(0)}$为任意给定的初始向量
2. $x^{(k+1)}=(I+\omega D^{-1}L)^{-1}[(1-\omega)I-\omega D^{-1}U]x^{(k)}+\omega(I+\omega D^{-1}L)^{-1}D^{-1}b$
其中$D$为$A$的对角线元素组成的对角矩阵，$L$为$A$的下三角矩阵，$U$为$A$的上三角矩阵。 
收敛条件：$A$为严格对角占优矩阵，且$0<\omega<2$
## 6.数值微分和数值积分
### 6.1数值微分
#### 6.1.1前向差分
$$f'(x_0)=\frac{f(x_0+h)-f(x_0)}{h}+O(h)$$
#### 6.1.2后向差分
$$f'(x_0)=\frac{f(x_0)-f(x_0-h)}{h}+O(h)$$
#### 6.1.3中心差分
$$f'(x_0)=\frac{f(x_0+h)-f(x_0-h)}{2h}+O(h^2)$$
#### 6.1.4插值型数值微分
$\cdots$
### 6.2Newton-Cotes数值积分
#### 6.2.1代数精度
记$[a,b]$上以$\{x_i,i=0,1,\cdots,n\}$为积分节点的数值积分公式为:
$$I_n(f)=\sum_{i=0}^n\alpha_if(x_i)$$
若$I_n(f)$满足:
$$E_n(x^k)=I(x^k)-I_n(x^k)=0,k=0,1,\cdots,m$$
而
$$E_n(x^{m+1})\neq0$$
则称$I_n(f)$具有$m$阶代数精度。
#### 6.2.1插值型数值积分
$\cdots$
#### 6.2.2梯形公式
$$I_1(f)=\frac{b-a}{2}[f(a)+f(b)]$$
截断误差:
$$E_1(f)=-\frac{1}{12}f''(\xi)(b-a)^3$$
#### 6.2.3Simpson公式
$$I_2(f)=\frac{b-a}{6}[f(a)+4f(\frac{a+b}{2})+f(b)]$$
截断误差:
$$E_2(f)=-\frac{1}{90}f^{(4)}(\xi)(b-a)^5$$
### 6.3复化Newton-Cotes公式
#### 6.3.1复化梯形公式
$$I_n(f)=\frac{h}{2}[f(x_0)+2\sum_{i=1}^{n-1}f(x_i)+f(x_n)]$$
截断误差:
$$E_n(f)=-\frac{b-a}{12}h^2f''(\xi)$$
#### 6.3.2复化Simpson公式
$$I_n(f)=\frac{h}{3}[f(x_0)+2\sum_{i=1}^{n-1}f(x_{2i})+4\sum_{i=1}^{n}f(x_{2i-1})+f(x_n)]$$
截断误差:
$$E_n(f)=-\frac{b-a}{180}h^4f^{(4)}(\xi)$$
## 7.常微分方程数值解
### 7.1常微分方程初值问题
#### 7.1.1向前Euler法
$$y_{n+1}=y_n+hf(x_n,y_n)$$
#### 7.1.2向后Euler法
$$y_{n+1}=y_n+hf(x_{n+1},y_{n+1})$$
#### 7.1.3局部截断误差与整体截断误差
局部截断:$$T_(n+1)=O(h^{p+1})$$
整体截断:$$e^{L(b-a)}\frac{T}{Lh}$$
其中$T=O(h^2)$为局部截断误差，$L$为Lipschitz常数
### 7.2线性多步法
用Runge-Kutta法起步，然后用前面的计算结果递推后面的计算结果
#### 7.2.1Adams-Bashforth公式
$$y_{n+1}=y_n+\frac{h}{2}[3f(x_n,y_n)-f(x_{n-1},y_{n-1})]$$
#### 7.2.2Adams-Moulton公式
$$y_{n+1}=y_n+\frac{h}{12}[5f(x_{n+1},y_{n+1})+8f(x_n,y_n)-f(x_{n-1},y_{n-1})]$$
$\cdots$
## 8.计算矩阵的特征值和特征向量
### 8.1幂法
#### 8.1.1幂法
$$x^{(k+1)}=\frac{Ax^{ (k) } }{||Ax^{(k)}||_\infty}$$
按模最大
#### 8.1.2反幂法
$$x^{(k+1)}=\frac{A^{-1}x^{ (k) } }{||A^{-1}x^{(k)}||_\infty}$$
按模最小
### 8.2实对称矩阵的Jacobi方法
#### 8.2.1Jacobi方法
把$A$分解为
$$A=QDQ^T$$
其中$Q$为正交矩阵，$D$为对角矩阵 
利用$Givens$旋转变换:
$$G(i,j,\theta)=\begin{bmatrix}1&\cdots&0&\cdots&0&\cdots&0\\
\vdots&\ddots&\vdots&&\vdots&&\vdots\\
0&\cdots&c&\cdots&s&\cdots&0\\
\vdots&&\vdots&\ddots&\vdots&&\vdots\\
0&\cdots&-s&\cdots&c&\cdots&0\\
\vdots&&\vdots&&\vdots&\ddots&\vdots\\
0&\cdots&0&\cdots&0&\cdots&1\end{bmatrix}$$
其中$c=\cos\theta,s=\sin\theta$，$G(i,j,\theta)$为正交矩阵，且
$$G(i,j,\theta)A=\begin{bmatrix}a_{11}&\cdots&a_{1i}&\cdots&a_{1j}&\cdots&a_{1n}\\
\vdots&\ddots&\vdots&&\vdots&&\vdots\\
a_{i1}&\cdots&a_{ii}&\cdots&a_{ij}&\cdots&a_{in}\\
\vdots&&\vdots&\ddots&\vdots&&\vdots\\
a_{j1}&\cdots&a_{ji}&\cdots&a_{jj}&\cdots&a_{jn}\\
\vdots&&\vdots&&\vdots&\ddots&\vdots\\
a_{n1}&\cdots&a_{ni}&\cdots&a_{nj}&\cdots&a_{nn}\end{bmatrix}$$
其中
$$\begin{cases}a_{ii}=c^2a_{ii}-2csa_{ij}+s^2a_{jj}\\
a_{jj}=s^2a_{ii}+2csa_{ij}+c^2a_{jj}\\
a_{ij}=a_{ji}=(c^2-s^2)a_{ij}+sc(a_{ii}-a_{jj})\end{cases}$$
故Jacobi方法为：
1. 选取非对角线按模最大元素$|a_{pq}|=\max_{i\neq j} |a_{ij}|$
2. 确定旋转角度$\theta$，其中$t=tan \theta$,令$s=\frac{a_{pq}-a_{pp} }{2a_{pq} },$
3. 若s=0，则$t=1$，否则$t=\max \{-s-\sqrt{s^2+1},-s+\sqrt{s^2+1} \}$
4. $cos\theta =\frac{1}{\sqrt{1+t^2}}$, $sin\theta =t\cos\theta$
5. 计算$Q(i,j,\theta)$
6. $A=Q^TAQ$
### 8.3QR方法
#### 8.3.1QR方法
1. $A_0=A$
2. $A_k=Q_kR_k$
3. $A_{k+1}=R_kQ_k$
4. $A_k$收敛到上三角矩阵
其中$Q_k$为正交矩阵，$R_k$为上三角矩阵
#### 8.3.2矩阵的QR分解
对于矩阵$A=[\alpha_1,\alpha_2,\cdots,\alpha_n]$,利用Householder变换
$$H=I-2uu^T$$
其中$u=\frac{\alpha_1-sgn(\alpha_1)||\alpha_1||_2e_1}{||\alpha_1-sgn(\alpha_1)||\alpha_1||_2e_1||_2}$ 
于是
$$H_1A=\begin{bmatrix}
*&*\\
0&\hat{A}\\
\end{bmatrix}$$
继续对$\hat{A}$进行Householder变换.. 
最终$$R=H_nH_{n-1}\cdots H_1A$$
$$Q^T=H_nH_{n-1}\cdots H_1$$
## 9.函数逼近
### 9.1赋范线性空间
#### 9.1.1范数
三要素:
1. 非负性: $||x||\geq 0$, $||x||=0$当且仅当$x=0$
2. 齐次性: $||\alpha x||=|\alpha|||x||$
3. 三角不等式: $||x+y||\leq ||x||+||y||$
#### 9.1.2距离
$$d(x,y)=||x-y||$$
#### 9.1.3内积
1. 正定性: $(x,x)\geq 0$, $(x,x)=0$当且仅当$x=0$
2. 线性性: $(\alpha x+y,z)=\alpha(x,z)+(y,z)$
3. 共轭对称性: $(x,y)=(y,x)^*$
4. 内积的范数: $||x||=(x,x)^{\frac{1}{2}}$
#### 9.1.4最佳逼近
如果有元素$m^*\in M$使得
$$||x-m^*||=inf_{m\in M}||x-m||=^{\Delta}d(x,M)$$
则称$m^*$为$x$在$M$中的最佳逼近元素 
#### 9.1.5Schmidt正交化
1. $x_1=\alpha_1$
2. $x_2=\alpha_2-\frac{(x_1,\alpha_2)}{(x_1,x_1)}x_1$
3. $x_k=\alpha_k-\sum_{i=1}^{k-1}\frac{(x_i,\alpha_k)}{(x_i,x_i)}x_i$
### 9.2最佳平方逼近
#### 9.2.1最佳平方逼近
如果有元素$m^*\in M$使得
$$||Ax-m^*||=inf_{m\in M}||Ax-m||=^{\Delta}d(Ax,M)$$
则称$m^*$为$Ax$在$M$中的最佳平方逼近元素，其中$||\cdot||$为内积范数，且内积为
$L_\rho ^2[a,b]$上的内积:
$$(f,g)=\int_a^b f(x)g(x)\rho(x)dx$$
### 9.3正交多项式
正交多项式基$\{g_l^*(x)\}_{l=0}^n$有递推关系:
$$g_0^*(0)=1,g_1^*(x)=x-\frac{(xg_0^*,g_0^*)}{(g_0^*,g_0^*)}$$
$$g_{l+1}^*(x)=(x-\frac{(xg_l^*,g_l^*)}{(g_l^*,g_l^*)})g_l^*(x)-\frac{(xg_l^*,g_{l-1}^*)}{(g_{l-1}^*,g_{l-1}^*)}g_{l-1}^*(x)$$
### 9.4离散傅里叶变换
#### 9.4.1离散傅里叶变换
$$\hat{f}_k=\sum_{j=0}^{n-1}f_je^{-\frac{2\pi i}{n}jk},k=0,1,\cdots,n-1$$
$$f_j=\frac{1}{n}\sum_{k=0}^{n-1}\hat{f}_ke^{\frac{2\pi i}{n}jk},j=0,1,\cdots,n-1$$
### 9.5最佳一致逼近多项式
#### 9.5.1最佳一致逼近多项式
如果有元素$p^*\in P_n$使得
$$||f-p^*||_{\infty}=inf_{p\in P_n}||f-p||_{\infty}=^{\Delta}d(f,P_n)$$
则称$p^*$为$f$在$P_n$中的最佳一致逼近多项式
#### 9.5.2最佳一致逼近多项式的存在唯一性
1. 存在性: $P_n$是有限维线性空间，$||\cdot||_{\infty}$是范数，所以$P_n$是赋范线性空间，由于$P_n$是有限维线性空间，所以$P_n$是完备的，所以由最佳逼近定理可知，$f$在$P_n$中存在最佳逼近多项式
2. 唯一性: 假设$f$在$P_n$中存在两个最佳逼近多项式$p_1^*,p_2^*$,则有
$$||f-p_1^*||_{\infty}=||f-p_2^*||_{\infty}=d(f,P_n)$$
于是
$$||p_1^*-p_2^*||_{\infty}=||p_1^*-f+f-p_2^*||_{\infty}\leq ||p_1^*-f||_{\infty}+||f-p_2^*||_{\infty}=0$$
所以$p_1^*=p_2^*$
#### 9.5.3交错点组
设函数$f\in C[a,b]$,则$p^*$是函数$f$在$P_n$中的最佳一致逼近多项式的充要条件是$f-p^*$在$[a,b]$上有$n+2$个交错点:$a\le x_0\le x_1\le \cdots \le x_{n+1}\le b$使得
$$f(x_i)-p^*(x_i)=(-1)^i\sigma ||f-p^*||_{\infty},i=0,1,\cdots,n+1$$
又有:若$p^*$是$f$的n次最佳一致逼近多项式,如果$f$在$[a,b]$上有$n+1$阶导数,且$f^(n+1)(x)$在$[a,b]$上不变号,则$f-p^*$在$[a,b]$上有$n+2$个交错点,包括端点。
## 10.最优化方法
### 10.1线性规划
#### 10.1.1线性规划的标准形式
$$max\quad c^Tx$$
$$s.t.\quad Ax=b$$
$$\quad\quad x_i\geq 0$$
其中$A\in R^{m\times n},b\in R^m,c\in R^n$
基解:
$$x=(x_{B_1},x_{B_2},\cdots,x_{B_m},0,0,\cdots,0)^T$$
$$x_{B_i}=\sum_{j=1}^m\alpha_{ij}b_j$$
$$\alpha_{ij}=\left \{\begin{matrix}
1,i=j\\
0,i\neq j
\end{matrix}\right.$$
基可行解:
$$x=(x_{B_1},x_{B_2},\cdots,x_{B_m},0,0,\cdots,0)^T$$
$$x_{B_i}=\sum_{j=1}^m\alpha_{ij}b_j$$
$$\alpha_{ij}\geq 0,i,j=1,2,\cdots,m$$
#### 10.1.2单纯形法
单纯性表:
$$\begin{bmatrix}
x_{B_1}&x_{B_2}&\cdots&x_{B_m}&x_{N_1}&x_{N_2}&\cdots&x_{N_{n-m}}&b\\
1&0&\cdots&0&-c_{N_1}&-c_{N_2}&\cdots&-c_{N_{n-m}}&0\\
0&1&\cdots&0&-c_{N_1}&-c_{N_2}&\cdots&-c_{N_{n-m}}&0\\
\vdots&\vdots&\ddots&\vdots&\vdots&\vdots&\ddots&\vdots&\vdots\\
0&0&\cdots&1&-c_{N_1}&-c_{N_2}&\cdots&-c_{N_{n-m}}&0\\
0&0&\cdots&0&-c_{N_1}&-c_{N_2}&\cdots&-c_{N_{n-m}}&0
\end{bmatrix}$$
### 10.2非线性规划
#### 10.2.1一维搜索
首先用进退法确定搜索区间$[a,b]$,然后用黄金分割法或斐波那契法确定搜索点
进退法:
$$\alpha_0=0,\alpha_1=\alpha_0+\Delta$$
$$\alpha_{k+1}=\alpha_k+\Delta$$
$$\alpha_{k+n}=\alpha_k+2^{n-1}\Delta$$
1. 黄金分割法:
$$\alpha_k=\frac{1}{2}(\frac{\sqrt{5}-1}{2})^k$$
1. 斐波那契法:
$$\alpha_k=\frac{F_{n-k} }{F_{n-k+1} }$$
#### 10.2.2梯度下降法
$$x^{(k+1)}=x^{(k)}-\alpha_k\nabla f(x^{(k)})$$
#### 10.2.3牛顿法
$$x^{(k+1)}=x^{(k)}-[\nabla^2f(x^{(k)})]^{-1}\nabla f(x^{(k)})$$
其中$\nabla^2f(x^{(k)})$为$f(x)$的海森矩阵
#### 10.2.4拟牛顿法
$$x^{(k+1)}=x^{(k)}-B_k^{-1}\nabla f(x^{(k)})$$
其中$B_k$为$f(x)$的拟海森矩阵
#### 10.2.5共轭梯度法
$$x^{(k+1)}=x^{(k)}+\alpha_kd^{(k)}$$
$$d^{(k+1)}=-\nabla f(x^{(k+1)})+\beta_kd^{(k)}$$
其中$\alpha_k,\beta_k$为确定的步长
