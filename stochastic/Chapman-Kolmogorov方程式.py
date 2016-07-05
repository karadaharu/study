
# coding: utf-8

# # Table of Contents
#  <p><div class="lev2"><a href="#Chapman-Kolmogorov方程式"><span class="toc-item-num">0.1&nbsp;&nbsp;</span>Chapman-Kolmogorov方程式</a></div><div class="lev2"><a href="#Chapman-Kolmogorov方程式の微分形"><span class="toc-item-num">0.2&nbsp;&nbsp;</span>Chapman-Kolmogorov方程式の微分形</a></div><div class="lev3"><a href="#準備"><span class="toc-item-num">0.2.1&nbsp;&nbsp;</span>準備</a></div><div class="lev3"><a href="#導出"><span class="toc-item-num">0.2.2&nbsp;&nbsp;</span>導出</a></div><div class="lev3"><a href="#解釈"><span class="toc-item-num">0.2.3&nbsp;&nbsp;</span>解釈</a></div><div class="lev4"><a href="#Jump-Processs:-マスター方程式"><span class="toc-item-num">0.2.3.1&nbsp;&nbsp;</span>Jump Processs: マスター方程式</a></div><div class="lev4"><a href="#Diffusion-Process:-Fokker-Planck方程式"><span class="toc-item-num">0.2.3.2&nbsp;&nbsp;</span>Diffusion Process: Fokker-Planck方程式</a></div><div class="lev4"><a href="#決定論的過程-Liouville方程式"><span class="toc-item-num">0.2.3.3&nbsp;&nbsp;</span>決定論的過程 Liouville方程式</a></div>

# ## Chapman-Kolmogorov方程式
# 確率過程に対して、マルコフ性が成り立つときに、確率過程の途中の経路をすべて積分すれば、任意の時刻についての、条件付き確率が得られる
# 
# $p(x_{1},t_{1}|x_{3},t_{3})=\int p(x_{1}, t_{1}| x_{2}, t_{2}) p(x_{2}, t_{2}| x_{3}, t_{3}) dx_{2}$

# ## Chapman-Kolmogorov方程式の微分形
# ### 準備
# 適切な仮定のもとでは、CK方程式は微分形に変形できる。以下の仮定をおく。
# 
# * i) $\lim_{\Delta t \rightarrow 0} p(x, t+\Delta t| z, t) / \Delta t = W(x | z, t)$  uniformly in x,z and t for $|x-z|\geq \varepsilon$
# * ii) $\lim_{\Delta t \rightarrow 0} \frac{1}{\Delta t} \int_{|x-z| < \varepsilon} dx(x_{i} - z_{i}) p(x, t+\Delta t| z, t) = A_{i}(z, t) + O(\varepsilon)$ : 微小時間後の座標についての期待値(距離εの中だけ)は、zとtの関数+誤差で表される？
# * ii) $\lim_{\Delta t \rightarrow 0} \frac{1}{\Delta t} \int_{|x-z| < \varepsilon} dx(x_{i} - z_{i}) (x_{j} - z_{j} )p(x, t+\Delta t| z, t) = B_{i, j}(z, t) + O(\varepsilon) $
# 
# ii)、ii)はz, ε, tについてuniform。
# 
# ここで、ii) ii) の形の3次以上の係数は消えることに注意。例えば、3次の項は
# 
# $\lim_{\Delta t \rightarrow 0} \frac{1}{\Delta t} \int_{|x-z| < \varepsilon} dx(x_{i} - z_{i}) (x_{j} - z_{j} )  (x_{k} - z_{k} ) p(x, t+\Delta t| z, t) = C_{i, j,k}(z, t) + O(\varepsilon) $
# 
# $C_{i,j,k}$は対称なので、
# 
# $\sum_{i,j,k}\alpha_{i}\alpha_{j}\alpha_{k}C_{i,j,k}(z,t)\equiv {\bar C}(\alpha, z,t)$
# 
# とおくと
# 
# $C_{i,j,k} = \frac{1}{3!}\frac{\partial^{3}}{\partial \alpha_{i} \alpha_{j} \alpha_{k}}{\bar C}(\alpha, z, t)$
# 
# です。すると、
# 
# $|{\bar C}(\alpha, z, t)| \leq \lim_{\Delta t \rightarrow 0} \frac{1}{\Delta t} \int_{|x-z| < \varepsilon} dx | \alpha (x - z)| [\alpha(x - z )]^{2} p(x, t+\Delta t| z, t) + O(\varepsilon) \\
# = \cdots \\
# = O(\varepsilon)$
# 
# となり、Cは0です。
# 
# もし$W(x| z, t)$がすべての$x \neq z$で消えるとき、その過程は連続な経路のみをもつことができます。(微小時間の極限をとったときに、違う値である確率は0→連続である）
# 
# 
# ### 導出
# 関数f(z)の期待値の時間発展を考えます。
# 
# $\partial_{t} \int dx f(x)p(x, t| y,t')\\
# =\lim_{\Delta \rightarrow 0} \{\int dx f(x) [p(x, t+\Delta t| y,t') - p(x,t| y, t') ]\} / \Delta t\\
# =\lim_{\Delta \rightarrow 0} \{\int dx \int dz f(x) [p(x, t+\Delta t| z,t)p(z,t | y,t') - \int dz f(z) p(z, t| y, t') ]\} / \Delta t\\$ (1)

# xについての積分を$|x-z| \geq \varepsilon$と$|x-z| < \varepsilon$のふたつの領域に分ける。$|x-z| < \varepsilon$のとき、f(z)は2階微分可能な連続関数であるという仮定をおいたので、
# 
# $f(x)=f(z)+\sum_{i} \frac{\partial f(z)}{\partial z_{i}}(x_{i}-z_{i})\\
# + \sum_{i,j} \frac{1}{2} \frac{\partial^{2}f(z)}{\partial z_{i} \partial z_{j} } (x_{i} - z_{i}) (x_{j} - z_{j}) + |x-z|^{2} R(x,z) $
# 
# と展開できる。式(1)に代入すると、
# 
# $(1) = \lim_{\Delta t \rightarrow 0} \frac{1}{\Delta t} \{ 
# \int \int_{|x-z| < \varepsilon} dxdz [\sum_{i} \frac{\partial f(z)}{\partial z_{i}}(x_{i}-z_{i}) + \sum_{i,j} \frac{1}{2} \frac{\partial^{2}f(z)}{\partial z_{i} \partial z_{j} } (x_{i} - z_{i}) (x_{j} - z_{j}) ] p(x, t+\Delta t | z, t) p(z, t| y, t') + \\
# +\int \int_{|x-z| < \varepsilon} dxdz  |x-z|^{2} R(x,z)p(x, t+\Delta t | z, t) p(z, t| y, t') \\
# + \int \int_{|x-z| \geq \varepsilon} dxdz f(x) p(x, t+\Delta t | z, t) p(z, t| y, t') \\
# +\int \int_{|x-z| < \varepsilon} dxdz f(z) p(x, t+\Delta t | z, t) p(z, t| y, t') \\
# -\int \int dxdz f(z) p(z, t+\Delta t| z, t) p(z,t| y, t') \}$

# この右辺は一行ずつ解釈することができます。
# 
# **1行目**
# 単一の収束を仮定すると、積分のなかの極限をとって、以下を得る(条件ii, 条件iiiを使って)
# 
# $\int dz[ \sum_{i} A_{i}(z)\frac{\partial f}{\partial z_{i}} + \frac{1}{2}\sum_{i,j} B_{i,j}(z) \frac{\partial^{2} f}{\partial z_{i} \partial z_{j}}]p(z,t|y,t')+O(\varepsilon)$
# 
# **2行目**  
# 残差項で$\varepsilon \rightarrow 0$で消えます。
# 
# $|\int \int_{|x-z| < \varepsilon} dx  |x-z|^{2} R(x,z)p(x, t+\Delta t | z, t) |  \leq [\int \int_{|x-z| < \varepsilon} dx  |x-z|^{2} p(x, t+\Delta t | z, t) ] {\rm Max}_{|x-z|< \varepsilon} |R(x,z)| \\
# \rightarrow [\sum_{i}B_{i,i}(z,t) + O(\varepsilon)] \{{\rm Max}_{|x-z|< \varepsilon} |R(x,z)|\}$
# 
# となり、波括弧の項は$\varepsilon \rightarrow 0$で消えるからです。
# 
# **3~5行目**  
# すべてをまとめると、
# 
# $\int\int_{|x-z|\geq \varepsilon} dxdz f(z) [W(z|x,t)p(x,t|y,t')-W(x|z,t)p(z,t|y,t')]$
# 
# $\varepsilon \rightarrow 0$の極限をとると、
# 
# $\partial_{t} \int dx f(x)p(x, t| y,t') = \int dz [\sum_{i} A_{i}(z)\frac{\partial f}{\partial z_{i}} + \frac{1}{2}\sum_{i,j} B_{i,j}(z) \frac{\partial^{2} f}{\partial z_{i} \partial z_{j}}]p(z,t|y,t')+\\
# \int dzf(z){\int' dx[W(x|z,t)p(x,t|y,t')-W(x|z,t)p(z,t|y,t')]}$ (2)
# 
# となります。ただしここで、
# 
# $\lim_{\varepsilon \rightarrow 0} \int_{|x-z|>\varepsilon dx F(x,z) \equiv \int'dxF(x,z)}$
# 
# とおきました。

# (2)について、左辺は時間微分を中に、右辺は部分積分をして、
# 
# $ \int dz f(z)\partial_{t}p(z, t| y,t') = \int dz f(z) \{ - \sum_{i} \frac{\partial}{\partial z_{i}} A_{i}(z) + \sum_{i,j} \frac{1}{2}\frac{\partial^{2} f}{\partial z_{i} \partial z_{j}}B_{i,j}(z) p(z,t|y,t')+\\
# \int dx[W(z|x,t)p(x,t|y,t')-W(x|z,t)p(z,t|y,t')]\}+{\rm surface \, terms}$
# 
# 

# **右辺の導出**
# 
# $\int dz A_{i} \frac{\partial}{\partial z_{i}}f(z) = [A_{i} f(z)] - \int dz f(z) \frac{\partial}{\partial z_{i}} A_{i} $ 
# 
# $[A_{i}f(z)]$の部分は境界線の値のみを考えればいい→surface項　　
# cf. ストークスの定理　内部の積分で全部打ち消し合って境界部分だけ考えればいい的な

# 両辺で$\int dzf(z)$は共通しているので、
# 
# **Chapman Kolmogorov方程式の微分形**  
# $\partial_{t}p(z,t|y,t') = - \sum_{i} \frac{\partial}{\partial z_{i}} A_{i}(z) p(z,t|y,t')+ \sum_{i,j} \frac{1}{2}\frac{\partial^{2} f}{\partial z_{i} \partial z_{j}}B_{i,j}(z) p(z,t|y,t')+\\
# \int dx[W(z|x,t)p(x,t|y,t')-W(x|z,t)p(z,t|y,t')]\}$

# ### 解釈
# #### Jump Processs: マスター方程式
# $A_{i}=B_{i,j}=0$のとき、マスター方程式を得ます。
# 
# $\partial_{t} (z,t|y,t') =\int dx [W(z,|x,t)p(z,t|y, t') - W(x|z,t)p(z,t|y,t')]$
# 
# $\Delta t$の1乗オーダーで近似すると **?**
# 
# $p(z, t+\Delta t| y, t) = \delta(y-z) [1 - \int dxW(x|y,t)\Delta t] + W(z|y,t)\Delta t$
# 
# つまり、$z=y$となる確率が第１項、そうでない確率が第２項、ということです。
# 
# #### Diffusion Process: Fokker-Planck方程式
# $W(z|x,t)=0$のとき、Fokker-Planck方程式になります。
# 
# $\partial_{t} (z,t|y,t') =  - \sum_{i} \frac{\partial}{\partial z_{i}} A_{i}(z) p(z,t|y,t')+ \sum_{i,j} \frac{1}{2}\frac{\partial^{2} f}{\partial z_{i} \partial z_{j}}B_{i,j}(z) p(z,t|y,t')$
# 
# これは、数学的には拡散過程です。第１項がドリフト項、第２項が拡散項です。(**後で拡散方程式について学ぶ**)
# 
# #### 決定論的過程 Liouville方程式
# Chapman-Kolomogorov方程式の微分形の第１項のみが非ゼロの場合、Liouville方程式の特別な場合を導けます。これは完全に決定論的な運動を記述することができます。

# In[ ]:



