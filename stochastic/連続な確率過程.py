
# coding: utf-8

# # Table of Contents
#  <p><div class="lev2"><a href="#Intro"><span class="toc-item-num">0.1&nbsp;&nbsp;</span>Intro</a></div><div class="lev2"><a href="#連続マルコフ過程の定義"><span class="toc-item-num">0.2&nbsp;&nbsp;</span>連続マルコフ過程の定義</a></div><div class="lev2"><a href="#マルコフ過程の定常解"><span class="toc-item-num">0.3&nbsp;&nbsp;</span>マルコフ過程の定常解</a></div>

# ## Intro
# 確率変数X(t)の取りうる値が連続な範囲をもつかどうか、というのと、その1経路が連続な関数か、というのはまったく別の問いです。
# 
# 例えば、空気中の分子の衝突を考えたとき、その速度V(t)はどんな値もとりうるので取りうる値の範囲は連続ですが、衝突が起きた時に速度は離散的に変わるので、その経路は不連続です。しかし、その位置X(t)については連続であると考えられます。
# 
# では、連続な経路をもつマルコフ過程というのは存在するのでしょうか？

# ## 連続マルコフ過程の定義
# マルコフ過程に対して、任意の$\varepsilon > 0$について、以下を満たせば、確率1で経路は連続であるといえます。
# 
# $\lim_{\Delta t\rightarrow 0} \frac{1}{\Delta t} \int_{|x-z|>\varepsilon} dx p(x, t+\Delta t| z,t) = 0$
# 
# つまり、$\Delta t$が小さくなるより速く、時刻$\Delta t +t$において$|x-z|>\varepsilon$となるような位置にいる確率が小さくなるということです。
# 
# **例**
# i) アインシュタインの解は連続な経路をもちます。
# 
# $p(x, t+\Delta t|z,t) = (4 \pi D \Delta t)^{1/2}\exp [-(x-z)^{2}/4D \Delta t]$
# 
# ii) コーシー列は連続な経路を持ちません。
# 
# $p(x, t+\Delta t|z,t) = \frac{\Delta t}{\pi} \frac{1}{(x-z)^{2}+\Delta t^{2}}$
# 
# いずれの場合でもChapman-Kolomogorov方程式は満たします。
# 

# ## マルコフ過程の定常解
# A,B,Wが時間に対して独立だとして、どんな条件の場合に確率過程は定常過程に近づくのでしょうか？
# 
# ２つの確率過程$p_{1}$、$p_{2}$に対して、次のようなリアプノフ関数Kを考えます。（KL情報量）
# 
# $K=\int dx p_{1}(x,t) \log [\frac{p_{1}(x,t)}{p_{2}(x,t)}]$
# 
# これに対して、$\frac{dK}{dt} \leq 0$、すなわKL情報量が単調減少であること、すなわち分布が時間無限大で一致することを示すことで、定常分布が存在することを示します。
# 
# $\frac{dK}{dt} = \int dx\{\frac{\partial p_{1}}{\partial t} [\log p_{1} + 1- \log p_{2} ] - \frac{\partial p_{2}}{\partial t} [\frac{p_{1}}{p_{2}}]\}$

# ここで、Chapman Kolmogorov方程式の微分形を$\frac{\partial p_{1}}{partial t}$、$\frac{\partial p_{2}}{partial t}$代入して、各項に分けて考えると、
# 
# $(\frac{dK}{dt})_{\rm drift} = \sum_{i} \int dx \frac{\partial}{\partial x_{i}} [-A_{i}p_{i} \log (p_{1}/p_{2})]$
# 
# 
# 

# In[ ]:



