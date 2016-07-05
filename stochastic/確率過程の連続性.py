
# coding: utf-8

# # Table of Contents
#  <p><div class="lev2"><a href="#Intro"><span class="toc-item-num">0.1&nbsp;&nbsp;</span>Intro</a></div><div class="lev2"><a href="#連続マルコフ過程の定義"><span class="toc-item-num">0.2&nbsp;&nbsp;</span>連続マルコフ過程の定義</a></div>

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
# 
# 

# In[ ]:



