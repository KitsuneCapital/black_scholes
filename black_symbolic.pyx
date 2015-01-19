import sympy as sy
from sympy.statistics import Normal as syNormal

S, K, vol, dte = sy.symbols('S,K,vol,dte',positive=True)
r,cp = sy.symbols('r,cp')

T = dte / 260.
N = syNormal(0.0, 1.0)
 
d1 = (sy.ln((S+.0000001) / K) + (0.5 * vol ** 2) * T) / (vol * sy.sqrt(T))
d2 = d1 - vol * sy.sqrt(T)
 
theo = sy.exp(-r * T) * (cp *S * N.cdf(cp*d1) - cp * K  * N.cdf(cp*d2))
 
#Black TV
tv_ = sy.lambdify((S, K, vol, dte, r,cp),theo)
 
#1st Order Greeks
delta_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(S))
vega_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(vol)/100.) #WATCH UNITS
theta_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(dte))
rho_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(r))
 
#2nd Order Greeks
gamma_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(S,S))
vanna_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(S,vol)/100.)
vomma_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(vol,vol)/1e4) #IN TICKS
charm_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(S,dte)) #DELTA DECAY
 
#3rd Order -- Who cares about anything about dGamma?
speed_ = sy.lambdify((S, K, vol, dte, r,cp),theo.diff(S,3))