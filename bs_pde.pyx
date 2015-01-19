import numpy as np


def black_pde(S,K,v,t,r,opttype,nt,ns,smax):
    """Calculates Simple Black Scholes theo via PDE 
         - S = Spot, K= Strike
         - v = vol,  r= interest rate
         - opttype = 1 for call -1 for put
         - nt = # of time steps
         - ns = # of spot steps
         - smax = Max underlying price (for boundary)
         
        Forward difference in time, central difference in space, Order O^2
    """
    s = np.linspace(0,smax,ns)
    dt = 1/(ns**2 * v**2)
    nt = int(t/dt)
    V = np.zeros((ns,nt))
    
    #IC + BC
    
    if opttype == 1:
        V[-1,-1] = np.max(S-K,0) #intrinsic for call
        V[:,0] = np.maximum(s-K,0)
        V[0,:] = 0
        V[-1,1:] = smax - K*np.exp(-r * np.linspace(0.0, t, nt-1))
    elif opttype == -1:
        V [-1,-1] = np.max(K - S,0) #intrinsic for put
        V[:,0] = np.maximum(K-s,0)
        V[0,:] = K*np.exp(-r * np.linspace(0.0, t, nt))
        V[-1,1:] = 0
    
    #step through
    
    dex_slice = np.arange(1,ns-1)
    
    c1 = 0.5 * dt * (v**2 * dex_slice**2)
    c2 = 1 - dt * ( v**2 * dex_slice**2 + r)
    c3 = 0.5 * dt * (v**2 * dex_slice**2)
    
    
    for n in range(1,nt):       
        V[1:-1,n] = c1 * V[0:-2,n-1] + c2 * V[1:-1,n-1] + c3 * V[2:,n-1]
        V[0,n] = (1-r*dt)* V[0,n-1] #2nd to last BC
        V[-1,n] = 2*V[-2,n] - V[-3,n] #last BC (unclear if I really want this..)
    V = np.fliplr(V)
    
    return V
    