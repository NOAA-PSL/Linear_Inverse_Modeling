'''
File Name    : LIM_analysis
Function Name: LIM  

Author       : Meg D. Fowler 
Date         : 8 Jan 2020 

Summary      : This function builds a linear inverse model (LIM) given a certain set of data. 
               That data should be supplied without an annual cycle. NOTE: there is no 
               check included here to make sure LIM is an appropriate choice for the dataset 
               supplied - that is left to the user. It is useful to check, for example, that
               the tau test is passed and that there are no Nyquist modes present. 

Inputs       : xDat - data to use in building LIM. Order should be [variables, time]. For example, 
                      if the goal is to use a year's worth of daily surface temperaure at 6 weather 
                      stations, the data should have dimensions of [6, 365]. 
               lag - the value supplied for Tau_0, the lag used in the lagged covariance matrix.  

Outputs      : b_alpha - values of Beta (not in a diagonal matrix as used in the calculation of L) 
               L       - the matrix of L values 
               Q       - the matrix of Q values 
               G       - the Green function 
               c0      - the contemporaneous covariance matrix 
               cT      - the lagged covariance matrix 
               normU   - the modes of G, normalized 
               v       - the adjoints of G 
               g       - the eigenvalues of G 
               periods - the periods of oscillations in LIM 
               decayT  - the decay times of the periods above 
               * NOTE: These can be easily reduced (or expanded) by altering the return statement at the bottom of the script) 

'''

def LIM(xDat,lag):  
    import numpy as np 
    from numpy import linalg as LA 

    # Take transpose of input data matrix 
    xDat_T = np.transpose(xDat) 

    # ------------------------------------------------------------------
    # STEP 1: Compute the lagged and contemporaneous covariance matrices 
    sizes = np.shape(xDat)    #Get size of matrix to determine how many data points and how many time records to consider 
    nDat = sizes[0]
    nT   = sizes[1]
  
    #Get the value of the data (xDat) at the specified lag to use in computing the lagged covariance matrix 
    xLagged = np.full([nDat,nT-lag],np.nan)  #Initialize matrix full of NaNs
    for iT in range(nT-lag):                 #Get the value of the data at the specified lag
        xLagged[:,iT] = xDat[:,iT+lag]

    # Initialize matrices full of NaNs 
    c0 = np.full([nDat, nDat], np.nan)    #Initialize matrix full of NaNs
    cT = np.full([nDat, nDat], np.nan)    #Initialize matrix full of NaNs
    
    # Compute covariance matrices for each data point 
    for iR in range(nDat):
        for iC in range(nDat):
            # Contemporaneous covariance matrix:
            c0[iR,iC] = np.nansum(xDat[iR,:]*xDat_T[:,iC]) / np.nansum(np.isfinite(xDat[iR,:]*xDat_T[:,iC]))
            # Lagged covariance matrix:
            cT[iR,iC] = np.nansum(xLagged[iR,:]*xDat_T[:-lag,iC]) / np.nansum(np.isfinite((xLagged[iR,:]*xDat_T[:-lag,iC])))

    # --------------------------------------------------------------------
    # STEP 2: Compute the Green function, defining its eigen values and vectors 
    
    G = cT.dot(LA.inv(c0))    #The Green function is defined as the product between covariance matrices 

    # Define the modes (u) and eigen-values (g) of G
    g, u = LA.eig(G)

    iSort = g.argsort()[::-1]    #Sort the eigen values and vectors in order 
    g     = g[iSort]
    u     = u[:,iSort] 

    # Define the adjoints (v) based on the transpose of G 
    eigVal_T, v = LA.eig(np.transpose(G))
    iSortT      = eigVal_T.argsort()[::-1]
    eigVal_T    = eigVal_T[iSortT]
    v           = v[:,iSortT] 
   
    # But modes should ultimately be sorted by decreasing decay time (i.e., decreasing values of 1/beta.real) 

    # Compute Beta  
    b_tau   = np.log(g)
    b_alpha = b_tau/lag

    # Sort data by decreasing decay time 
    sortVal = -1/b_alpha.real              #Decay time 

    iSort2 = sortVal.argsort()[::-1]      #Sorted indices 
    u      = u[:,iSort2]
    v      = v[:,iSort2]
    g      = g[iSort2]
    b_alpha = b_alpha[iSort2]

    # Make diagonal array of Beta (values should be negative)
    beta = np.zeros((nDat, nDat), complex)
    np.fill_diagonal(beta, b_alpha)

 
    #Need to normalize u so that u_transpose*v = identitity matrix, and u*v_transpose = identity matrix as well 
    normFactors = np.dot(np.transpose(u),v)
    normU       = np.dot(u,LA.inv(normFactors))

    # --------------------------------------------------------------------
    # STEP 3: Compute L and Q matrices 

    # Compute L matrix as normU * beta * v_transpose 
    L = np.dot(normU, np.dot(beta, np.transpose(v)))

    # Compute Q matrix 
    Q_negative = np.dot(L, c0) + np.dot(c0, np.transpose(L))
    Q = -Q_negative 

    # Also define the periods and decay times 
    periods = (2 * np.pi) / b_alpha.imag 
    decayT  = -1 / b_alpha.real 

    # --------------------------------------------------------------------
    # RETURN statement 
    return(b_alpha, L, Q, G, c0, cT, normU, v, g, periods, decayT) 


