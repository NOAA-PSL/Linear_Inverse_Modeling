# LinearInverseModeling
Scripts to build a linear inverse model (LIM) based on an input pickle file.

These scripts were initially created for analysis carried out in Fowler and Penland (2020; in prep): "Predictability of soil moisture in Northern California." They have, however, been modified to be as versatile as possible for any given user (at least in the case of LinearInverseModel.ipynb). 

Scripts: 
1. <b>PrepareSoilData.ipynb</b>: This case-specific file prepares data for input to the linear inverse model script. We provide this code as an example, using in-situ soil temperature and moisture observations. Data from multiple files (saved on a yearly basis) are first combined into a single array at each observting station. The average annual cycle, including the mean, is then removed linearly from each time series. These anomalies are then cast in terms of their z-scores so that they vary on similar scales. We regress the soil temperature onto moisture and subtract that signal from temperature to avoid analysis of redundant data in this particular case.   

   The process *you* take to prepare your own data may differ, but the key steps will be to center the data and remove the    annual cycle (how you do this is up to you). Then save that data as a pickle file (https://wiki.python.org/moin/UsingPickle), with an array of dimensions [nData, nTime]. 

2. <b>LinearInverseModel.ipynb</b>: This is the crux of the work. Here, a linear inverse model is built based on a user-supplied pickle file. Other inputs required are a list of stationIDs (i.e., how you would identify the different components of your multi-variable timeseries) and a value for tau_0, the lag at which LIM will estimate its parameters. The LIM is built in a step-by-step process here, with a key section (Section #5) devoted to testing the validity of applying LIM as originally specified. *This is a critical piece to consider.* If your data gives rise to Nyquist modes (see Penland 2019), the lag you've chosne needs to be reduced. If the tau test fails or if the eigenvalues of Q are not all positive, then it is possible that LIM is not applicable. 


3. <b>MSE_SoilMoistureHindcasts.ipynb</b>: This script is as generic as possible, though it is noted throughout that you may want to change a number of options. The overall purpose here is to compare the mean squared error of LIM-generated hindcasts (based on results from the above script) to errors associated with an AR1 process. This is also where we define "forecasts of opportunity" based on how strong the projection of the optimal structure is on initial conditions. 

