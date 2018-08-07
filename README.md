# Nonlinear-Systems
A program that finds equilibrium immiscibility limits for a liquid-liquid system that obeys the Margules activity coefficient model.


The objective of this project is to write a program that takes an array of initial guesses of compositions x1alpha and x1beta, respectively, and predict the real composition (equilibrium immiscibility limits) of two partially miscible liquids of the first component in the α and β phases of the miscible liquid mixture.

The liquid-liquid system in this projects is assumed to obey the Margules activity coefficient model with A = 3 and B = 2.

The program uses newton raphson to solve a system of non linear equations: 
𝑓 = ln(𝛾_1_α/𝛾_1_)+ ln(x_1_α/x_1_β)
and 
𝑔 = ln(𝛾_2_α/𝛾_2_β)+ ln((1-x_1_α)/(1-x_1_β))

These non linear equations are based on the the non-ideal behavior for a given liquid-liquid system that are described by the Margules activity coefficient model:
ln(𝛾_1)= (𝑥_2)^2*(A+2(B-A)𝑥_1);
ln(𝛾_2) = (𝑥_1)^2*(B + 2(A-B)𝑥_2);
where 𝛾_1 𝑎𝑛𝑑 𝛾_2 are the activity coefficients of liquid 1 and liquid 2, respectively. 𝑥_1 and 𝑥_2 are the mole fractions of each liquid. A and B are coefficients that depend on the temperature only.

The choice of the objective functions 𝑓 = ln(𝛾_1_α/𝛾_1_β)+ ln(x_1_α/x_1_β) and 𝑔 = ln(𝛾_2_α/𝛾_2_β)+ ln((1-x_1_α)/(1-x_1_β)) comes from the change of the Gibbs free energy of mixing equation in a Margules fit model which is a measure of spontaneity of the mixing . The values of Δ G_mix/RT should reach their minima at equilibrium according to thermodynamic laws.
ΔG_mix/RT=(𝑥_1)ln(𝑥_1.𝛾_1) + (𝑥_2)ln(𝑥_2.𝛾_2) is the model of the energy of mixing.

In this model if ΔG_mix is less than 0 it means the forward reaction is favored (mixing is spontaneous). If Δ G_mix is larger than 0 it means that splitting between the two phases occurs (backward reaction is favored). At equilibrium where the rate of the backward reaction is equal to the rate of the forward reaction happens when Δ G_mix = 0.

(x_1_α.𝛾_1_α)=(x_1_β.𝛾_1_β) AND (x_2_α.𝛾_2_α)=(x_2_β.𝛾_2_β) 

If we take the natural log to both sides for both equations and rearrange them, we get f=0 and g=0.


The Newton-Raphson method was used to solve the nonlinear system of 2 equations and 2 unknowns for x_1_α and x_1_β. This was done using Gauss-Jordan elimination with scaled row pivoting and warn the user if the RCON number obtained is near or smaller than 0.001. 

In this program, the jacobian was found using the ‘jaco’ function. After that it was required to solve for [J^-1][F]. The preferred method of solving this system was using the gauss-jordan function ‘gsrp’ with scaled row pivoting on the system J x=F instead of finding the inverse. The reason behind this, is in order to find the inverse of a matrix, the operations done on the matrix and the identity matrix are the same as the operations done on the matrix in the elimination process. In other words, while we are finding the inverse we are solving the system at the same time. Thus, our choice of using gauss-jordan elimination is based on the importance of the efficiency of our program. So that the extra step of multiplying the inverse of the matrix by the b matrix won’t be needed.

The initial guess fed to the program is an important criterion for the success of the program, because the newton raphson is the most dependent function on the initial guess. The values of x1alpha and x1beta should be apart by a reasonable value that makes the jacobian obtained at their respective values solvable with high accuracy. Another reason for the importance of the initial guess, is newton raphson could diverge and give complex values or values outside the range of 0 and 1.
For example, if you feed the program an initial guess of 0.07 and 0.8 for x1alpha and x1beta, respectively, newton raphson converges and gives values of 0.0719 for x1alpha and 0.79 for x1beta, respectively.
On the other hand, if we feed the program 0.2 for x1alpha and 0.9 for x1beta, at one of the iterations the calculated jacobian will have an RCON number<0.001 which shows as an error when we run the program. This means that using newton raphson with the obtained jacobian won’t give accurate approximations of the new x1alpha and x1beta.
Finally, if we use x1alpha=0 and x1beta=0.9 for example, newton raphson will diverge from the first iteration.

HOW TO RUN?

Run script YousefKhalilProject4(xg) where xg is the initial guess of composition

