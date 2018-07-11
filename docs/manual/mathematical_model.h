/**

@page model Equations

@tableofcontents

@section Governing equations

The experimental data concerning the relaxation times of different processes in reacting mixtures show that in many cases of practical interest the following relation is valid:

\f$\tau_{el}<\tau_{rot}\ll\tau_{vibr}<\tau_{react}\sim\theta\f$

\f$\tau_{el}\f$, \f$\tau_{rot}\f$, \f$\tau_{vibr}\f$, \f$\tau_{react}\f$ are respectively the mean times of translational, rotational, vibrational relaxation and chemical reactions and $\theta$ is the macroscopic gasdynamic time. Translational energy distribution is known to equilibrate fast and the rotational relaxation time is of the same order as the translational one and much smaller in comparison to the vibrational and chemical relaxation time. Therefore, processes of translational and rotational relaxation may be considered as rapid processes and on the contrary vibrational and chemical relaxation as the slow ones. The mean time of slow processes is comparable with the macroscopic time and these processes are strongly non-equilibrium. The condition given in Eq. [] provides the so-called level approach in non-equilibrium gas dynamics which describes the simultaneous processes of the vibrational and chemical relaxation. In this case the macroscopic conservation equations for mass, momentum and total energy should be considered together with the equations for level populations of different chemical species since no quasi-stationary vibrational distributions exist. The zero and first order one-particle distribution functions for mixture molecules over the chemical species \f$c\f$, vibrational \f$i\f$ and rotational \f$j\f$ energy levels in the phase space of the velocity \f$\boldsymbol{u}\f$, coordinates \f$\boldsymbol{r}\f$ and time \f$t\f$, respectively, \f$f_{cij}^{\left(0\right)}\left(\boldsymbol{r},\boldsymbol{u},t\right)\f$ and \f$f_{cij}^{\left(1\right)}\left(\boldsymbol{r},\boldsymbol{u},t\right)\f$, are obtained under the form of Maxwell-Boltzmann distribution over velocities and rotational energy depending on the non-equilibrium vibrational level populations and chemical species concentrations by: 

\f$f_{cij}^{\left(0\right)}=\left(\frac{m_{c}}{2\pi kT}\right)^{\frac{3}{2}}\frac{n_{ci}s_{j}^{ci}}{Z_{rot}^{ci}}\exp\left(-\frac{m_{c}c_{c}^{2}}{2kT}-\frac{\varepsilon_{j}^{ci}}{kT}\right)\f$

\f$f_{cij}^{\left(1\right)}=f_{cij}^{\left(0\right)}\left(-\frac{1}{n}\boldsymbol{A}_{cij}\cdot\nabla\ln T-\frac{1}{n}\sum_{dk}\boldsymbol{D}_{cij}^{dk}\cdot\boldsymbol{d}_{dk}-\frac{1}{n}\boldsymbol{B}_{cij}:\nabla\boldsymbol{v}-\frac{1}{n}F_{cij}\nabla\cdot\boldsymbol{v}-\frac{1}{n}G_{cij}\right)\f$

where \f$m_{c}\f$ is the molecular mass of chemical species \f$c\f$, \f$k=1.3807\text{·}10^{-23}\frac{J}{K}\f$ is the Boltzmann constant, \f$T\f$ is the static gas temperature, \f$n_{ci}\f$ are the non-equilibrium vibrational level populations, \f$\boldsymbol{c}_{c}=\boldsymbol{u}_{c}-\boldsymbol{v}_{c}\f$ is the peculiar velocity, \f$\boldsymbol{v}\f$ is the macroscopic gas velocity vector, \f$s_{j}^{ci}=2j+1\f$, \f$\varepsilon_{j}^{ci}\f$ are the statistic weight (degeneracy) and energy of the \f$j\f$th rotational state corresponding to the vibrational level $i$ of chemical species \f$c\f$, \f$Z_{rot}^{ci}\f$ is the rotational partition function, \f$n\f$ is the total number density, \f$\boldsymbol{d}_{ci}\f$ are the diffusion driving forces for each chemical and vibrational species. Functions \f$\boldsymbol{A}_{cij}\f$, \f$\boldsymbol{B}_{cij}\f$, \f$\boldsymbol{D}_{cij}^{dk}\f$, \f$F_{cij}\f$ and \f$G_{cij}\f$ are found from the linear integral equations.
A closed set of macroscopic equations for the macroscopic parameters \f$n_{ci}\left(\boldsymbol{r},t\right), \boldsymbol{v}\left(\boldsymbol{r},t\right)\f$ and \f$T\left(\boldsymbol{r},t\right)\f$ is obtained and consists of the conservation equations for the momentum and total energy coupled with the equations of detailed vibration-chemical kinetics for the vibrational level populations. It reads:

\f$\frac{dn}{dt}+n\nabla\cdot\boldsymbol{v}=0\f$

\f$\rho\frac{d\boldsymbol{v}}{dt}+\nabla\cdot\boldsymbol{P}=0\f$

\f$\rho\frac{dU}{dt}+\nabla\cdot\boldsymbol{q}+\boldsymbol{P}:\nabla\boldsymbol{v}=0\f$

\f$\frac{dn_{ci}}{dt}+n_{ci}\nabla\cdot\boldsymbol{v}+\nabla\cdot\left(n_{ci}\boldsymbol{V}_{ci}\right)=R_{ci}\f$ 

\f$\qquad c=1,...,L,\quad i=0....,L_{c}\f$

where \f$\rho\f$ is the density of the mixture and the population of molecular species \f$c\f$ for the vibrational level \f$i\f$ per unit volume (or the number density of molecules of species $c$ at the vibrational level \f$i\f$) is given by the expression:

\f$n_{ci}\left(\boldsymbol{r},t\right)=\sum_{j}\int f_{cij}\left(\boldsymbol{r},\boldsymbol{u},t\right)d\boldsymbol{u}_{c}\f$ 

whereas the number density of molecular species \f$c\f$ is: 

\f$n_{c}\left(\boldsymbol{r},t\right)=\sum_{ij}\int f_{cij}\left(\boldsymbol{r},\boldsymbol{u},t\right)d\boldsymbol{u}_{c}=\sum_{i}n_{ci}\f$

and the number density of a gas mixture is given by:

\f$n\left(\boldsymbol{r},t\right)=\sum_{cij}\int f_{cij}\left(\boldsymbol{r},\boldsymbol{u},t\right)d\boldsymbol{u}_{c}=\sum_{c}n_{c}\f$

\f$\boldsymbol{P}\f$ is the stress tensor, \f$\boldsymbol{q}\f$ the total energy flux and \f$\boldsymbol{V}_{ci}\f$ the diffusion velocities of molecular species \f$c\f$ on the vibrational level \f$i\f$. The specific total energy \f$U\f is given by the following expression:

\f$\rho U\left(\boldsymbol{r},t\right)=\frac{3}{2}nkT+\sum_{ci}\langle\varepsilon^{ci}\rangle n_{ci}+\sum_{ci}\varepsilon_{i}^{c}n_{ci}+\sum_{c}\varepsilon_{c}n_{c}\f$

\f$\langle\varepsilon^{ci}\rangle\f$ is the average rotational energy of molecule of species \f$c\f$ at vibrational level \f$i\f$; \f$\varepsilon_{i}^{c}\f$ is the vibrational energy of level $i$ of molecular species \f$c\f$;
\f$\varepsilon_{c}\f$ is the formation energy of species \f$c\f$. The average rotational energy is defined as: 

\f$\langle\varepsilon^{ci}\rangle_{rot}=\frac{1}{\sigma}\frac{1}{Z_{rot}^{ci}}\sum_{j}s_{j}^{ci}\varepsilon_{j}^{ci}\exp\left(-\frac{\varepsilon_{j}^{ci}}{kT}\right)\f$

where

\f$\sigma\f$ is a symmetry factor, equal to 1 for heteronuclear and 2 for homonuclear molecules and the rotational partition function:

\f$Z_{rot}^{ci}=\frac{1}{\sigma}\sum_{j}s_{j}^{ci}\exp\left(-\frac{\varepsilon_{j}^{ci}}{kT}\right)\f$

The generalized Chapman-Enskog method allows to express, in any approximation, the transport and relaxation terms as functions of the chosen macroscopic parameters

and thus to close completely the set of governing equations. Specifically, it is known that the zero-order approximation describes detailed state-to-state vibrational and chemical kinetics in an inviscid non-conductive gas
mixture flow (Euler) while the first-order approximation corresponds to the Navier-Stokes approximation. 

The first order transport terms are determined by the first order distribution functions. The pressure tensor \f$\boldsymbol{P}\f$ has been obtained in the following form:

\f$\boldsymbol{P}=\left(p-p_{rel}\right)\boldsymbol{I}-2\eta\boldsymbol{S}-\zeta\nabla\cdot\boldsymbol{v}\boldsymbol{I}\f$

here \f$\boldsymbol{I}\f$ is the unit tensor, \f$\boldsymbol{S}\f$ is the tensor of deformation velocities, \f$p_{rel}\f$ is the relaxation pressure, \f$\eta\f$ and \f$\zeta\f$ are the coefficients of shear and bulk viscosity, respectively:

\f$\eta=\frac{kT}{10}[\boldsymbol{B},\boldsymbol{B}],\quad\zeta=kT[F,F],\quad p_{rel}=kT[F,G]\f$

The bracket integrals \f$\left[A,B\right]\f$ which appear in are introduced in and depend on the cross-sections of the most frequent collisions, i.e. the elastic collisions and those
leading to the rotational energy exchange.

The diffusion velocity \f$\boldsymbol{V}_{ci}\f$ of molecular components \f$c\f$ at the vibrational level \f$i\f$ is specified by:

\f$\boldsymbol{V}_{ci}=-\sum_{dk}D_{cidk}\boldsymbol{d}_{dk}-D_{Tci}\nabla\ln T\f$ 

where \f$D_{cidk}\f$ and \f$D_{Tci}\f$ are the multicomponent diffusion and thermal diffusion coefficients for each chemical and vibrational species:

\f$D_{cidk}=\frac{1}{3n}[\boldsymbol{D}^{ci},\boldsymbol{D}^{dk}],\quad D_{Tci}=\frac{1}{3n}[\boldsymbol{D}^{ci},\boldsymbol{A}]\f$

and the species specific driving forces are defined as: 

\f$\boldsymbol{d}_{ci}=\nabla\left(\frac{n_{ci}}{n}\right)+\left(\frac{n_{ci}}{n}-\frac{\rho_{ci}}{\rho}\right)\nabla\ln p\f$

where \f$\rho_{ci}\f$ is the mass density of molecular species \f$c\f$ on the vibrational level \f$i\f$ and \f$p=nkT\f$ is the static gas mixture pressure. The expression for the total heat flux in the first order
approximation takes the following form:

\f$q=-\lambda^{'}\nabla T-p\sum_{ci}D_{Tci}\boldsymbol{d}_{ci}+\sum_{ci}\left(\frac{5}{2}kT+\langle\varepsilon_{j}^{ci}\rangle_{rot}+\varepsilon_{i}^{c}+\varepsilon_{c}\right)n_{ci}\boldsymbol{V}_{ci}\f$

and the thermal conductivity coefficient are found in the form:

\f$\lambda^{'}=\lambda_{tr}+\lambda_{rot}=\frac{k}{3}[\boldsymbol{A},\boldsymbol{A}]\f$
 
It is important to emphasize that vibrational energy does not contribute to the thermal conductivity in this approach. We also observe that the expressions for the heat flux and diffusion velocity contain not
only the gradient of the gas temperature but also the gradients of all level populations of different molecular species and number densities of atoms with corresponding diffusion and thermal diffusion coefficients.

Equations describe the non-equilibrium vibrational and chemical kinetics in a gas flow. The source terms which appear therein, characterize the variation of the vibrational level populations
and atomic number densities caused by vibrational energy exchanges and chemical reactions and are expressed via the integral operators of slow processes proceeding on the gasdynamic time scale:

\f$R_{ci}=R_{ci}^{vibr}+R_{ci}^{react}\f$

We finally note that the vibrational level populations are included to the set of main macroscopic parameters and the equations for their calculation are coupled to the equations of gas dynamics. For this
reason we refer to vibrational state-to-state approach.
*/
