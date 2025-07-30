#include "ccp-1d.hpp"
#include "Const.hpp"
#include <iostream>

void FVM::Face::getFlux()
{
    using namespace Const;
    using namespace Tools;

    for (int iT = 0; iT < numT; ++iT)
    {
        // Obtain the gradients of variables at faces based on the finite difference method
        const scalar dNe    = getDiff(cellL.Ne(iT),  cellR.Ne(iT));
        const scalar dNi    = getDiff(cellL.Ni(iT),  cellR.Ni(iT));
        const scalar dTe    = getDiff(cellL.Te(iT),  cellR.Te(iT));
        const scalar dPhi   = getDiff(cellL.Phi(iT), cellR.Phi(iT));

        // Get upwinded variables based on direction of electric potential 
        const scalar NeC    = getUpwind(NeL(iT), NeR(iT), iT);
        const scalar NiC    = getDownwind(NiL(iT), NiR(iT), iT);
        const scalar TeC    = getUpwind(TeL(iT), TeR(iT), iT);  
        const scalar PhiC   = getUpwind(PhiL(iT), PhiR(iT), iT);

        // Electric field intensity at the face
        Ef(iT)              = - dPhi;
        // Drift-diffusion flux for electron and ion 
        fluxNe(iT)          = - De * dNe + muE * NeC * dPhi;
        fluxNi(iT)          = - Di * dNi - muI * NiC * dPhi;
        fluxPhi(iT)         = - Ef(iT) * eps0;
        // Energy transport equation for electron
        // This flux is not good. Next time we need to use fluxNe
        const scalar qE    = - 5.0 / 3.0 * De * NeC * dTe 
                           +   5.0 / 3.0 * TeC * fluxNe(iT);
        const scalar dPhiJ = - e * PhiC * fluxNe(iT);
        const scalar dJ    = e * fluxNe(iT);

        fluxEe(iT)         = qE + dPhiJ;
        fluxEeJoule(iT)    = dJ;
    }
}


void FVM::Face::getFluxLeftBC()
{
    using namespace Const;
    scalar Je;
    for (label iT = 0; iT < numT; ++iT)
    {
        // Electric field at the face
        Ef(iT)             = - (cellR.Phi(iT) - PhiL(iT)) / dist;

        // Coefficient related to the thermal velocity
        const scalar k     = 0.25 * sqrt(TeR(iT));

        // Restriction towards flux of electron flux
        Je                 = (Ef(iT) > 0 ? - muE * NeR(iT) * Ef(iT) : 0.0);

        // Fluxes
        fluxNi(iT)         = (Ef(iT) <= 0 ? muI * NiR(iT) * Ef(iT) : 0.0);
        fluxNe(iT)         = (- k  * NeR(iT)  - Gam * fluxNi(iT) + Je);
        fluxPhi(iT)        = - Ef(iT) * eps0;

        const scalar qE    = 5.0 / 3.0 *  TeR(iT) * fluxNe(iT);
        const scalar phiJ  = - e * PhiR(iT) * fluxNe(iT);
        const scalar dJ    = e * fluxNe(iT);

        fluxEe(iT)         = qE + phiJ;
        fluxEeJoule(iT)    = dJ;
    }
}


void FVM::Face::getFluxRightBC()
{
    using namespace Const;
    scalar Je;

    for (int iT = 0; iT < Const::numT; ++iT)
    {
        // Electric field at the face        
        Ef(iT)             = - (PhiR(iT) - cellL.Phi(iT)) / dist; 

        // Coefficient related to the thermal velocity
        const scalar k     = 0.25 * sqrt(TeL(iT));

        Je                 = (Ef(iT) < 0 ? -muE * NeL(iT) * Ef(iT) : 0.0);

        // Fluxes
        fluxNi(iT)         = Ef(iT) >= 0 ? muI * NiL(iT) * Ef(iT) : 0.0;
        fluxNe(iT)         =  k * NeL(iT) - Gam * fluxNi(iT) + Je;
        fluxPhi(iT)        = - Ef(iT) * eps0;

        const scalar qE    = 5.0 / 3.0 * TeL(iT) * fluxNe(iT);

        const scalar phiJ  = - e * PhiL(iT) * fluxNe(iT);
        const scalar dJ    = e * fluxNe(iT);

        fluxEe(iT)          = qE + phiJ;
        fluxEeJoule(iT) = dJ;
    }
}