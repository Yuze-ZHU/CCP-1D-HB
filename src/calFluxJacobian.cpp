#include "ccp-1d.hpp"
#include "Const.hpp"
#include <cstring>
#include <cmath> 
#include <petsc.h>
#include <iomanip>

void FVM::Face::getFluxJacobianFCI()
{
    using namespace Const;
    for (label iT = 0; iT < numT; ++numT)
    {
        jacobianFluxL[iT].setZero();
        jacobianFluxR[iT].setZero();

        const scalar dPhi   = getDiff(cellL.Phi(iT), cellR.Phi(iT));
        const scalar dNe    = getDiff(cellL.Ne(iT),  cellR.Ne(iT));
        if(PhiR(iT) > PhiL(iT))
        {
            jacobianFluxL[iT](0, 0) = De / dist + muE * dPhi;
            jacobianFluxR[iT](0, 0) = - De / dist;

            jacobianFluxL[iT](0, 3) = - muE * cellL.Ne(iT) / dist;
            jacobianFluxR[iT](0, 3) = muE * cellL.Ne(iT) / dist;

            jacobianFluxL[iT](1, 1) = Di / dist;
            jacobianFluxR[iT](1, 1) = - Di / dist - muI * dPhi;

            jacobianFluxL[iT](1, 3) = muI * cellR.Ni(iT) / dist;
            jacobianFluxR[iT](1, 3) = - muI * cellR.Ni(iT) / dist;

            jacobianFluxL[iT](2, 0) = 5.0 / 3.0 * De / dist * (cellL.Ee(iT) * cellR.Ne(iT) 
                                    / std::pow(cellL.Ne(iT), 2.0) - cellR.Ee(iT) / cellR.Ne(iT));
            jacobianFluxL[iT](2, 0) -= e * cellL.Phi(iT) * (De / dist + muE * dPhi);
            jacobianFluxR[iT](2, 0) = 5.0 / 3.0 * De / dist * (cellR.Ee(iT) * cellL.Ne(iT) 
                                    / std::pow(cellR.Ne(iT), 2.0) - cellL.Ee(iT) / cellL.Ne(iT));
            jacobianFluxR[iT](2, 0) += e * cellL.Phi(iT) * De / dist;

            jacobianFluxL[iT](2, 2) = 5.0 / 3.0 * De / dist * (2.0 * cellL.Ne(iT) - cellR.Ne(iT)) 
                                    / cellL.Ne(iT) + 5.0 / 3.0 * muE * dPhi;
            jacobianFluxR[iT](2, 2) = - 5.0 / 3.0 * De / dist * cellL.Ne(iT) / cellR.Ne(iT);

            jacobianFluxL[iT](2, 3) = - 5.0 / 3.0 * muE * cellL.Ee(iT) / dist + e * De * dNe;
            jacobianFluxL[iT](2, 3) -= e * muE * cellR.Ne(iT) * (cellR.Phi(iT) - 2.0 * cellL.Phi(iT)) / dist;
            jacobianFluxR[iT](2, 3) = 5.0 / 3.0 * muE * cellL.Ee(iT) / dist;
            jacobianFluxR[iT](2, 3) -= e * muE * cellL.Ne(iT) * cellL.Phi(iT) / dist;

            jacobianFluxL[iT](3, 3) = -eps0 / dist;
            jacobianFluxR[iT](3, 3) = eps0 / dist;

        }
        else
        {
            jacobianFluxL[iT](0, 0) = De / dist;
            jacobianFluxR[iT](0, 0) = - De / dist + muE * dPhi;
            jacobianFluxL[iT](0, 3) = - muE * cellR.Ne(iT) / dist;
            jacobianFluxR[iT](0, 3) = muE * cellR.Ne(iT) / dist;


            jacobianFluxL[iT](1, 1) = Di / dist - muI * dPhi;
            jacobianFluxR[iT](1, 1) = - Di / dist;

            jacobianFluxL[iT](1, 3) = muI * cellL.Ni(iT) / dist;
            jacobianFluxR[iT](1, 3) = - muI * cellL.Ni(iT) / dist;

            jacobianFluxL[iT](2, 0) = 5.0 / 3.0 * De / dist * (cellR.Ee(iT) / cellR.Ne(iT) 
                                - cellL.Ee(iT) * cellR.Ne(iT) / std::pow(cellL.Ne(iT), 2.0));
            jacobianFluxL[iT](2, 0) -= e * De * cellR.Phi(iT) / dist;

            jacobianFluxR[iT](2, 0) = 5.0 / 3.0 * De / dist * (cellL.Ee(iT) / cellL.Ne(iT) 
                                - cellR.Ee(iT) * cellL.Ne(iT) / std::pow(cellR.Ne(iT), 2.0));
            jacobianFluxR[iT](2, 0) -= e * cellR.Phi(iT)  * (-De / dist + muE * dPhi);

            jacobianFluxL[iT](2, 2) = 5.0 / 3.0 * De / dist * cellR.Ne(iT) / cellL.Ne(iT);
            jacobianFluxR[iT](2, 2) = - 5.0 / 3.0 * De / dist * (2.0 * cellR.Ne(iT) 
                                - cellL.Ne(iT)) / cellR.Ne(iT) + 5.0 / 3.0 * muE * dPhi;

            jacobianFluxL[iT](2, 3) = -5.0 / 3.0 * muE * cellR.Ee(iT) / dist;
            jacobianFluxL[iT](2, 3) += e * muE * cellR.Phi(iT) * cellR.Ne(iT) / dist;
            jacobianFluxR[iT](2, 3) = 5.0 / 3.0 * muE * cellR.Ee(iT) / dist + e * De * dNe;
            jacobianFluxR[iT](2, 3) -= e * muE * cellR.Ne(iT) * (2 * cellR.Phi(iT) - cellL.Phi(iT)) / dist;

            jacobianFluxL[iT](3, 3) = - eps0 / dist;
            jacobianFluxR[iT](3, 3) = eps0 / dist;
        }
    }

}


void FVM::Face::getFluxJacobianLeftBCFCI()
{
    using namespace Const;

    for (label iT = 0; iT < numT; ++iT)
    {
        const scalar dPhiNeg   = - getDiff(PhiL(iT), cellR.Phi(iT));
        jacobianFluxL[iT].setZero();
        jacobianFluxR[iT].setZero();

        if (PhiR(iT) > PhiL(iT))
        {
            jacobianFluxR[iT](0, 0) = - 1.0 / 8.0 * sqrt(cellR.Ee(iT) / cellR.Ne(iT));

            jacobianFluxR[iT](0, 1) = - Gam * muI * dPhiNeg;

            jacobianFluxR[iT](0, 2) = - 1.0 / 8.0 * sqrt(cellR.Ne(iT) / cellR.Ee(iT));

            jacobianFluxR[iT](0, 3) = Gam * muI *cellR.Ni(iT) / dist;

            jacobianFluxR[iT](1, 1) = muI * dPhiNeg;

            jacobianFluxR[iT](1, 3) = - muI * cellR.Ni(iT) / dist;

            jacobianFluxR[iT](2, 0) = 5.0 / 24.0 * std::pow(cellR.Ee(iT), 1.5) / std::pow(cellR.Ne(iT), 1.5);
            jacobianFluxR[iT](2, 0) += 5.0 / 3.0 * Gam * muI * cellR.Ni(iT) 
                                    * cellR.Ee(iT) / std::pow(cellR.Ne(iT), 2.0) * dPhiNeg;
            jacobianFluxR[iT](2, 0) += e * cellR.Phi(iT) / 8.0 * sqrt(cellR.Ee(iT) / cellR.Ne(iT));

            jacobianFluxR[iT](2, 1) = -5.0 / 3.0 * muI * Gam * cellR.Ee(iT) / cellR.Ne(iT) * dPhiNeg;
            jacobianFluxR[iT](2, 1) += e * cellR.Phi(iT) * Gam * muI * dPhiNeg;

            jacobianFluxR[iT](2, 2) = -5.0 / 8.0 * sqrt(cellR.Ee(iT) / cellR.Ne(iT));
            jacobianFluxR[iT](2, 2) -= 5.0 / 3.0 * Gam * muI * cellR.Ni(iT) / cellR.Ne(iT) * dPhiNeg;
            jacobianFluxR[iT](2, 2) += e * cellR.Phi(iT) / 8.0 * sqrt(cellR.Ne(iT) / cellR.Ee(iT));

            jacobianFluxR[iT](2, 3) = 5.0 / 3.0 * cellR.Ee(iT) * Gam * muI * cellR.Ni(iT) /cellR.Ne(iT) / dist;
            jacobianFluxR[iT](2, 3) += e * 1.0 / 4.0 * sqrt(cellR.Ee(iT) * cellR.Ne(iT));
            jacobianFluxR[iT](2, 3) += e * Gam * muI * cellR.Ni(iT) * dPhiNeg;
            jacobianFluxR[iT](2, 3) -= e * cellR.Phi(iT) * Gam * muI * cellR.Ni(iT) / dist;

            jacobianFluxR[iT](3, 3) = eps0 / dist;
        }
        else
        {
        
            jacobianFluxR[iT](0, 0) = -1.0 / 8.0 * sqrt(cellR.Ee(iT) / cellR.Ne(iT));

            jacobianFluxR[iT](0, 2) = -1.0 / 8.0 * sqrt(cellR.Ne(iT) / cellR.Ee(iT));

            jacobianFluxR[iT](2, 0) = 5.0 / 24.0 * std::pow(cellR.Ee(iT), 1.5) / std::pow(cellR.Ne(iT), 1.5);
            jacobianFluxR[iT](2, 0) += e * cellR.Phi(iT) / 8.0 * sqrt(cellR.Ee(iT) / cellR.Ne(iT));

            jacobianFluxR[iT](2, 2) = -5.0 / 8.0 * sqrt(cellR.Ee(iT) / cellR.Ne(iT));
            jacobianFluxR[iT](2, 2) += e * cellR.Phi(iT) / 8.0 * sqrt(cellR.Ne(iT) / cellR.Ee(iT));

            jacobianFluxR[iT](2, 3) = e * 1.0 / 4.0 * sqrt(cellR.Ee(iT) * cellR.Ne(iT));

            jacobianFluxR[iT](3, 3) = eps0 / dist;
        }
    }

}

void FVM::Face::getFluxJacobianRightBCFCI()
{
    using namespace Const;

    for (label iT = 0; iT < numT; ++iT)
    {
        const scalar dPhi   = getDiff(cellL.Phi(iT), PhiR(iT));
        jacobianFluxL[iT].setZero();
        jacobianFluxR[iT].setZero();

        if (PhiR(iT) > PhiL(iT))
        {
            jacobianFluxL[iT](0, 0)  = 1.0 / 8.0 * sqrt(cellL.Ee(iT) / cellL.Ne(iT));

            jacobianFluxL[iT](0, 2)  = 1.0 / 8.0 * sqrt(cellL.Ne(iT) / cellL.Ee(iT));

            jacobianFluxL[iT](2, 0)  = - 5.0 / 24.0 * std::pow(cellL.Ee(iT), 1.5) / std::pow(cellL.Ne(iT), 1.5);
            jacobianFluxL[iT](2, 0) -= e * cellL.Phi(iT) / 8.0 * sqrt(cellL.Ee(iT) / cellL.Ne(iT));

            jacobianFluxL[iT](2, 2)  = 5.0 / 8.0 * sqrt(cellL.Ee(iT) / cellL.Ne(iT));
            jacobianFluxL[iT](2, 2) -= e * cellL.Phi(iT) / 8.0 * sqrt(cellL.Ne(iT) / cellL.Ee(iT));

            jacobianFluxL[iT](2, 3)  = - e * 1.0 / 4.0 * sqrt(cellL.Ee(iT) * cellL.Ne(iT));

            jacobianFluxL[iT](3, 3)  = - eps0 / dist;

        }
        else
        {
        
            jacobianFluxL[iT](0, 0)  = 1.0 / 8.0 * sqrt(cellL.Ee(iT) / cellL.Ne(iT));

            jacobianFluxL[iT](0, 1)  = Gam * muI * dPhi;

            jacobianFluxL[iT](0, 2)  = 1.0 / 8.0 * sqrt(cellL.Ne(iT) / cellL.Ee(iT));

            jacobianFluxL[iT](0, 3)  = - Gam * muI * cellL.Ni(iT) / dist;

            jacobianFluxL[iT](1, 1)  = -muI * dPhi;

            jacobianFluxL[iT](1, 3)  = muI * cellL.Ni(iT) / dist;


            jacobianFluxL[iT](2, 0)  = - 5.0 / 24.0 * std::pow(cellL.Ee(iT), 1.5) / std::pow(cellL.Ne(iT), 1.5);
            jacobianFluxL[iT](2, 0) -= 5.0 / 3.0 * Gam * muI * cellL.Ni(iT) 
                                     * cellL.Ee(iT) / std::pow(cellL.Ne(iT), 2.0) * dPhi;
            jacobianFluxL[iT](2, 0) -= e * cellL.Phi(iT) / 8.0 * sqrt(cellL.Ee(iT) / cellL.Ne(iT));

            jacobianFluxL[iT](2, 1)  = 5.0 / 3.0 * muI * Gam * cellL.Ee(iT) / cellL.Ne(iT) * dPhi;
            jacobianFluxL[iT](2, 1) -= e * cellL.Phi(iT) * Gam * muI * dPhi;

            jacobianFluxL[iT](2, 2)  =  5.0 / 8.0 * sqrt(cellL.Ee(iT) / cellL.Ne(iT));
            jacobianFluxL[iT](2, 2) += 5.0 / 3.0 * Gam * muI * cellL.Ni(iT) / cellL.Ne(iT) * dPhi;
            jacobianFluxL[iT](2, 2) -= e * cellL.Phi(iT) / 8.0 * sqrt(cellL.Ne(iT) / cellL.Ee(iT));

            jacobianFluxL[iT](2, 3)  = - 5.0 / 3.0 * cellL.Ee(iT) * Gam * muI * cellL.Ni(iT) / cellL.Ne(iT) / dist;
            jacobianFluxL[iT](2, 3) -= e * 1.0 / 4.0 * sqrt(cellL.Ee(iT) * cellL.Ne(iT));
            jacobianFluxL[iT](2, 3) -= e * Gam * muI * cellL.Ni(iT) * dPhi;
            jacobianFluxL[iT](2, 3) += e * cellL.Phi(iT) * Gam * muI * cellL.Ni(iT) / dist;

            jacobianFluxL[iT](3, 3)  = - eps0 / dist;
        }
    }
}

void FVM::Solver::getCsJacobian()
{
    using namespace Const;
    for (label i = 0; i < numCells; ++i)
    {
        for (label iT = 0; iT < numT; ++iT)
        {
            cells[i].jacobianCs[iT].setZero();
            cells[i].jacobianCs[iT](0, 0) = cells[i].kl(iT) * N * cells[i].vol;
            cells[i].jacobianCs[iT](1, 0) = cells[i].kl(iT) * N * cells[i].vol;
            cells[i].jacobianCs[iT](2, 0) = - cells[i].kl(iT) * Hl * N * cells[i].vol;
            cells[i].jacobianCs[iT](3, 0) = e * cells[i].vol;
            cells[i].jacobianCs[iT](3, 1) = - e * cells[i].vol;
        }
    }
}

void FVM::Solver::getCsJacobianNeg()
{
    using namespace Const;
    for (label i = 0; i < numCells; ++i)
    {
        for (label iT = 0; iT < numT; ++iT)
        {
            cells[i].jacobianCs[iT] = - cells[i].jacobianCs[iT];
        }
    }
}

void FVM::Face::getFluxJacobianNeg()
{
    using namespace Const;
    for (label iT = 0; iT < numT; ++iT)
    {
        jacobianFluxLNeg[iT].setZero();
        jacobianFluxRNeg[iT].setZero();

        jacobianFluxLNeg[iT] = - jacobianFluxL[iT];
        jacobianFluxRNeg[iT] = - jacobianFluxR[iT];
    }
}


void FVM::Solver::getFluxJouleJacobian()
{
    using namespace Const;
    for (label iT = 0; iT < numT; ++iT)
    {
        for (label i = 0; i < numCells; ++i)
        {
            jacobianFluxJouleLCell[iT][i].setZero();
            jacobianFluxJouleCCell[iT][i].setZero();
            jacobianFluxJouleRCell[iT][i].setZero();

            //- Left cell
            jacobianFluxJouleLCell[iT][i](2, 0)  = - e * cells[i].Phi(iT) * faces[i].jacobianFluxL[iT](0, 0);
            jacobianFluxJouleLCell[iT][i](2, 1)  = - e * cells[i].Phi(iT) * faces[i].jacobianFluxL[iT](0, 1);
            jacobianFluxJouleLCell[iT][i](2, 2)  = - e * cells[i].Phi(iT) * faces[i].jacobianFluxL[iT](0, 2);
            jacobianFluxJouleLCell[iT][i](2, 3)  = - e * cells[i].Phi(iT) * faces[i].jacobianFluxL[iT](0, 3);

            //- Centeral cell
            jacobianFluxJouleCCell[iT][i](2, 0)  = e * cells[i].Phi(iT) * (faces[i + 1].jacobianFluxL[iT](0, 0) 
                                                 - faces[i].jacobianFluxR[iT](0, 0));
            jacobianFluxJouleCCell[iT][i](2, 1)  = e * cells[i].Phi(iT) * (faces[i + 1].jacobianFluxL[iT](0, 1) 
                                                 - faces[i].jacobianFluxR[iT](0, 1));
            jacobianFluxJouleCCell[iT][i](2, 2)  = e * cells[i].Phi(iT) * (faces[i + 1].jacobianFluxL[iT](0, 2) 
                                                 - faces[i].jacobianFluxR[iT](0, 2));
            jacobianFluxJouleCCell[iT][i](2, 3)  = e * (faces[i + 1].fluxNe(iT) - faces[i].fluxNe(iT));
            jacobianFluxJouleCCell[iT][i](2, 3) += e * cells[i].Phi(iT) * (faces[i+1].jacobianFluxL[iT](0, 3) 
                                                 - faces[i].jacobianFluxR[iT](0, 3));
            
            //- Right cell
            jacobianFluxJouleRCell[iT][i](2, 0)  = e * cells[i].Phi(iT) * faces[i + 1].jacobianFluxR[iT](0,0);
            jacobianFluxJouleRCell[iT][i](2, 1)  = e * cells[i].Phi(iT) * faces[i + 1].jacobianFluxR[iT](0,1); 
            jacobianFluxJouleRCell[iT][i](2, 2)  = e * cells[i].Phi(iT) * faces[i + 1].jacobianFluxR[iT](0,2);
            jacobianFluxJouleRCell[iT][i](2, 3)  = e * cells[i].Phi(iT) * faces[i + 1].jacobianFluxR[iT](0,3);           
        }
    }
}
