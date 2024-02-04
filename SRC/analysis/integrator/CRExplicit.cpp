/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

#include <CRExplicit.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include "fullGEN/FullGenLinSOE.h"
#include "fullGEN/FullGenLinLapackSolver.h"
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#define OPS_Export


void *    OPS_CRExplicit(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 0 && argc != 1) {
        opserr << "WARNING - incorrect number of args want CRExplicit <-updateElemDisp>\n";
        return 0;
    }
    bool updElemDisp = false;

    if (argc == 1) {
        const char* argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-updateElemDisp") == 0)
            updElemDisp = true;
    }

    theIntegrator = new CRExplicit(updElemDisp);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating CRExplicit integrator\n";
    
    return theIntegrator;
}


CRExplicit::CRExplicit( )
    : TransientIntegrator(INTEGRATOR_TAGS_CRExplicit),
    deltaT(0.0), 
    alpha1(0), Mhat(0),
    updateCount(0), initAlphaMatrices(1),
    c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    U(0), Udot(0), Udotdot(0),
    Utdothat(0)
{
    
}


CRExplicit::CRExplicit(bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_CRExplicit),
    deltaT(0.0),
    alpha1(0), Mhat(0),
    updateCount(0), initAlphaMatrices(1),
    c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    U(0), Udot(0), Udotdot(0),
    Utdothat(0)
{
    
}


CRExplicit::~CRExplicit()
{
    // clean up the memory created
    if (alpha1 != 0)
        delete alpha1;
    if (Mhat != 0)
        delete Mhat;
    if (Ut != 0)
        delete Ut;
    if (Utdot != 0)
        delete Utdot;
    if (Utdotdot != 0)
        delete Utdotdot;
    if (U != 0)
        delete U;
    if (Udot != 0)
        delete Udot;
    if (Udotdot != 0)
        delete Udotdot;
    if (Utdothat != 0)
        delete Utdothat;
}


int CRExplicit::newStep(double _deltaT)
{
    updateCount = 0;
    
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING CRExplicit::newStep() - no AnalysisModel set\n";
        return -2;
    }
    
    if (initAlphaMatrices || _deltaT != deltaT)  {
        
        // update time step increment
        deltaT = _deltaT;
        if (deltaT <= 0.0)  {
            opserr << "WARNING CRExplicit::newStep() - error in variable\n";
            opserr << "dT = " << deltaT << endln;
            return -3;
        }
        
        // get the LinearSOE and the ConvergenceTest so we can switch back later
        LinearSOE *theLinSOE = this->getLinearSOE();
        ConvergenceTest *theTest = this->getConvergenceTest();
        
        // set up the FullLinearSOE (needed to compute the alpha matrices)
        int size = theLinSOE->getNumEqn();
        FullGenLinSolver *theFullLinSolver = new FullGenLinLapackSolver();
        LinearSOE *theFullLinSOE = new FullGenLinSOE(size, *theFullLinSolver);
        if (theFullLinSOE == 0)  {
            opserr << "WARNING CRExplicit::newStep() - failed to create FullLinearSOE\n";
            return -4;
        }
        theFullLinSOE->setLinks(*theModel);
        
        // now switch the SOE to the FullLinearSOE
        this->IncrementalIntegrator::setLinks(*theModel, *theFullLinSOE, theTest);
        
        // get a pointer to the A matrix of the FullLinearSOE
        const Matrix *tmp = theFullLinSOE->getA();
        if (tmp == 0)  {
            opserr << "WARNING CRExplicit::newStep() - ";
            opserr << "failed to get A matrix of FullGeneral LinearSOE\n";
            return -5;
        }
        
        // calculate the integration parameter matrices
        // c1:K  c2:C  c3:M
        c1 = deltaT*deltaT;
        c2 = 2.0*deltaT;
        c3 = 4.0;
        //this->TransientIntegrator::formTangent(INITIAL_TANGENT); bug fix by zxh
        this->TransientIntegrator::formTangent(CURRENT_TANGENT);
        Matrix A(*tmp);
        
        c1 = 0.0;
        c2 = 0.0;
        c3 = 4.0;
        //this->TransientIntegrator::formTangent(INITIAL_TANGENT); bug fix by zxh
        this->TransientIntegrator::formTangent(CURRENT_TANGENT);
        Matrix B1(*tmp);
        
        // solve [4M + 2*deltaT*C + deltaT^2*K]*[alpha1] = [4*M]
        // for alpha1
        A.Solve(B1, *alpha1);
        
        
        // calculate the effective mass matrix Mhat
        Mhat->addMatrix(0.0, B1, 0.25);
        
        // switch the SOE back to the user specified one
        this->IncrementalIntegrator::setLinks(*theModel, *theLinSOE, theTest);
        
        initAlphaMatrices = 0;
    }
    
    if (U == 0)  {
        opserr << "WARNING CRExplicit::newStep() - domainChange() failed or hasn't been called\n";
        return -6;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // determine new response at time t+deltaT
    // alpha1=alpha2
    Utdothat->addMatrixVector(0.0, *alpha1, *Utdotdot, deltaT);
    
    Udot->addVector(1.0, *Utdothat, 1.0);
    
    U->addVector(1.0, *Utdot, deltaT);
    U->addVector(1.0, *Utdothat, deltaT);

    
    // set the trial response quantities
    Udotdot->addVector(0, *Utdotdot, 0);
    theModel->setResponse(*U, *Udot, *Udotdot);
    //theModel->setDisp(*U);
    //theModel->setVel(*Udot);
    
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "WARNING CRExplicit::newStep() - failed to update the domain\n";
        return -7;
    }
    
    return 0;
}


int CRExplicit::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int CRExplicit::formTangent(int statFlag)
{
    statusFlag = statFlag;
    
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING CRExplicit::formTangent() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    theLinSOE->zeroA();
    
    int size = theLinSOE->getNumEqn();
    ID id(size);
    for (int i=1; i<size; i++)  {
        id(i) = id(i-1) + 1;
    }
    if (theLinSOE->addA(*Mhat, id) < 0)  {
        opserr << "WARNING CRExplicit::formTangent() - ";
        opserr << "failed to add Mhat to A\n";
        return -2;
    }
    
    return 0;
}


int CRExplicit::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)
        theEle->addKtToTang(c1);
    else if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(c1);
    
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
    
    return 0;
}


int CRExplicit::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int CRExplicit::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // create the new Matrix and Vector objects
    if (Ut == 0 || Ut->Size() != size)  {
        
        // delete the old
        if (alpha1 != 0)
            delete alpha1;
        if (Mhat != 0)
            delete Mhat;
        if (Ut != 0)
            delete Ut;
        if (Utdot != 0)
            delete Utdot;
        if (Utdotdot != 0)
            delete Utdotdot;
        if (U != 0)
            delete U;
        if (Udot != 0)
            delete Udot;
        if (Udotdot != 0)
            delete Udotdot;
        if (Utdothat != 0)
            delete Utdothat;
        
        // create the new
        alpha1 = new Matrix(size,size);
        Mhat = new Matrix(size,size);
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Utdothat = new Vector(size);
        
        // check we obtained the new
        if (alpha1 == 0 || alpha1->noRows() != size || alpha1->noCols() != size ||
            Mhat == 0 || Mhat->noRows() != size || Mhat->noCols() != size ||
            Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Utdothat == 0 || Utdothat->Size() != size)  {
            
            opserr << "WARNING CRExplicit::domainChanged() - ";
            opserr << "ran out of memory\n";
            
            // delete the old
            if (alpha1 != 0)
                delete alpha1;
            if (Mhat != 0)
                delete Mhat;
            if (Ut != 0)
                delete Ut;
            if (Utdot != 0)
                delete Utdot;
            if (Utdotdot != 0)
                delete Utdotdot;
            if (U != 0)
                delete U;
            if (Udot != 0)
                delete Udot;
            if (Udotdot != 0)
                delete Udotdot;
            if (Utdothat != 0)
                delete Utdothat;
            
            alpha1 = 0;; Mhat = 0;
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Utdothat = 0;
            
            return -1;
        }
    }
    
    // now go through and populate U, Udot and Udotdot by iterating through
    // the DOF_Groups and getting the last committed velocity and accel
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0)  {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        
        int i;
        const Vector &disp = dofPtr->getCommittedDisp();
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*U)(loc) = disp(i);
            }
        }
        
        const Vector &vel = dofPtr->getCommittedVel();
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Udot)(loc) = vel(i);
            }
        }
        
        const Vector &accel = dofPtr->getCommittedAccel();
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Udotdot)(loc) = accel(i);
            }
        }
    }
    
    // recalculate integration parameter matrices b/c domain changed
    initAlphaMatrices = 1;
    
    return 0;
}


int CRExplicit::update(const Vector &aiPlusOne)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING CRExplicit::update() - called more than once -";
        opserr << " CRExplicit integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING CRExplicit::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING CRExplicit::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check aiPlusOne is of correct size
    if (aiPlusOne.Size() != U->Size())  {
        opserr << "WARNING CRExplicit::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << aiPlusOne.Size() << endln;
        return -4;
    }
    
    //  determine the response at t+deltaT
    //U->addVector(1.0, aiPlusOne, c1);  // c1 = 0.0
    
    //Udot->addVector(1.0, aiPlusOne, c2);  // c2 = 0.0
    
    Udotdot->addVector(0.0, aiPlusOne, 1.0);
    
    // update the response at the DOFs
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "CRExplicit::update() - failed to update the domain\n";
        return -5;
    }
    // do not update displacements in elements only at nodes
    theModel->setDisp(*U);



    
    return 0;
}


int CRExplicit::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING CRExplicit::commit() - no AnalysisModel set\n";
        return -1;
    }
    
    // have set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    //time += 0.0*deltaT;
    theModel->setCurrentDomainTime(time);
    
    // update the displacements in the elements
    if (updElemDisp == true)
        theModel->updateDomain();
    
    return theModel->commitDomain();
}

const Vector &
CRExplicit::getVel()
{
  return *Udot;
}

int CRExplicit::sendSelf(int cTag, Channel &theChannel)
{

    Vector data(1);

    if (updElemDisp == false)
        data(0) = 0.0;
    else
        data(0) = 1.0;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING KRAlphaExplicit::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int CRExplicit::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(1);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "WARNING KRAlphaExplicit::recvSelf() - could not receive data\n";
        return -1;
    }

    if (data(0) == 0.0)
        updElemDisp = false;
    else
        updElemDisp = true;
    
    return 0;
}


void CRExplicit::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "CRExplicit - currentTime: " << currentTime << endln ;
        //s << "  alpha1=alpha2: " << alpha1 << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        if (updElemDisp)
            s << "  updateElemDisp: yes\n";
        else
            s << "  updateElemDisp: no\n";
    } else
        s << "CRExplicit - no associated AnalysisModel\n";
}
