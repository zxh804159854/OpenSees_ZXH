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

// $Revision$
// $Date$
// $URL$

#ifndef GuiLambdaExplicit_h
#define GuiLambdaExplicit_h

// Developed: Gui Y, Wang J T, Jin F, et al. (jinfeng@tsinghua.edu.cn)
// Implemented: Xiaohang Zhang (xiaohangzhang@tju.edu.cn) 
// Implemented: Ning Li (neallee@tju.edu.cn)
// Created: 01/24
// Revision: A
//
// Description: This file contains the class definition for GuiLambdaExplicit.
// GuiLambdaExplicit is an algorithmic class for performing a transient analysis
//
// Reference: Gui Y, Wang J T, Jin F, et al.
// Development of a family of explicit algorithms for structural dynamics with unconditional stability.
// Nonlinear Dynamics, 2014, 77(4) : 1157 - 1170.




#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;
class Matrix;

class GuiLambdaExplicit : public TransientIntegrator
{
public:
    // constructors
    GuiLambdaExplicit();
    GuiLambdaExplicit(double lambda ,
        bool updElemDisp = false);
    
    // destructor
    ~GuiLambdaExplicit();
    
    // method to set up the system of equations
    int formTangent(int statFlag);
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    // methods to update the domain
    int domainChanged(void);
    int newStep(double deltaT);
    int revertToLastStep(void);
    int update(const Vector &aiPlusOne);
    int commit(void);

    const Vector &getVel(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:
    
private:
    double lambda;
    bool updElemDisp;  // a flag indicating if element displacements are updated during commit
    double deltaT;
    
    Matrix *alpha1;  // integration parameter matrices, alpha2 = *alpha1
    Matrix *Mhat;             // effective mass matrix for linear SOE
    
    int updateCount;                            // method should only have one update per step
    int initAlphaMatrices;                      // a flag to initialize the alpha matrices
    double c1, c2, c3;                          // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;              // response quantities at time t
    Vector *U, *Udot, *Udotdot;                 // response quantities at time t + deltaT
    Vector *Utdothat;                           // velocity-like vector at time t
};

#endif
