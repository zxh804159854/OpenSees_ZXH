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

#ifndef CRExplicit_h
#define CRExplicit_h

// Developed: Cheng Chen (chc4@lehigh.edu)
// Implemented: Xiaohang Zhang (xiaohangzhang@tju.edu.cn) 
// Implemented: Ning Li (neallee@tju.edu.cn)
// Created: 01/24
// Revision: A
//
// Description: This file contains the class definition for CRExplicit.
// CRExplicit is an algorithmic class for performing a transient analysis
// using the explicit Chen-Ricles integration scheme based on the mapping rule of poles.
//
// Reference: Chen C, Ricles J M. 
// Development of Direct Integration Algorithms for Structural Dynamics Using Discrete Control Theory.
// Journal of Engineering Mechanics, 2008, 134(8): 676-683.
// doi:10.1061/(ASCE)0733-9399(2008)134:8(676)



#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;
class Matrix;

class CRExplicit : public TransientIntegrator
{
public:
    // constructors
    CRExplicit();
    CRExplicit(bool updElemDisp = false);
    
    // destructor
    ~CRExplicit();
    
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
