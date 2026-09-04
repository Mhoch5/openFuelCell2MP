/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::sigmaModels::resistanceSigma

Description
    Electrical conductivity based on experimentally derived electrical resistance.
    The conductivity is calculated as follows.
    \f[
        \sigma = thickness / resistance
    \f]
    
    where
    \vartable
        thickness           | Thickness of the component, e.g. GDL [m]
        resistance          | Electrical resistance of the component [Omega m^2]
    \endvartable

SourceFiles
    resistanceSigma.C

\*---------------------------------------------------------------------------*/

#ifndef resistanceSigma_H
#define resistanceSigma_H

#include "sigmaModel.H"
#include "List.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaModels
{

/*---------------------------------------------------------------------------*\
                           Class resistanceSigma Declaration
\*---------------------------------------------------------------------------*/

class resistanceSigma
:
    public sigmaModel
{
    // Private data

        //- thickness
        scalar thickness_;

        //- resistance
        scalar resistance_;

        //- catalyst
        bool catalyst_;

public:

    //- Runtime type information
    TypeName("resistanceSigma");


    // Constructors

        //- Construct from components
        resistanceSigma
        (
            const fvMesh& mesh,
            const dictionary& sigmaDictionary
        );


    //- Destructor
    virtual ~resistanceSigma();


    // Member Functions

        //- Return the sigma field
        virtual void correct(volScalarField& sigmaField) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sigmaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
