/*
 * GiNaCRA - GiNaC Real Algebra package
 * Copyright (C) 2010-2012  Ulrich Loup, Joachim Redies, Sebastian Junges
 *
 * This file is part of GiNaCRA.
 *
 * GiNaCRA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GiNaCRA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GiNaCRA.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


/**
 * Example file for testing the GiNaC Real Algebra package.
 *
 * Compile with:
 * g++ -lcln -lginac -lginacra -o example_IR_0 example_IR_0.cpp
 * @author Ulrich Loup
 * @since 2010-12-18
 * @version 2010-12-25
 */

#include <ctime>
#include <iostream>

using namespace std;

#include <ginacra/ginacra.h>

using namespace GiNaC;
using namespace GiNaCRA;

const string   VERSION        = "2010-12-25";
const string   SUPPORT        = "Ulrich Loup <loup* @cs.rwth-aachen.de>";
const unsigned TESTCOUNT      = 100;    // number of repetitions
const unsigned TERMCOUNT      = 62;    // number of terms in the polynomial except for the constant part
const unsigned MAXDEGREE      = 24;    // maximum degree (+1) of polynomials
const unsigned MAXCOEFFICIENT = 1001;    // maximum coefficient (+1) of polynomials

//////////////////
// Main program //
//////////////////

int main( int argc, char** argv )
{
    cout << endl;
    cout << " GiNaCRA Example on RealAlgebraicNumberIR Version " << VERSION << endl << "";
    cout << " Support: " << SUPPORT << endl << endl;

    GiNaC::symbol x( "x" );    // main variable of the GiNaCRA polynomial
    cout << "Starting " << TESTCOUNT << " real root computations of a template polynomial with random degrees, random coefficients and "
         << (TERMCOUNT + 1) << " terms." << endl << endl;
    for( unsigned i = 0; i < TESTCOUNT; ++i )
    {
        // construct random polynomial
        const int constant = rand() % MAXCOEFFICIENT;
        ex pEx( constant );    // constant part
        for( unsigned t = 0; t < TERMCOUNT; ++t )
        {
            int coefficient = rand() % MAXCOEFFICIENT;
            int exponent    = rand() % MAXDEGREE;
            pEx += coefficient * pow( x, GiNaC::abs( exponent ));
        }
        RationalUnivariatePolynomial pGiNaCRA( pEx, x );

        cout << endl << "------------------------------- Test " << i << " -------------------------------" << endl;
        cout << "**** GiNaCRA ****" << endl;
        cout << "GiNaCRA Example polynomial: " << pGiNaCRA << endl;
        // compute real roots
        double                       startTime = clock();
        list<RealAlgebraicNumberPtr> roots     = RealAlgebraicNumberFactory::realRoots( pGiNaCRA );
        // access the list of real roots
        cout << "Has " << roots.size() << " real root(s):" << endl;
        for( list<RealAlgebraicNumberPtr>::const_iterator root = roots.begin(); root != roots.end(); ++root )
        {
            cout << " " << *root << endl;
            // access individual roots
            cout << " Verified root of the polynomial: " << ((*root)->sgn( pGiNaCRA ) == 0 ? "True" : "False") << endl;
        }
        cout << "**** Test " << i << " took " << (clock() - startTime) / CLOCKS_PER_SEC << " sec with GiNaCRA" << endl;
    }
    cout << "Example finished." << endl << endl;
    return 0;
}
