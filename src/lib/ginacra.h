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
 * @mainpage GiNaCRA - GiNaC Real Algebra package
 *
 * Welcome to the GiNaC Real Algebra project! In the following we will give some short abstracts on what to expect from this package.
 * These abstracts are excerpts from published papers on GiNaCRA.
 *
 * @section intro Introduction at NFM in 2011
 *
 * GiNaCRA provides efficient and easy-to-integrate data structures and methods for real algebra.
 * It is based on the C++ library GiNaC, supporting the symbolic representation and manipulation of polynomials. In contrast to other
 * similar tools, our open source library aids exact, real algebraic computations based on an appropriate data type representing real zeros of polynomials.
 * The only non-standard library GiNaCRA depends on is GiNaC, which makes the installation and usage of our library simple. Our longterm goal is
 * to integrate decision procedures for real algebra within the Satisfiability-Modulo-Theories (SMT) context and thereby provide tool support for many
 * applied formal methods.
 *
 * @see http://www.springerlink.com/content/x444nt72x4236073/
 */

#ifndef GINACRA_H
#define GINACRA_H

/**
 * Collection containing all classes and useful algorithms for real algebra.
 *
 * @author Ulrich Loup
 * @since 2010-08-10
 * @version 2011-12-06
 * @see ISBN 0-387-94090-1 and ISBN-13: 978-3642069642
 *
 * Notation is following http://www.possibility.com/Cpp/CppCodingStandard.html.
 */

#include "settings.h"
#include "constants.h"
#include "utilities.h"
#include "Polynomial.h"
#include "Constraint.h"
#include "UnivariatePolynomial.h"
#include "RationalUnivariatePolynomial.h"
#include "OpenInterval.h"
#include "MultivariateMonomialMR.h"
#include "MultivariateTermMR.h"
#include "MultivariatePolynomialMR.h"
#include "Groebner.h"
//#include "MultivariatePolynomialFactory.h"
//#include "SpecialQuotientRingMultiplicationTable.h"

#include "RealAlgebraicNumber.h"
#include "RealAlgebraicNumberNR.h"
#include "RealAlgebraicNumberIR.h"
#include "RealAlgebraicNumberFactory.h"
#include "UnivariatePolynomialSet.h"
#include "tree.h"
#include "RealAlgebraicPoint.h"
#include "CAD.h"
#include "operators.h"

#endif // GINACRA_H





