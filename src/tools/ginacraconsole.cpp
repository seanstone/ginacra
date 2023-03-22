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
 * Testing shell for the GiNaC Real Algebra package.
 * @author Ulrich Loup
 * @author Bateer Siqin
 * @since 2010-09-10
 * @version 2012-02-02
 */

#include "ginacraconsole.h"

using namespace GiNaC;

///////////////////////
// Auxiliary methods //
///////////////////////

static char** auto_completion( const char* text, int start, int end )
{
    char** matches;
    matches = (char**)NULL;
    if( start == 0 )
        matches = rl_completion_matches( (char*)text, &auto_generator );
    else
        rl_bind_key( '\t', rl_insert );    // rl_abort is illegal and compiler-dependent
    return (matches);
}

char* auto_generator( const char* text, int state )
{
    static unsigned list_index, len;
    string          name;
    if( !state )
    {
        list_index = 0;
        len        = strlen( text );
    }
    while( list_index < COMMAND_COUNT )
    {
        name = COMMANDS[list_index++];
        if( strncmp( name.c_str(), text, len ) == 0 )
            return (dupstr( name.c_str() ));
    }
    // If no names matched, then return NULL.
    return ((char*)NULL);
}

char* dupstr( const char* s )
{
    char* r;
    r = (char*)xmalloc( (strlen( s ) + 1) );
    strcpy( r, s );
    return (r);
}

void* xmalloc( int size )
{
    void* buf;
    buf = malloc( size );
    if( !buf )
    {
        fprintf( stderr, "Error: Out of memory. Exiting.'n" );
        exit( 1 );
    }
    return buf;
}

///////////////////////////////////////////////////////////////////////////

void printHelp( string context ) throw ( invalid_argument )
{
    if( context.empty() )
    {
        cout << endl;
        cout << "   h[elp]                                  Display this message." << endl;
        cout << "   h[elp] [command]                        Help about the specific command." << endl << endl;
        cout << "   q[uit], exit                            Abort the interpreter." << endl << endl;
        cout << "   warranty                                Display Disclaimer of Warranty." << endl;
        cout << "   conditions                              Display Terms and Conditions." << endl << endl;
        cout << "   realRoots arg1 arg2 [arg3]              Compute the real roots of polynomial arg1 w.r.t. the variable arg2 [up to a precision of arg3]."
             << endl << endl;
        cout << "   cad arg1 arg2 arg3                      Compute a sample satisfying the given sign condition (arg1, arg2). arg1 specifies the polynomials"
             << endl
             << "                                           (comma-separated list of polynomials without white spaces) w.r.t. the variables arg3 (comma-separated"
             << endl
             << "                                           list of symbols without white spaces). arg2 specifies the sign condition (comma-separated list of numbers -1,0,1)."
             << endl;
    }
    else if( context.compare( CMD_H ) == 0 || context.compare( CMD_HELP ) == 0 || context.compare( CMD_Q ) == 0 || context.compare( CMD_QUIT ) == 0
             || context.compare( CMD_EXIT ) == 0 || context.compare( CMD_WARRANTY ) == 0 || context.compare( CMD_CONDITIONS ) == 0 )
        cout << "No arguments needed for the " << context << " command." << endl;
    else if( context.compare( CMD_REALROOTSUNIVARIATE ) == 0 )
    {
        cout << "The format of input should be \"realRoots arg1 arg2 [arg3]\"" << endl;
        cout << "Here arg1 represents the polynomial, arg2 represents the variable with which you want to compute the real roots, arg3 represents the width of the interval of the results, which is optional.";
        cout << endl << "Example: " << endl << PROMPT << "realRoots x^5-3*x^4*x^3-x^2+2*x-2 x" << endl;
    }
    else if( context.compare( CMD_CAD ) == 0 )
    {
        cout << "The format of input should be \"cad arg1 arg2 arg3\"" << endl;
        cout << "arg1 specifies the input polynomials (comma-separated list of polynomials without white spaces) w.r.t. the variables arg3 (comma-separated"
             << "list of symbols without white spaces). arg2 specifies the sign condition (comma-separated list of numbers -1,0,1).";
        cout << endl << "Example: " << endl << PROMPT << "cad (x-3)^2+(y-3)^2-r1,(x+3)^2+(y+3)^2-r2,r1-3,r2-3 -1,-1,1,1 r1,r2,x,y" << endl;
    }
    else
        throw invalid_argument( "Command \"" + context + "\" not found" );
    cout << endl;
    context = "";
}

void printError( string e )
{
    cout << endl;
    cout << "ERROR: " << e << "." << endl << endl;
}

//////////////////
// Test program //
//////////////////

int main( int argc, char** argv )
{
    cout << endl;
    cout << PACKAGE_NAME << " Version " << VERSION << endl << "";
    cout << "      Command Line Interpreter " << endl << endl;
    cout << WELCOME << endl << endl;
    cout << COPYRIGHT << endl << endl;
    cout << "Support: " << SUPPORT << endl << endl;
    cout << "Usage: Type a command and press enter." << endl;
    cout << "       Use the command \"help\" for more information." << endl << endl;
    string      cmd, arg1, arg2, arg3;
    const char* c;
    char*       buf;
#ifdef HAVE_READLINE
    rl_attempted_completion_function = auto_completion;
#endif

    while( (buf = readline( PROMPT.c_str() )) != NULL )
    {
        if( strlen( buf ) == 0 )
            continue;
        int pos = 0;
        int index;
        string str = "";
        cmd        = "";
        arg1       = "";
        arg2       = "";
        arg3       = "";
        string::size_type opos;

        //auto-complete
        rl_bind_key( '\t', rl_complete );
        //read input: cmd arg1 arg2 ...

        str  = buf;
        opos = str.find( " ", 0 );
        if( opos != string::npos )
        {
            index = str.find( " " );
            cmd   = str.substr( pos, index );
            pos   += cmd.size() + 1;
            str   = str.substr( pos );
            pos   = 0;
            opos  = str.find( " ", 0 );
            if( opos != string::npos )
            {
                index = str.find( " " );
                arg1  = str.substr( pos, index );
                pos   += arg1.size() + 1;
                str   = str.substr( pos );
                pos   = 0;
                opos  = str.find( " ", 0 );
                if( opos != string::npos )
                {
                    index = str.find( " " );
                    arg2  = str.substr( pos, index );
                    pos   += arg2.size() + 1;
                    str   = str.substr( pos );
                    arg3  = str;
                }
                else
                    arg2 = str;
            }
            else
                arg1 = str;
        }
        else
            cmd = str;
        try
        {
#ifdef HAVE_READLINE
            add_history( buf );
#endif
            if( cmd.compare( CMD_EXIT ) == 0 || cmd.compare( CMD_QUIT ) == 0 || cmd.compare( CMD_Q ) == 0 )
            {
                free( buf );
                return 0;
            }
            else if( cmd.compare( CMD_HELP ) == 0 || cmd.compare( CMD_H ) == 0 )
            {
                printHelp( arg1 );    // help with possible context
            }
            else if( cmd.compare( CMD_WARRANTY ) == 0 )
            {
                cout << WARRANTY << endl;
            }
            else if( cmd.compare( CMD_CONDITIONS ) == 0 )
            {
                cout << CONDITIONS << endl;
            }
            else if( cmd.compare( CMD_REALROOTSUNIVARIATE ) == 0 )
            {
                if( arg1.empty() )
                    throw invalid_argument( "No polynomial specified" );
                if( arg2.empty() )
                    throw invalid_argument( "No variable specified" );
                // in this case, cmd should look like "realRootsUnivariate x^2+1 x"
                if( arg3.empty() )
                {
                    cout << "Computing real roots of " << arg1 << " w.r.t. the variable " << arg2 << "..." << endl;
                }
                else
                {
                    cout << "Computing real roots of " << arg1 << " w.r.t. the variable " << arg2 << " up to a precision of " << arg3 << "..."
                         << endl;
                }
                parser reader;
                ex e         = reader( arg1 );
                symtab table = reader.get_syms();
                symtab::const_iterator sIt = table.find( arg2 );
                if( sIt == table.end() )
                    throw invalid_argument( "Specified variable not used in the polynomial" );
                symbol s = ex_to<symbol>( sIt->second );
                RationalUnivariatePolynomial p( e, s );
                list<RealAlgebraicNumberPtr> l = RealAlgebraicNumberFactory::realRoots( p );
                if( !l.empty() )
                {
                    if( arg3.empty() )    // no refinement
                    {
                        for( list<RealAlgebraicNumberPtr>::const_iterator i = l.begin(); i != l.end(); ++i )
                            print( *i, cout << "  " );
                        cout << endl;
                    }
                    else    // refinement specified
                    {
                        c = arg3.c_str();
                        for( list<RealAlgebraicNumberPtr>::iterator i = l.begin(); i != l.end(); ++i )
                        {
                            RealAlgebraicNumberIRPtr iIR = std::tr1::dynamic_pointer_cast<RealAlgebraicNumberIR>( *i );
                            if( iIR != 0 )
                                iIR->refine( atof( c ));
                            print( *i, cout << "  " );
                        }
                        cout << endl;
                    }
                }
                else
                    cout << "The polynomial has no real roots." << endl;
            }
            else if( cmd.compare( CMD_CAD ) == 0 )
            {
                if( arg1.empty() )
                    throw invalid_argument( "No polynomials specified" );
                if( arg2.empty() )
                    throw invalid_argument( "No variables specified" );
                // in this case, cmd should look like "realRootsUnivariate x^2+1 x"
                if( arg3.empty() )
                    throw invalid_argument( "No signs specified" );
                cout << "Computing a CAD for the sign condition (" << arg1 << " | " << arg2 << ") w.r.t. the variables " << arg3 << "..." << endl;
                parser reader;
                // parse variables
                vector<symbol> variables = vector<symbol>();
                size_t pos    = 0;
                size_t posOld = 0;
                while( true )
                {
                    pos = arg3.find( ",", posOld );
                    variables.push_back( ex_to<symbol>( reader( arg3.substr( posOld, pos - posOld ))));
                    if( pos == string::npos )
                        break;
                    ++pos;
                    posOld = pos;
                }
                // parse polynomials
                list<ex> polynomialExs = list<ex>();
                pos    = 0;
                posOld = 0;
                while( true )
                {
                    pos = arg1.find( ",", posOld );
                    polynomialExs.push_back( reader( arg1.substr( posOld, pos - posOld )));
                    if( pos == string::npos )
                        break;
                    ++pos;
                    posOld = pos;
                }
                // parse sign conditions
                list<int> signs = list<int>();
                pos    = 0;
                posOld = 0;
                while( true )
                {
                    pos = arg2.find( ",", posOld );
                    signs.push_back( ex_to<numeric>( reader( arg2.substr( posOld, pos - posOld ))).to_int() );
                    if( pos == string::npos )
                        break;
                    ++pos;
                    posOld = pos;
                }
                UnivariatePolynomialSet polynomialSet = UnivariatePolynomialSet();
                vector<Constraint>        constraints = vector<Constraint>();
                list<int>::const_iterator s           = signs.begin();
                for( list<ex>::const_iterator p = polynomialExs.begin(); s != signs.end() && p != polynomialExs.end(); ++p )
                {
                    polynomialSet.insert( UnivariatePolynomial( *p, variables.front() ));
                    sign sg = *s > 0 ? POSITIVE_SIGN : *s < 0 ? NEGATIVE_SIGN : ZERO_SIGN;
                    constraints.push_back( Constraint( Polynomial( *p ), sg, variables ));
                    ++s;
                }
                cout << "Solving the system " << endl;
                for( vector<Constraint>::const_iterator i = constraints.begin(); i != constraints.end(); ++i )
                    cout << "  " << *i << endl;
                cout << "by computing a CAD..." << endl;
                cout << "Initializing (elimination processing)...";
                CAD cad = CAD( polynomialSet, variables, UnivariatePolynomial::univariatePolynomialIsLessDeg );
                cout << " done." << endl;
                cout << "Elimination sets: " << endl;
                vector<vector<UnivariatePolynomial> > elimSets = cad.eliminationSets();
                for( unsigned i = 0; i != elimSets.size(); ++i )
                {
                    cout << "  Level " << i << ": ";
                    for( vector<UnivariatePolynomial>::const_iterator j = elimSets[i].begin(); j != elimSets[i].end(); ++j )
                        cout << *j << "   ";
                    cout << endl;
                }
                RealAlgebraicPoint r = RealAlgebraicPoint();
                cout << "Checking sign condition..." << endl;
                cout << "Result: " << cad.check( constraints, r ) << endl;
                cout << "Sample: " << r << endl;
            }
            else
                throw 1;
        }
        catch( invalid_argument e )
        {
            printError( e.what() );
        }
        catch( int e )
        {
            if( e == 1 )
                printError( "Command \"" + cmd + "\" not found" );
            else
                printError( "Error in (sub)command \"" + cmd + "\"." );
        }
    }
    free( buf );
    return 0;
}
