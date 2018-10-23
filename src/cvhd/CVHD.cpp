/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "function/Function.h"
#include "core/ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION CVHD
/*
Convert a (set of) CV(s) into the CV of the CVHD method.

In the simplest use case (one argument), the CV is calculated as follows for
values between 0 and CUTOFF:
\f[
C=\frac{1}{2} \left ( 1 - \cos left ( \pi \frac{x}{x_{cut}} \right ) \right )
\f]
Effectively, any CV is projected on an interval between 0 and 1.

If multiple arguments are supplied, these are first combined as follows, using
the power supplied through the POWER keyword.
\f[
x=\left ( \sum_{i} x_{i}^{p} \right )^{1/p}
\f]

The coefficients c, the parameters a and the powers p are provided as vectors.

Typically, this function is used with the BONDDISTORTION or GLOBALDISTORTION
CVs, as an argument for a metadynamics or VES bias in CVHD mode.

Notice that CVHD is not able to predict which will be periodic domain
of the computed value automatically. The user is thus forced to specify it
explicitly. Use PERIODIC=NO if the resulting variable is not periodic,
and PERIODIC=A,B where A and B are the two boundaries if the resulting variable
is periodic.

*/
//+ENDPLUMEDOC


class CVHD :
  public Function
{
  double cutoff;
  double power;
public:
  explicit CVHD(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(CVHD,"CVHD")

void CVHD::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG"); keys.use("PERIODIC");
  keys.add("compulsory","CUTOFF","1.0","the cutoff distance");
  keys.add("compulsory","POWER","1.0","the powers to which you are raising each of the arguments in your function");
}

CVHD::CVHD(const ActionOptions&ao):
Action(ao),
Function(ao)
{
  cutoff = 1.0;
  power = 1.0;
  parse("CUTOFF",cutoff);
  parse("POWER",power);

  addValueWithDerivatives();
  checkRead();

  log.printf("  with cutoff:");
  log.printf(" %f",cutoff);
  log.printf("\n");
  log.printf("  with power:");
  log.printf(" %f",power);
  log.printf("\n");
}

void CVHD::calculate(){
  const double pi = 3.141592653589793;
  const double cutoff2 = cutoff*cutoff;
  double combine = 0.0;
  double val, prefact, arg, pnorm;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    double cv = getArgument(i);
    combine += pow(cv, power);
  };
  pnorm = pow(combine,2.0/power);
  arg = pnorm/cutoff2;
  if(arg < 1.0 && arg > 0.0){
    val = 0.5*(1.0-cos(pi*arg));
    prefact = pow(combine, 2.0/power-1.0)*pi*sin(pi*arg)/cutoff2;
  } else if (arg >= 1.0) {
    val = 1.0;
    prefact = 0.0;
  } else {
    val = 0.0;
    prefact = 0.0;
  }
  for(unsigned i=0;i<getNumberOfArguments();++i){
    double cv = getArgument(i);
    setDerivative(i,prefact*pow(cv,power-1.0));
  };
  setValue(val);
}

}
}
