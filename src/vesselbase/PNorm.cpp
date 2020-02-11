/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "VesselRegister.h"
#include "FunctionVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class PNorm : public FunctionVessel {
private:
  double power;
  double shift;
  double scale;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit PNorm( const VesselOptions& da );
  std::string value_descriptor();
  double calcTransform( const double& val, double& dv ) const ;
  double finalTransform( const double& val, double& dv );
};

PLUMED_REGISTER_VESSEL(PNorm,"PNORM")

void PNorm::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords( keys );
  keys.add("compulsory","P","the value of P");
  keys.add("compulsory","SHIFT","the value by which values are shifted before taking the norm");
  keys.add("compulsory","SCALE","the value by which values (after shifting) are multiplied");
}

void PNorm::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","PNORM","calculate the PNorm value.");
  keys.addOutputComponent("pnorm","PNORM","the PNorm value.");

}

PNorm::PNorm( const VesselOptions& da ) :
  FunctionVessel(da)
{
  if( getAction()->isPeriodic() ) error("PNORM is not a meaningful option for periodic variables");
  parse("P",power);
  parse("SHIFT",shift);
  parse("SCALE",scale);
  if( diffweight ) error("can't calculate PNORM if weight is differentiable");
}

std::string PNorm::value_descriptor() {
  std::string str_power; Tools::convert( power, str_power );
  return "the PNORM value. P is equal to " + str_power;
}

double PNorm::calcTransform( const double& val, double& dv ) const {
  double temp = val-shift;
  temp *= scale;
  double f = pow(temp, power); dv=power*f*scale/temp; return f;
}

double PNorm::finalTransform( const double& val, double& dv ) {
  double dist=pow(val, 1.0/power);
  dv = dist/(val*power); return dist;
}

}
}

