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
#include "core/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR BONDDISTORTION
/*
Calculate the global bond distortion

A bond distortion is defined as the deviation (typically due to stretching) of a bond length from
a reference value. A list of these bonds is generated at the first simulation step by storing all
atom pairs separated by a distance shorter than a user-defined cutoff value.
More details can be found in \cite Tiwary2013 and \cite Bal2015.

To obtain a global bond distortion, i.e., a measure of what is the strongest bond distortion in the
whole system, a p-norm is used:

\f[
s = \left ( \sum_i \left ( \frac{r_i - r_i^{min}}{r_i^{min}} \right )^p \right )^{1/p}
\f]

When taking larger values of the exponent, the larger the contribution of the most strongly distorted
bond in the system to the CV value. This way, the CV is suited to describe the breaking of any bond in
the system, without having to specify or identify these in advance.

It is also possible to "reset" the CV, meaning that a new list of bonds is created if the CV value is
above a certain threshold for a specified number of steps. This option is primarily intended to be used
in combination with the CVHD option in \ref METAD.

\par Examples

The following example instructs plumed to calculate the bond distortion function of atoms in two groups.
Atoms 1-20 in GROUPA have a reference bond length of 0.15 nm, and all interatomic contacts shorter than
0.18 nm at the first step are considered as bonds, whereas these values are 0.07 nm and 0.1 nm, respectively,
for atoms 21-60 in GROUPB. By default, bonds between atoms in GROUPA and GROUPB will also be considered,
where the RMIN and RCUT values are the arithmetic mean of the GROUPA and GROUPB values (which, in this case,
would be an RMIN of 0.11 nm and and RCUT of 0.14 nm). If, however, you wish to use other values for these 
cross-terms (or ignore them altogether), you can use the RMINAB and RCUTAB keywords.
\verbatim
BONDDISTORTION GROUPA=1-20 GROUPB=21-60 P=8 RMINA=0.15 RCUTA=0.18 RMINB=0.07 RCUTB=0.1 LABEL=ch
PRINT ARG=ch STRIDE=1
\endverbatim
(See also \ref PRINT)

It is not mandatory to use two groups. For instance, if we're only interested in the atoms 1-20 of the
previous example, an input could look like this:
\verbatim
BONDDISTORTION GROUPA=1-20 P=8 RMINA=0.15 RCUTA=0.18 LABEL=cc
\endverbatim

Further expanding on the previous example, if the CV value remains above 0.4 for 1000 simulation steps, we
instruct plumed to reset and generate a new reference list of bonds. Please note that these 1000 steps
correspond to 1000 steps at which the CV is actually calculated: if no bias (such as the CVHD method) is
applied, the CV is only calculated when it is printed.
\verbatim
BONDDISTORTION GROUPA=1-20 P=8 RMINA=0.15 RCUTA=0.18 RESET_BONDS RESET_MAXDIST=0.4 RESET_TIME=1000 LABEL=cc
\endverbatim

*/
//+ENDPLUMEDOC
   
class GlobalDistortion : public Colvar {
  bool pbc, twogroups, reBuildReflist;
  bool reset_ref, useSwitch, use_cutoff;
  bool do_bondform;
  double ref_val, ref_min, ref_max;
  double reset_maxdist;
  unsigned power;
  unsigned reset_time, reset_wait;
  unsigned nl_stride, nl_wait;
  unsigned num_a, num_b;
  SwitchingFunction switchingFunction;
  vector<AtomNumber> ga_lista, gb_lista;
  std::vector<AtomNumber> fullatomlist;
  vector<double> reflist;
  vector<bool> checklist;
  void buildReflist();
  void buildNeigbourlist();
  double getRefval(double val);
  double pairterm(unsigned i, unsigned j, double ref, vector<Vector>& der, Tensor& viral);
public:
  explicit GlobalDistortion(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(GlobalDistortion,"GLOBALDISTORTION")

void GlobalDistortion::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.addFlag("RESET_REF",false,"Reset the computation by generating a new reference list, as specified by the RESET_MAXDIST and RESET_TIME keywords");
  keys.addFlag("USE_CUTOFF",false,"Use the REF_MIN and REF_MAX keywords to limited the number of contacts consider");
  keys.addFlag("DO_BONDFORM",false,"bias bond formation. Instead of specifying a reference bond length, define a reference nonbonded contact.");
  keys.add("compulsory","P", "8","The p parameter for computing the norm");
  keys.add("compulsory","REF","The reference value from which its multiples the distortions are calculated");
  keys.add("optional","REF_MIN","Don't consider contacts of which the interaction parameter is smaller than this value");
  keys.add("optional","REF_MAX","Don't consider contacts of which the interaction parameter is larger than this value");
  keys.add("optional","NL_STRIDE","Update the neigbour list every this many steps (based on REF_MIN and REF_MAX)");
  keys.add("optional","RESET_MAXDIST","Maximal distortion before the reference list is reset (if RESET_BONDS is enabled)");
  keys.add("optional","RESET_TIME","Number of steps to wait before resetting the reference list");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available."); 
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, only contacts within GROUPA are considered)");
}

GlobalDistortion::GlobalDistortion(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  vector<AtomNumber> ga_lista, ga_listb;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  num_a = ga_lista.size();
  num_b = num_a + gb_lista.size();
  fullatomlist=ga_lista;
  fullatomlist.insert(fullatomlist.end(),gb_lista.begin(),gb_lista.end());

  reBuildReflist = true;
  twogroups = true;
  reset_ref = false;
  use_cutoff = false;
  do_bondform = false;
  reset_time = 0;
  reset_wait = 0;
  nl_stride = 0;
  nl_wait = 0;
  if(gb_lista.size() == 0) twogroups = false;
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  parseFlag("RESET_REF",reset_ref);
  parse("RESET_MAXDIST",reset_maxdist);
  parse("RESET_TIME",reset_time);

  parse("REF",ref_val);
  parse("P",power);

  if(ref_val < 0.0) error("REF should be explicitly specified and positive");
  if(power < 0.0) error("P should be explicitly specified and positive");

  parseFlag("USE_CUTOFF",use_cutoff);
  parse("REF_MIN",ref_min);
  parse("REF_MAX",ref_max);

  parse("NL_STRIDE",nl_stride);

  parseFlag("DO_BONDFORM",do_bondform);

  useSwitch = false;
  string sw,errors;
  parse("SWITCH",sw);
  if(sw.length()>0){
    useSwitch = true;
    if( do_bondform ) error("it does not make sense to use DO_BONDFORM with SWITCH");
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  }

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0;i<gb_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  checkRead();
  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(fullatomlist);
}


// calculator
void GlobalDistortion::calculate(){
  if (reBuildReflist) buildReflist();

  int k = 0;
  double pairsum = 0.0;
  double prefactor;
  double value;
  Vector distance;
  vector<Vector> deriv(getNumberOfAtoms());
  Tensor virial;

  if(!twogroups) {
    for(unsigned i = 0; i < num_a; i++){
      for(unsigned j = i+1; j < num_a; j++){
        if(checklist[k]) pairsum += pairterm(i, j, reflist[k], deriv, virial);
        k++;
      }
    }
  } else {
    for(unsigned i = 0; i < num_a; i++){
      for(unsigned j = num_a; j < num_b; j++){
        if(checklist[k]) pairsum += pairterm(i, j, reflist[k], deriv, virial);
        k++;
      }
    }
  }

  if (pairsum != 0.0){
    value = pow(pairsum, 1.0/power);
    prefactor = pow(pairsum, 1.0/power-1.0);
  } else {
    value = 0.0;
    prefactor = 0.0;
  }

  // bookkeeping. if we've exceeded our distortion limit, prepare for rebuild next time
  unsigned currstep = getStep();
  if(reset_ref){
    if(value < reset_maxdist) reset_wait=currstep;
    if(currstep-reset_wait >= reset_time){
      reset_wait = currstep;
      reBuildReflist = true;
    }
  }

  // neigbour list stuff. should be used when using full connectivity matrix in a large system
  if(nl_stride>0 && currstep-nl_wait>=nl_stride){
    if(!(reset_ref && value>=reset_maxdist)){
      buildNeigbourlist();
      nl_wait=currstep;
    }
  }

  for(unsigned i=0;i<deriv.size();++i) setAtomsDerivatives(i,deriv[i]*prefactor);
  setBoxDerivatives  (virial*prefactor);
  setValue           (value);
}

void GlobalDistortion::buildReflist(){
  reBuildReflist = false;
  reflist.clear();
  checklist.clear();

  Vector distance;
  double val, dummy;
  bool consider;

  if(!twogroups) {
    for(unsigned i = 0; i < num_a; i++){
      for(unsigned j = i+1; j < num_a; j++){
        if(pbc){
          distance = pbcDistance(getPosition(i),getPosition(j));
        } else {
          distance = delta(getPosition(i),getPosition(j));
        }
        val = distance.modulo();
        if(useSwitch) val = switchingFunction.calculate(val,dummy);
        consider = true;
        if(use_cutoff && (val < ref_min || val > ref_max)) consider = false;
        reflist.push_back(getRefval(val));
        checklist.push_back(consider);
      }
    }
  } else {
    for(unsigned i = 0; i < num_a; i++){
      for(unsigned j = num_a; j < num_b; j++){
        if(pbc){
          distance = pbcDistance(getPosition(i),getPosition(j));
        } else {
          distance = delta(getPosition(i),getPosition(j));
        }
        val = distance.modulo();
        if(useSwitch) val = switchingFunction.calculate(val,dummy);
        consider = true;
        if(use_cutoff && (val < ref_min || val > ref_max)) consider = false;
        reflist.push_back(getRefval(val));
        checklist.push_back(consider);
      }
    }
  }
}

void GlobalDistortion::buildNeigbourlist(){
  checklist.clear();

  Vector distance;
  double val, dummy;
  bool consider;

  if(!twogroups) {
    for(unsigned i = 0; i < num_a; i++){
      for(unsigned j = i+1; j < num_a; j++){
        if(pbc){
          distance = pbcDistance(getPosition(i),getPosition(j));
        } else {
          distance = delta(getPosition(i),getPosition(j));
        }
        val = distance.modulo();
        if(useSwitch) val = switchingFunction.calculate(val,dummy);
        consider = true;
        if(use_cutoff && (val < ref_min || val > ref_max)) consider = false;
        checklist.push_back(consider);
      }
    }
  } else {
    for(unsigned i = 0; i < num_a; i++){
      for(unsigned j = num_a; j < num_b; j++){
        if(pbc){
          distance = pbcDistance(getPosition(i),getPosition(j));
        } else {
          distance = delta(getPosition(i),getPosition(j));
        }
        val = distance.modulo();
        if(useSwitch) val = switchingFunction.calculate(val,dummy);
        consider = true;
        if(use_cutoff && (val < ref_min || val > ref_max)) consider = false;
        checklist.push_back(consider);
      }
    }
  }
}

double GlobalDistortion::getRefval(double val){
  if (useSwitch) {
    return ref_val*round(val/ref_val);
  } else { return ref_val;
  }
}

double GlobalDistortion::pairterm(unsigned i, unsigned j, double ref, vector<Vector>& der, Tensor& virial){
  double r;
  double stretch;
  double val;
  double dcoord;
  Vector distance;
  Vector dval;

  if(pbc){
    distance=pbcDistance(getPosition(i),getPosition(j));
  } else {
    distance=delta(getPosition(i),getPosition(j));
  }

  r = distance.modulo();
  dcoord = 1.0/r;
  if(useSwitch) { 
    r = switchingFunction.calculate(r,dcoord);
  } else if ((do_bondform && r > ref) || (!do_bondform && r < ref)) {
    return 0.0;
  }
  stretch = (r-ref)/ref_val;
  val = pow(stretch, power);
  dval = dcoord*distance*pow(stretch, power-1.0)/ref_val;
  Tensor vv(dval,distance);

  der[i] -= dval;
  der[j] += dval;
  virial -= vv;
  return val;
}

}
}



