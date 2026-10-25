////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/atomdb/NCAtomDB.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace AtomDB {

    //We keep the database in two internal layers. One is essentially a static
    //and memory efficient list containing all the raw numbers for elements and
    //isotopes known to NCrystal, and the other is a standard CachedFactoryBase
    //which delivers the actual OptionalAtomDataSP objects, with on-demand
    //initialisation and cache-clearance as usual.

    namespace internal {

      using ScatLenFM = AtomData::ScatLenFM;

      class AtomDBKey {
      public:
        //Key is essentially (Z,A) where A=0 indicates natural element. For
        //efficiency this is packed into a single integer (properties similar to
        //std::pair<Z,A> but more efficient).
        AtomDBKey(unsigned Z, unsigned A) : m_key((Z<<16)+A) { nc_assert(isZAValid(Z,A)); }
        static bool isZAValid(unsigned Z, unsigned A) { return Z>0 && Z<150 && (A==0||A>=Z) && A<10000; }
        unsigned Z() const { return m_key>>16; }
        unsigned A() const { return m_key & 0xFFFF; }
        bool operator<( AtomDBKey o) const { return m_key < o.m_key; }
        bool operator==( AtomDBKey o) const { return m_key == o.m_key; }
        bool operator!=( AtomDBKey o) const { return m_key != o.m_key; }
      private:
        uint32_t m_key;
      };

      //Entry in internal list of raw numbers (MoveOnly to make sure that vector<Entry> is newer copied).
      struct Entry : private MoveOnly {
        unsigned Z() const { return m_key.Z(); }
        unsigned A() const { return m_key.A(); }
        AtomDBKey key() const { return m_key; }
        Entry(unsigned Z, unsigned A, double mm, ScatLenFM csl, double ixs, double axs)
          : m_key(Z,A), m_mass_amu(mm), m_coh_sl(csl), m_inc_xs(ixs), m_abs_xs(axs) {}
        bool operator<(const Entry&o) const {
          nc_assert( m_key!=o.m_key || this==&o );//check uniqueness of (Z,A) values during sort!
          return m_key < o.m_key;
        }
        bool operator<( AtomDBKey cmpkey ) const { return m_key < cmpkey; }

        AtomDataSP createAtomDataSP() const
        {
          return std::make_shared<AtomData>( SigmaBound{m_inc_xs},
                                             m_coh_sl,
                                             SigmaAbsorption{m_abs_xs},
                                             AtomMass{m_mass_amu},
                                             this->Z(),
                                             this->A() );
        }
        void modifyIncXS(double newval) { m_inc_xs = newval; }//for fixing Xenon
        void modifyCohScatLen(ScatLenFM newval) { m_coh_sl = newval; }//for fixing Calcium
        std::string getAtomDBLine() const {
          std::ostringstream ss;
          ss << elementZToName(Z());
          if ( A()!=0 )
            ss << A();
          ss << ' ' << fmt(m_mass_amu) << "u " << fmt(m_coh_sl.value) << "fm "
             << fmt(m_inc_xs) << "b " << fmt(m_abs_xs) << 'b';
          return ss.str();
        }
      private:
        AtomDBKey m_key;
        double m_mass_amu;
        ScatLenFM m_coh_sl;
        double m_inc_xs, m_abs_xs;
      };

      std::vector<Entry> setupDBValues();//Only call from inside internalDB()
      const std::vector<Entry>& internalDB() {
        //Wrapping internal DB in this function, to avoid static initialision
        //order fiasco and MT issues):
        static std::vector<Entry> s_db = setupDBValues();
        return s_db;
      }

      const Entry* lookupEntry( AtomDBKey key )
      {
        const auto& db = internalDB();
        auto it = std::lower_bound(db.begin(),db.end(),key);
        if ( it == db.end() || it->key() != key )
          return nullptr;
        return &(*it);
      }

      //This template argument of base class is true, i.e. factory keeps strong
      //not weak references.
      class StdAtomDataFactory : public NC::CachedFactoryBase<AtomDBKey,AtomData,CachedFactory_KeepAllStrongRefs> {
      public:
        const char* factoryName() const final { return "StdAtomDataFactory"; }
        std::string keyToString( const AtomDBKey& key ) const final
        {
          const unsigned Z = key.Z();
          const unsigned A = key.A();
          std::ostringstream ss;
          ss<<"(Z="<<Z;
          if (A==0) {
            ss<<";natural)";
          } else {
            ss<<";A="<<A<<")";
          }
          return ss.str();
        }

      protected:
        virtual OptionalAtomDataSP actualCreate( const AtomDBKey& key ) const final
        {
          auto entry = lookupEntry(key);
          if ( !entry )
            return nullptr;
          nc_assert( key.Z() == entry->Z() && key.A() == entry->A() );
          return entry->createAtomDataSP();
        }
      };

      static StdAtomDataFactory& getStdAtomDBFact()
      {
        static StdAtomDataFactory s_stdAtomDBFact;
        return s_stdAtomDBFact;
      }

    }
  }
}

NC::OptionalAtomDataSP NC::AtomDB::getNaturalElement( unsigned Z )
{
  if (!internal::AtomDBKey::isZAValid(Z,0))
    return nullptr;
  return internal::getStdAtomDBFact().create(internal::AtomDBKey(Z,0));
}

NC::OptionalAtomDataSP NC::AtomDB::getNaturalElement( const std::string& name )
{
  //Could use AtomSymbol class, but too simple:
  unsigned Z = elementNameToZ(name);
  if (Z==0)
    return nullptr;
  return internal::getStdAtomDBFact().create(internal::AtomDBKey(Z,0));
}

NC::OptionalAtomDataSP NC::AtomDB::getIsotope( unsigned Z, unsigned A )
{
  if (!internal::AtomDBKey::isZAValid(Z,A))
    return nullptr;
  return A>=Z ? internal::getStdAtomDBFact().create(internal::AtomDBKey(Z,A)) : nullptr;
}

NC::OptionalAtomDataSP NC::AtomDB::getIsotope( const std::string& name )
{
  AtomSymbol atomsymbol(name);
  if ( !atomsymbol.isIsotope() )
    return nullptr;
  return getIsotope(atomsymbol.Z(),atomsymbol.A());
}

NC::OptionalAtomDataSP NC::AtomDB::getIsotopeOrNatElem( unsigned Z, unsigned A )
{
  if (!internal::AtomDBKey::isZAValid(Z,A))
    return nullptr;
  return internal::getStdAtomDBFact().create(internal::AtomDBKey(Z,A));
}

NC::OptionalAtomDataSP NC::AtomDB::getIsotopeOrNatElem( const std::string& name )
{
  AtomSymbol atomsymbol(name);
  if ( ! (atomsymbol.isIsotope() || atomsymbol.isElement() ) )
    return nullptr;
  unsigned Z = atomsymbol.Z();
  unsigned A = atomsymbol.A();
  if (!internal::AtomDBKey::isZAValid(Z,A))
    return nullptr;

  return internal::getStdAtomDBFact().create(internal::AtomDBKey(Z,A));
}

unsigned NC::AtomDB::getAllEntriesCount()
{
  return internal::internalDB().size();
}

std::vector<std::pair<unsigned,unsigned>> NC::AtomDB::getAllEntries()
{
  const auto& db = internal::internalDB();
  std::vector<std::pair<unsigned,unsigned>> result;
  result.reserve(db.size());
  for (const auto& e : db )
    result.emplace_back(e.Z(),e.A());
  return result;
}

std::vector<NC::AtomDB::internal::Entry> NC::AtomDB::internal::setupDBValues()
{
  constexpr unsigned numberOfEntriesCountInDB = 348;//Exact number of isotopes and elements added below.
  std::vector<Entry> result;
  result.reserve(numberOfEntriesCountInDB);

  auto addNatElem = [&result](unsigned Z, double mass_amu, ScatLenFM coh_sl, double inc_xs, double abs_xs)
  {
    result.emplace_back(Z,0,mass_amu,coh_sl,inc_xs,abs_xs);
  };

  auto addIsotope = [&result](unsigned Z, unsigned A, double mass_amu, ScatLenFM coh_sl, double inc_xs, double abs_xs)
  {
    result.emplace_back(Z,A,mass_amu,coh_sl,inc_xs,abs_xs);
  };

  //Note that in addition to elements/isotopes added inside the next piece of
  //auto-generated code, we perform a few manual additions/fixes below the
  //auto-generated code, so don't forget to check the bottom of the file as
  //well!!

  // ----- Autogenerated code begin (no manual edits inside please!!!) ----- //

  /////////////////////////////////////////////////////////////////////////////////////////////
  //                                                                                         //
  // Neutron scattering lengths and cross sections are taken from recommended                //
  // values in [1] according to [2,3].                                                       //
  // [1] "Table of coherent scattering lengths and cross sections", available at:            //
  //      http://www.ati.ac.at/~neutropt/scattering/Scattering_lengths_table_20010419.pdf    //
  //      NB: The values for this file were parsed from the HTML table at (some values might //
  //      differ from the PDF above, but hopefully nothing significant...):                  //
  //      http://www.ati.ac.at/~neutropt/scattering&cross/data.html                          //
  // [2] LANDOLT-BORNSTEIN, New Series I/16A (Ed. H. Schopper) Chap.6, Springer, Berlin 2000 //
  // [3] Atomic Data Nuclear Data Tables 49 (1991) 65                                        //
  //                                                                                         //
  // Isotope masses and natural abundances are from:                                         //
  //                                                                                         //
  //       "The Ame2016 atomic mass evaluation (I)"                                          //
  //       W.J.Huang, et.al., Chinese Physics C41 030002, March 2017.                        //
  //       doi:10.1088/1674-1137/41/3/030003.                                                //
  //       http://amdc.impcas.ac.cn/masstables/Ame2016/mass16.txt                            //
  //                                                                                         //
  // Masses of natural elements are for consistency (and simplicity) calculated directly     //
  // from isotope masses and abundances. This can of course give slight deviations from      //
  // official CIAAW values, with the uncertainty almost exclusively attributed to the        //
  // abundance values. Most users of NCrystal will likely not care about these small         //
  // differences, and those that do should probably provide their own specific values        //
  // directly (using the @ATOMDB section in NCMAT files, or the atomdb cfg parameter).       //
  //                                                                                         //
  // Where available natural abundances of isotopes have been added (commented out           //
  // for now, but /* */ chars can be removed if we ever need them).                          //
  //                                                                                         //
  /////////////////////////////////////////////////////////////////////////////////////////////

  addNatElem(1, 1.007975973752, ScatLenFM{-3.739}, 80.26, 0.3326);//H
  addIsotope(1, 1, 1.00782503224, ScatLenFM{-3.7406}, 80.27, 0.3326 /*,0.99985*/);
  addIsotope(1, 2, 2.01410177811, ScatLenFM{6.671}, 2.05, 0.000519 /*,0.00015*/);
  addIsotope(1, 3, 3.01604928199, ScatLenFM{4.792}, 0.14, 0);
  addNatElem(2, 4.0026018729265, ScatLenFM{3.26548}, 0, 0.00747);//He
  addIsotope(2, 3, 3.01602932265, ScatLenFM{5.74}, 1.6, 5333 /*,1.4e-06*/);
  addIsotope(2, 4, 4.00260325413, ScatLenFM{3.26}, 0, 0 /*,0.999999*/);
  addNatElem(3, 6.94093739547, ScatLenFM{-1.9}, 0.92, 70.5);//Li
  addIsotope(3, 6, 6.01512288742, ScatLenFM{2}, 0.46, 940 /*,0.075*/);
  addIsotope(3, 7, 7.01600343666, ScatLenFM{-2.22}, 0.78, 0.0454 /*,0.925*/);
  addNatElem(4, 9.012183066, ScatLenFM{7.79}, 0.0018, 0.0076);//Be
  addIsotope(4, 9, 9.012183066, ScatLenFM{7.79}, 0.0018, 0.0076 /*,1*/);
  addNatElem(5, 10.8100315052, ScatLenFM{5.3}, 1.7, 767);//B
  addIsotope(5, 10, 10.012936862, ScatLenFM{-0.1}, 3, 3835 /*,0.2*/);
  addIsotope(5, 11, 11.009305166, ScatLenFM{6.65}, 0.21, 0.0055 /*,0.8*/);
  addNatElem(6, 12.0110369031873, ScatLenFM{6.646}, 0.001, 0.0035);//C
  addIsotope(6, 12, 12., ScatLenFM{6.6511}, 0, 0.00353 /*,0.989*/);
  addIsotope(6, 13, 13.00335483521, ScatLenFM{6.19}, 0.034, 0.00137 /*,0.011*/);
  addNatElem(7, 14.00676303357, ScatLenFM{9.36}, 0.5, 1.9);//N
  addIsotope(7, 14, 14.00307400446, ScatLenFM{9.37}, 0.5, 1.91 /*,0.9963*/);
  addIsotope(7, 15, 15.00010889894, ScatLenFM{6.44}, 5e-05, 2.4e-05 /*,0.0037*/);
  addNatElem(8, 15.999304712099, ScatLenFM{5.803}, 0.0008, 0.00019);//O
  addIsotope(8, 16, 15.9949146196, ScatLenFM{5.803}, 0, 0.0001 /*,0.99762*/);
  addIsotope(8, 17, 16.99913175664, ScatLenFM{5.78}, 0.004, 0.236 /*,0.00038*/);
  addIsotope(8, 18, 17.99915961284, ScatLenFM{5.84}, 0, 0.00016 /*,0.002*/);
  addNatElem(9, 18.99840316288, ScatLenFM{5.654}, 0.0008, 0.0096);//F
  addIsotope(9, 19, 18.99840316288, ScatLenFM{5.654}, 0.0008, 0.0096 /*,1*/);
  addNatElem(10, 20.17944669657, ScatLenFM{4.566}, 0.008, 0.039);//Ne
  addIsotope(10, 20, 19.99244017619, ScatLenFM{4.631}, 0, 0.036 /*,0.9051*/);
  addIsotope(10, 21, 20.993846685, ScatLenFM{6.66}, 0.05, 0.67 /*,0.0027*/);
  addIsotope(10, 22, 21.991385109, ScatLenFM{3.87}, 0, 0.046 /*,0.0922*/);
  addNatElem(11, 22.98976928199, ScatLenFM{3.63}, 1.62, 0.53);//Na
  addIsotope(11, 23, 22.98976928199, ScatLenFM{3.63}, 1.62, 0.53 /*,1*/);
  addNatElem(12, 24.305051619, ScatLenFM{5.375}, 0.08, 0.063);//Mg
  addIsotope(12, 24, 23.985041697, ScatLenFM{5.66}, 0, 0.05 /*,0.7899*/);
  addIsotope(12, 25, 24.985836964, ScatLenFM{3.62}, 0.28, 0.19 /*,0.1*/);
  addIsotope(12, 26, 25.982592971, ScatLenFM{4.89}, 0, 0.0382 /*,0.1101*/);
  addNatElem(13, 26.981538408, ScatLenFM{3.449}, 0.0082, 0.231);//Al
  addIsotope(13, 27, 26.981538408, ScatLenFM{3.449}, 0.0082, 0.231 /*,1*/);
  addNatElem(14, 28.0855085183, ScatLenFM{4.1491}, 0.004, 0.171);//Si
  addIsotope(14, 28, 27.97692653499, ScatLenFM{4.107}, 0, 0.177 /*,0.9223*/);
  addIsotope(14, 29, 28.97649466525, ScatLenFM{4.7}, 0.001, 0.101 /*,0.0467*/);
  addIsotope(14, 30, 29.973770136, ScatLenFM{4.58}, 0, 0.107 /*,0.031*/);
  addNatElem(15, 30.97376199863, ScatLenFM{5.13}, 0.005, 0.172);//P
  addIsotope(15, 31, 30.97376199863, ScatLenFM{5.13}, 0.005, 0.172 /*,1*/);
  addNatElem(16, 32.0643885891, ScatLenFM{2.847}, 0.007, 0.53);//S
  addIsotope(16, 32, 31.97207117443, ScatLenFM{2.804}, 0, 0.54 /*,0.9502*/);
  addIsotope(16, 33, 32.97145890985, ScatLenFM{4.74}, 0.3, 0.54 /*,0.0075*/);
  addIsotope(16, 34, 33.967867012, ScatLenFM{3.48}, 0, 0.227 /*,0.0421*/);
  addIsotope(16, 36, 35.967080699, ScatLenFM{3}, 0, 0.15 /*,0.0002*/);
  addNatElem(17, 35.4527378823, ScatLenFM{9.577}, 5.3, 33.5);//Cl
  addIsotope(17, 35, 34.968852694, ScatLenFM{11.65}, 4.7, 44.1 /*,0.7577*/);
  addIsotope(17, 37, 36.965902584, ScatLenFM{3.08}, 0.001, 0.433 /*,0.2423*/);
  addNatElem(18, 39.94766073951, ScatLenFM{1.909}, 0.225, 0.675);//Ar
  addIsotope(18, 36, 35.967545105, ScatLenFM{24.9}, 0, 5.2 /*,0.00337*/);
  addIsotope(18, 38, 37.962732104, ScatLenFM{3.5}, 0, 0.8 /*,0.00063*/);
  addIsotope(18, 40, 39.96238312378, ScatLenFM{1.83}, 0, 0.66 /*,0.996*/);
  addNatElem(19, 39.09829991492, ScatLenFM{3.67}, 0.27, 2.1);//K
  addIsotope(19, 39, 38.96370648661, ScatLenFM{3.74}, 0.25, 2.1 /*,0.93258*/);
  addIsotope(19, 40, 39.963998166, ScatLenFM{3}, 0.5, 35 /*,0.00012*/);
  addIsotope(19, 41, 40.96182525796, ScatLenFM{2.69}, 0.3, 1.46 /*,0.0673*/);
  addNatElem(20, 40.0780225128, ScatLenFM{4.7}, 0.05, 0.43);//Ca
  addIsotope(20, 40, 39.962590865, ScatLenFM{4.8}, 0, 0.41 /*,0.96941*/);
  addIsotope(20, 42, 41.958617828, ScatLenFM{3.36}, 0, 0.68 /*,0.00647*/);
  addIsotope(20, 43, 42.95876643, ScatLenFM{-1.56}, 0.5, 6.2 /*,0.00135*/);
  addIsotope(20, 44, 43.955481543, ScatLenFM{1.42}, 0, 0.88 /*,0.02086*/);
  addIsotope(20, 46, 45.95368799, ScatLenFM{3.6}, 0, 0.74 /*,4e-05*/);
  addIsotope(20, 48, 47.952522904, ScatLenFM{0.39}, 0, 1.09 /*,0.00187*/);
  addNatElem(21, 44.955907503, ScatLenFM{12.29}, 4.5, 27.5);//Sc
  addIsotope(21, 45, 44.955907503, ScatLenFM{12.29}, 4.5, 27.5 /*,1*/);
  addNatElem(22, 47.8684394371, ScatLenFM{-3.438}, 2.87, 6.09);//Ti
  addIsotope(22, 46, 45.952626856, ScatLenFM{4.93}, 0, 0.59 /*,0.082*/);
  addIsotope(22, 47, 46.951757752, ScatLenFM{3.63}, 1.5, 1.7 /*,0.074*/);
  addIsotope(22, 48, 47.947940932, ScatLenFM{-6.08}, 0, 7.84 /*,0.738*/);
  addIsotope(22, 49, 48.947864627, ScatLenFM{1.04}, 3.3, 2.2 /*,0.054*/);
  addIsotope(22, 50, 49.944785839, ScatLenFM{6.18}, 0, 0.179 /*,0.052*/);
  addNatElem(23, 50.941464864, ScatLenFM{-0.3824}, 5.08, 5.08);//V
  addIsotope(23, 50, 49.947155845, ScatLenFM{7.6}, 0.5, 60 /*,0.0025*/);
  addIsotope(23, 51, 50.943956867, ScatLenFM{-0.402}, 5.07, 4.9 /*,0.9975*/);
  addNatElem(24, 51.995920918, ScatLenFM{3.635}, 1.83, 3.05);//Cr
  addIsotope(24, 50, 49.946041443, ScatLenFM{-4.5}, 0, 15.8 /*,0.0435*/);
  addIsotope(24, 52, 51.940504992, ScatLenFM{4.92}, 0, 0.76 /*,0.8379*/);
  addIsotope(24, 53, 52.940646961, ScatLenFM{-4.2}, 5.93, 18.1 /*,0.095*/);
  addIsotope(24, 54, 53.938878012, ScatLenFM{4.55}, 0, 0.36 /*,0.0236*/);
  addNatElem(25, 54.938043172, ScatLenFM{-3.73}, 0.4, 13.3);//Mn
  addIsotope(25, 55, 54.938043172, ScatLenFM{-3.73}, 0.4, 13.3 /*,1*/);
  addNatElem(26, 55.847211691, ScatLenFM{9.45}, 0.4, 2.56);//Fe
  addIsotope(26, 54, 53.939608306, ScatLenFM{4.2}, 0, 2.25 /*,0.058*/);
  addIsotope(26, 56, 55.934935617, ScatLenFM{9.94}, 0, 2.59 /*,0.917*/);
  addIsotope(26, 57, 56.935392134, ScatLenFM{2.3}, 0.3, 2.48 /*,0.022*/);
  addIsotope(26, 58, 57.933273738, ScatLenFM{15}, 0, 1.28 /*,0.003*/);
  addNatElem(27, 58.933193656, ScatLenFM{2.49}, 4.8, 37.18);//Co
  addIsotope(27, 59, 58.933193656, ScatLenFM{2.49}, 4.8, 37.18 /*,1*/);
  addNatElem(28, 58.68788578, ScatLenFM{10.3}, 5.2, 4.49);//Ni
  addIsotope(28, 58, 57.93534178, ScatLenFM{14.4}, 0, 4.6 /*,0.6827*/);
  addIsotope(28, 60, 59.930785256, ScatLenFM{2.8}, 0, 2.9 /*,0.261*/);
  addIsotope(28, 61, 60.931054945, ScatLenFM{7.6}, 1.9, 2.5 /*,0.0113*/);
  addIsotope(28, 62, 61.928344871, ScatLenFM{-8.7}, 0, 14.5 /*,0.0359*/);
  addIsotope(28, 64, 63.927966341, ScatLenFM{-0.37}, 0, 1.52 /*,0.0091*/);
  addNatElem(29, 63.545639907, ScatLenFM{7.718}, 0.55, 3.78);//Cu
  addIsotope(29, 63, 62.929597236, ScatLenFM{6.43}, 0.006, 4.5 /*,0.6917*/);
  addIsotope(29, 65, 64.927789487, ScatLenFM{10.61}, 0.4, 2.17 /*,0.3083*/);
  addNatElem(30, 65.396361173, ScatLenFM{5.68}, 0.077, 1.11);//Zn
  addIsotope(30, 64, 63.929141772, ScatLenFM{5.22}, 0, 0.93 /*,0.486*/);
  addIsotope(30, 66, 65.926033704, ScatLenFM{5.97}, 0, 0.62 /*,0.279*/);
  addIsotope(30, 67, 66.927127482, ScatLenFM{7.56}, 0.28, 6.8 /*,0.041*/);
  addIsotope(30, 68, 67.924844291, ScatLenFM{6.03}, 0, 1.1 /*,0.188*/);
  addIsotope(30, 70, 69.925319181, ScatLenFM{6}, 0, 0.092 /*,0.006*/);
  addNatElem(31, 69.723226004, ScatLenFM{7.288}, 0.16, 2.75);//Ga
  addIsotope(31, 69, 68.925573531, ScatLenFM{7.88}, 0.091, 2.18 /*,0.601*/);
  addIsotope(31, 71, 70.924702536, ScatLenFM{6.4}, 0.084, 3.61 /*,0.399*/);
  addNatElem(32, 72.632248855, ScatLenFM{8.185}, 0.18, 2.2);//Ge
  addIsotope(32, 70, 69.924248706, ScatLenFM{10}, 0, 3 /*,0.205*/);
  addIsotope(32, 72, 71.922075826, ScatLenFM{8.51}, 0, 0.8 /*,0.274*/);
  addIsotope(32, 73, 72.923458956, ScatLenFM{5.02}, 1.5, 15.1 /*,0.078*/);
  addIsotope(32, 74, 73.921177762, ScatLenFM{7.58}, 0, 0.4 /*,0.365*/);
  addIsotope(32, 76, 75.921402726, ScatLenFM{8.2}, 0, 0.16 /*,0.078*/);
  addNatElem(33, 74.921594562, ScatLenFM{6.58}, 0.06, 4.5);//As
  addIsotope(33, 75, 74.921594562, ScatLenFM{6.58}, 0.06, 4.5 /*,1*/);
  addNatElem(34, 78.993277226, ScatLenFM{7.97}, 0.32, 11.7);//Se
  addIsotope(34, 74, 73.922475935, ScatLenFM{0.8}, 0, 51.8 /*,0.009*/);
  addIsotope(34, 76, 75.919213704, ScatLenFM{12.2}, 0, 85 /*,0.09*/);
  addIsotope(34, 77, 76.91991415, ScatLenFM{8.25}, 0.05, 42 /*,0.076*/);
  addIsotope(34, 78, 77.917309243, ScatLenFM{8.24}, 0, 0.43 /*,0.235*/);
  addIsotope(34, 80, 79.916521785, ScatLenFM{7.48}, 0, 0.61 /*,0.496*/);
  addIsotope(34, 82, 81.916699537, ScatLenFM{6.34}, 0, 0.044 /*,0.094*/);
  addNatElem(35, 79.903527044, ScatLenFM{6.795}, 0.1, 6.9);//Br
  addIsotope(35, 79, 78.918337601, ScatLenFM{6.8}, 0.15, 11 /*,0.5069*/);
  addIsotope(35, 81, 80.916288206, ScatLenFM{6.79}, 0.05, 2.7 /*,0.4931*/);
  addNatElem(36, 83.8000174955, ScatLenFM{7.81}, 0.01, 25);//Kr
  addIsotope(36, 86, 85.91061062627, ScatLenFM{8.1}, 0, 0.003 /*,0.173*/);
  addNatElem(37, 85.4676635954, ScatLenFM{7.09}, 0.5, 0.38);//Rb
  addIsotope(37, 85, 84.9117897376, ScatLenFM{7.03}, 0.5, 0.48 /*,0.7217*/);
  addIsotope(37, 87, 86.909180531, ScatLenFM{7.23}, 0.5, 0.12 /*,0.2783*/);
  addNatElem(38, 87.6166442801, ScatLenFM{7.02}, 0.06, 1.28);//Sr
  addIsotope(38, 84, 83.91341912, ScatLenFM{7}, 0, 0.87 /*,0.0056*/);
  addIsotope(38, 86, 85.90926072631, ScatLenFM{5.67}, 0, 1.04 /*,0.0986*/);
  addIsotope(38, 87, 86.90887749615, ScatLenFM{7.4}, 0.5, 16 /*,0.07*/);
  addIsotope(38, 88, 87.90561225561, ScatLenFM{7.15}, 0, 0.058 /*,0.8258*/);
  addNatElem(39, 88.905841205, ScatLenFM{7.75}, 0.15, 1.28);//Y
  addIsotope(39, 89, 88.905841205, ScatLenFM{7.75}, 0.15, 1.28 /*,1*/);
  addNatElem(40, 91.2190408226, ScatLenFM{7.16}, 0.02, 0.185);//Zr
  addIsotope(40, 90, 89.904698758, ScatLenFM{6.4}, 0, 0.011 /*,0.5145*/);
  addIsotope(40, 91, 90.905640223, ScatLenFM{8.7}, 0.15, 1.17 /*,0.1132*/);
  addIsotope(40, 92, 91.905035322, ScatLenFM{7.4}, 0, 0.22 /*,0.1719*/);
  addIsotope(40, 94, 93.906312524, ScatLenFM{8.2}, 0, 0.0499 /*,0.1728*/);
  addIsotope(40, 96, 95.908277621, ScatLenFM{5.5}, 0, 0.0229 /*,0.0276*/);
  addNatElem(41, 92.906373161, ScatLenFM{7.054}, 0.0024, 1.15);//Nb
  addIsotope(41, 93, 92.906373161, ScatLenFM{7.054}, 0.0024, 1.15 /*,1*/);
  addNatElem(42, 95.9312871581, ScatLenFM{6.715}, 0.04, 2.48);//Mo
  addIsotope(42, 92, 91.906807155, ScatLenFM{6.91}, 0, 0.019 /*,0.1484*/);
  addIsotope(42, 94, 93.905083592, ScatLenFM{6.8}, 0, 0.015 /*,0.0925*/);
  addIsotope(42, 95, 94.905837442, ScatLenFM{6.91}, 0.5, 13.1 /*,0.1592*/);
  addIsotope(42, 96, 95.904674774, ScatLenFM{6.2}, 0, 0.5 /*,0.1668*/);
  addIsotope(42, 97, 96.906016903, ScatLenFM{7.24}, 0.5, 2.5 /*,0.0955*/);
  addIsotope(42, 98, 97.905403608, ScatLenFM{6.58}, 0, 0.127 /*,0.2413*/);
  addIsotope(42, 100, 99.907467976, ScatLenFM{6.73}, 0, 0.4 /*,0.0963*/);
  addNatElem(44, 101.070135, ScatLenFM{7.03}, 0.4, 2.56);//Ru
  addNatElem(45, 102.90549407, ScatLenFM{5.88}, 0.3, 144.8);//Rh
  addIsotope(45, 103, 102.90549407, ScatLenFM{5.88}, 0.3, 144.8 /*,1*/);
  addNatElem(46, 106.41532788, ScatLenFM{5.91}, 0.093, 6.9);//Pd
  addIsotope(46, 102, 101.90563206, ScatLenFM{7.7}, 0, 3.4 /*,0.0102*/);
  addIsotope(46, 104, 103.9040304, ScatLenFM{7.7}, 0, 0.6 /*,0.1114*/);
  addIsotope(46, 105, 104.90507949, ScatLenFM{5.5}, 0.8, 20 /*,0.2233*/);
  addIsotope(46, 106, 105.90348029, ScatLenFM{6.4}, 0, 0.304 /*,0.2733*/);
  addIsotope(46, 108, 107.90389181, ScatLenFM{4.1}, 0, 8.55 /*,0.2646*/);
  addIsotope(46, 110, 109.90517287, ScatLenFM{7.7}, 0, 0.226 /*,0.1172*/);
  addNatElem(47, 107.8683298, ScatLenFM{5.922}, 0.58, 63.3);//Ag
  addIsotope(47, 107, 106.90509153, ScatLenFM{7.555}, 0.13, 37.6 /*,0.5183*/);
  addIsotope(47, 109, 108.90475577, ScatLenFM{4.165}, 0.32, 91 /*,0.4817*/);
  addNatElem(48, 112.410057988, ScatLenFM{4.87}, 3.46, 2520);//Cd
  addIsotope(48, 106, 105.9064598, ScatLenFM{5}, 0, 1 /*,0.0125*/);
  addIsotope(48, 108, 107.90418359, ScatLenFM{5.4}, 0, 1.1 /*,0.0089*/);
  addIsotope(48, 110, 109.90300746, ScatLenFM{5.9}, 0, 11 /*,0.1251*/);
  addIsotope(48, 111, 110.90418377, ScatLenFM{6.5}, 0.3, 24 /*,0.1281*/);
  addIsotope(48, 112, 111.902763883, ScatLenFM{6.4}, 0, 2.2 /*,0.2413*/);
  addIsotope(48, 113, 112.904408097, ScatLenFM{-8}, 0.3, 20600 /*,0.1222*/);
  addIsotope(48, 114, 113.90336499, ScatLenFM{7.5}, 0, 0.34 /*,0.2872*/);
  addIsotope(48, 116, 115.90476323, ScatLenFM{6.3}, 0, 0.075 /*,0.0747*/);
  addNatElem(49, 114.817886585, ScatLenFM{4.065}, 0.54, 193.8);//In
  addIsotope(49, 113, 112.904060448, ScatLenFM{5.39}, 3.7e-05, 12 /*,0.043*/);
  addIsotope(49, 115, 114.903878773, ScatLenFM{4.01}, 0.55, 202 /*,0.957*/);
  addNatElem(50, 118.68540724, ScatLenFM{6.225}, 0.022, 0.626);//Sn
  addIsotope(50, 112, 111.904824877, ScatLenFM{6}, 0, 1 /*,0.01*/);
  addIsotope(50, 114, 113.902780132, ScatLenFM{6.2}, 0, 0.114 /*,0.007*/);
  addIsotope(50, 115, 114.903344697, ScatLenFM{6}, 0.3, 30 /*,0.004*/);
  addIsotope(50, 116, 115.901742824, ScatLenFM{5.93}, 0, 0.14 /*,0.147*/);
  addIsotope(50, 117, 116.90295402, ScatLenFM{6.48}, 0.3, 2.3 /*,0.077*/);
  addIsotope(50, 118, 117.90160661, ScatLenFM{6.07}, 0, 0.22 /*,0.243*/);
  addIsotope(50, 119, 118.90331122, ScatLenFM{6.12}, 0.3, 2.2 /*,0.086*/);
  addIsotope(50, 120, 119.90220187, ScatLenFM{6.49}, 0, 0.14 /*,0.324*/);
  addIsotope(50, 122, 121.903444, ScatLenFM{5.74}, 0, 0.18 /*,0.046*/);
  addIsotope(50, 124, 123.90527669, ScatLenFM{5.97}, 0, 0.133 /*,0.056*/);
  addNatElem(51, 121.75798257, ScatLenFM{5.57}, 0.007, 4.91);//Sb
  addIsotope(51, 121, 120.90381009, ScatLenFM{5.71}, 0.0003, 5.75 /*,0.573*/);
  addIsotope(51, 123, 122.90421402, ScatLenFM{5.38}, 0.001, 3.8 /*,0.427*/);
  addNatElem(52, 127.58579825, ScatLenFM{5.8}, 0.09, 4.7);//Te
  addIsotope(52, 120, 119.90405951, ScatLenFM{5.3}, 0, 2.3 /*,0.00096*/);
  addIsotope(52, 122, 121.90304343, ScatLenFM{3.8}, 0, 3.4 /*,0.026*/);
  addIsotope(52, 123, 122.90426975, ScatLenFM{-0.05}, 0.52, 418 /*,0.00908*/);
  addIsotope(52, 124, 123.90281706, ScatLenFM{7.96}, 0, 6.8 /*,0.04816*/);
  addIsotope(52, 125, 124.9044299, ScatLenFM{5.02}, 0.008, 1.55 /*,0.0714*/);
  addIsotope(52, 126, 125.90331087, ScatLenFM{5.56}, 0, 1.04 /*,0.1895*/);
  addIsotope(52, 128, 127.90446131, ScatLenFM{5.89}, 0, 0.215 /*,0.3169*/);
  addIsotope(52, 130, 129.906222747, ScatLenFM{6.02}, 0, 0.29 /*,0.338*/);
  addNatElem(53, 126.90447184, ScatLenFM{5.28}, 0.31, 6.15);//I
  addIsotope(53, 127, 126.90447184, ScatLenFM{5.28}, 0.31, 6.15 /*,1*/);
  addNatElem(54, 131.29308175, ScatLenFM{4.92}, 0, 23.9);//Xe
  addNatElem(55, 132.905451961, ScatLenFM{5.42}, 0.21, 29);//Cs
  addIsotope(55, 133, 132.905451961, ScatLenFM{5.42}, 0.21, 29 /*,1*/);
  addNatElem(56, 137.326671887, ScatLenFM{5.07}, 0.15, 1.1);//Ba
  addIsotope(56, 130, 129.90632087, ScatLenFM{-3.6}, 0, 30 /*,0.0011*/);
  addIsotope(56, 132, 131.9050611, ScatLenFM{7.8}, 0, 7 /*,0.001*/);
  addIsotope(56, 134, 133.904508399, ScatLenFM{5.7}, 0, 2 /*,0.0242*/);
  addIsotope(56, 135, 134.905688606, ScatLenFM{4.67}, 0.5, 5.8 /*,0.0659*/);
  addIsotope(56, 136, 135.904575959, ScatLenFM{4.91}, 0, 0.68 /*,0.0785*/);
  addIsotope(56, 137, 136.905827375, ScatLenFM{6.83}, 0.5, 3.6 /*,0.1123*/);
  addIsotope(56, 138, 137.905247229, ScatLenFM{4.84}, 0, 0.27 /*,0.717*/);
  addNatElem(57, 138.90545949, ScatLenFM{8.24}, 1.13, 8.97);//La
  addIsotope(57, 138, 137.90711783, ScatLenFM{8}, 0.5, 57 /*,0.0009*/);
  addIsotope(57, 139, 138.9063588, ScatLenFM{8.24}, 1.13, 8.93 /*,0.9991*/);
  addNatElem(58, 140.1148724, ScatLenFM{4.84}, 0.001, 0.63);//Ce
  addIsotope(58, 136, 135.90712944, ScatLenFM{5.8}, 0, 7.3 /*,0.0019*/);
  addIsotope(58, 138, 137.9059887, ScatLenFM{6.7}, 0, 1.1 /*,0.0025*/);
  addIsotope(58, 140, 139.90544642, ScatLenFM{4.84}, 0, 0.57 /*,0.8848*/);
  addIsotope(58, 142, 141.90924988, ScatLenFM{4.75}, 0, 0.95 /*,0.1108*/);
  addNatElem(59, 140.9076584, ScatLenFM{4.58}, 0.015, 11.5);//Pr
  addIsotope(59, 141, 140.9076584, ScatLenFM{4.58}, 0.015, 11.5 /*,1*/);
  addNatElem(60, 144.24064436, ScatLenFM{7.69}, 9.2, 50.5);//Nd
  addIsotope(60, 142, 141.90772889, ScatLenFM{7.7}, 0, 18.7 /*,0.2716*/);
  addIsotope(60, 143, 142.90981989, ScatLenFM{14}, 55, 337 /*,0.1218*/);
  addIsotope(60, 144, 143.91009286, ScatLenFM{2.8}, 0, 3.6 /*,0.238*/);
  addIsotope(60, 145, 144.9125792, ScatLenFM{14}, 5, 42 /*,0.0829*/);
  addIsotope(60, 146, 145.9131225, ScatLenFM{8.7}, 0, 1.4 /*,0.1719*/);
  addIsotope(60, 148, 147.91689909, ScatLenFM{5.7}, 0, 2.5 /*,0.0575*/);
  addIsotope(60, 150, 149.92090153, ScatLenFM{5.3}, 0, 1.2 /*,0.0563*/);
  addNatElem(62, 150.35023831, ScatLenFM{0.8}, 39, 5922);//Sm
  addIsotope(62, 144, 143.91200637, ScatLenFM{-3}, 0, 0.7 /*,0.031*/);
  addIsotope(62, 147, 146.91490406, ScatLenFM{14}, 143, 57 /*,0.151*/);
  addIsotope(62, 148, 147.91482901, ScatLenFM{-3}, 0, 2.4 /*,0.113*/);
  addIsotope(62, 149, 148.91719137, ScatLenFM{-19.2}, 137, 42080 /*,0.139*/);
  addIsotope(62, 150, 149.91728219, ScatLenFM{14}, 0, 104 /*,0.074*/);
  addIsotope(62, 152, 151.91973904, ScatLenFM{-5}, 0, 206 /*,0.266*/);
  addIsotope(62, 154, 153.92221616, ScatLenFM{9.3}, 0, 8.4 /*,0.226*/);
  addNatElem(63, 151.96457732, ScatLenFM{7.22}, 2.5, 4530);//Eu
  addIsotope(63, 151, 150.91985686, ScatLenFM{6.13}, 3.1, 9100 /*,0.478*/);
  addIsotope(63, 153, 152.92123704, ScatLenFM{8.22}, 1.3, 312 /*,0.522*/);
  addNatElem(64, 157.2510281, ScatLenFM{6.5}, 151, 49700);//Gd
  addIsotope(64, 152, 151.91979882, ScatLenFM{10}, 0, 735 /*,0.002*/);
  addIsotope(64, 154, 153.9208734, ScatLenFM{10}, 0, 85 /*,0.021*/);
  addIsotope(64, 155, 154.9226298, ScatLenFM{6}, 25, 61100 /*,0.148*/);
  addIsotope(64, 156, 155.92213056, ScatLenFM{6.3}, 0, 1.5 /*,0.206*/);
  addIsotope(64, 157, 156.92396787, ScatLenFM{-1.14}, 394, 259000 /*,0.157*/);
  addIsotope(64, 158, 157.92411165, ScatLenFM{9}, 0, 2.2 /*,0.248*/);
  addIsotope(64, 160, 159.92706154, ScatLenFM{9.15}, 0, 0.77 /*,0.218*/);
  addNatElem(65, 158.92535393, ScatLenFM{7.38}, 0.004, 23.4);//Tb
  addIsotope(65, 159, 158.92535393, ScatLenFM{7.38}, 0.004, 23.4 /*,1*/);
  addNatElem(66, 162.49453743, ScatLenFM{16.9}, 54.4, 994);//Dy
  addIsotope(66, 156, 155.92428404, ScatLenFM{6.1}, 0, 33 /*,0.0006*/);
  addIsotope(66, 158, 157.9244146, ScatLenFM{6}, 0, 43 /*,0.001*/);
  addIsotope(66, 160, 159.92520324, ScatLenFM{6.7}, 0, 56 /*,0.0234*/);
  addIsotope(66, 161, 160.92693909, ScatLenFM{10.3}, 3, 600 /*,0.19*/);
  addIsotope(66, 162, 161.92680417, ScatLenFM{-1.4}, 0, 194 /*,0.255*/);
  addIsotope(66, 163, 162.92873688, ScatLenFM{5}, 0.21, 124 /*,0.249*/);
  addIsotope(66, 164, 163.92918047, ScatLenFM{49.4}, 0, 2840 /*,0.281*/);
  addNatElem(67, 164.93032805, ScatLenFM{8.01}, 0.36, 64.7);//Ho
  addIsotope(67, 165, 164.93032805, ScatLenFM{8.01}, 0.36, 64.7 /*,1*/);
  addNatElem(68, 167.26221528, ScatLenFM{7.79}, 1.1, 159);//Er
  addIsotope(68, 162, 161.92878696, ScatLenFM{8.8}, 0, 19 /*,0.0014*/);
  addIsotope(68, 164, 163.92920739, ScatLenFM{8.2}, 0, 13 /*,0.0156*/);
  addIsotope(68, 166, 165.93029902, ScatLenFM{10.6}, 0, 19.6 /*,0.334*/);
  addIsotope(68, 167, 166.93205412, ScatLenFM{3}, 0.13, 659 /*,0.229*/);
  addIsotope(68, 168, 167.93237619, ScatLenFM{7.4}, 0, 2.74 /*,0.271*/);
  addIsotope(68, 170, 169.93547067, ScatLenFM{9.6}, 0, 5.8 /*,0.149*/);
  addNatElem(69, 168.93421835, ScatLenFM{7.07}, 0.1, 100);//Tm
  addIsotope(69, 169, 168.93421835, ScatLenFM{7.07}, 0.1, 100 /*,1*/);
  addNatElem(70, 173.0333950863, ScatLenFM{12.43}, 4, 34.8);//Yb
  addIsotope(70, 168, 167.93388911, ScatLenFM{-4.07}, 0, 2230 /*,0.0014*/);
  addIsotope(70, 170, 169.934767245, ScatLenFM{6.77}, 0, 11.4 /*,0.0306*/);
  addIsotope(70, 171, 170.936331517, ScatLenFM{9.66}, 3.9, 48.6 /*,0.143*/);
  addIsotope(70, 172, 171.936386658, ScatLenFM{9.43}, 0, 0.8 /*,0.219*/);
  addIsotope(70, 173, 172.938216215, ScatLenFM{9.56}, 3.5, 17.1 /*,0.161*/);
  addIsotope(70, 174, 173.938867548, ScatLenFM{19.3}, 0, 69.4 /*,0.318*/);
  addIsotope(70, 176, 175.942574708, ScatLenFM{8.72}, 0, 2.85 /*,0.127*/);
  addNatElem(71, 174.96692728, ScatLenFM{7.21}, 0.7, 74);//Lu
  addIsotope(71, 175, 174.94077731, ScatLenFM{7.24}, 0.6, 21 /*,0.9739*/);
  addIsotope(71, 176, 175.94269181, ScatLenFM{6.1}, 1.2, 2065 /*,0.0261*/);
  addNatElem(72, 178.48778639, ScatLenFM{7.7}, 2.6, 104.1);//Hf
  addIsotope(72, 174, 173.94004848, ScatLenFM{10.9}, 0, 561 /*,0.002*/);
  addIsotope(72, 176, 175.94140991, ScatLenFM{6.61}, 0, 23.5 /*,0.052*/);
  addIsotope(72, 177, 176.94323032, ScatLenFM{0.8}, 0.1, 373 /*,0.186*/);
  addIsotope(72, 178, 177.94370846, ScatLenFM{5.9}, 0, 84 /*,0.271*/);
  addIsotope(72, 179, 178.94582584, ScatLenFM{7.46}, 0.14, 41 /*,0.137*/);
  addIsotope(72, 180, 179.94655967, ScatLenFM{13.2}, 0, 13.04 /*,0.352*/);
  addNatElem(73, 180.94787927, ScatLenFM{6.91}, 0.01, 20.6);//Ta
  addIsotope(73, 180, 179.94746839, ScatLenFM{7}, 0.5, 563 /*,0.00012*/);
  addIsotope(73, 181, 180.94799933, ScatLenFM{6.91}, 0.011, 20.5 /*,0.99988*/);
  addNatElem(74, 183.85009188, ScatLenFM{4.86}, 1.63, 18.3);//W
  addIsotope(74, 180, 179.94671343, ScatLenFM{5}, 0, 30 /*,0.001*/);
  addIsotope(74, 182, 181.94820572, ScatLenFM{6.97}, 0, 20.7 /*,0.263*/);
  addIsotope(74, 183, 182.9502245, ScatLenFM{6.53}, 0.3, 10.1 /*,0.143*/);
  addIsotope(74, 184, 183.95093326, ScatLenFM{7.48}, 0, 1.7 /*,0.307*/);
  addIsotope(74, 186, 185.95436521, ScatLenFM{-0.72}, 0, 37.9 /*,0.286*/);
  addNatElem(75, 186.20670735, ScatLenFM{9.2}, 0.9, 89.7);//Re
  addIsotope(75, 185, 184.95295834, ScatLenFM{9}, 0.5, 112 /*,0.374*/);
  addIsotope(75, 187, 186.95575229, ScatLenFM{9.3}, 1, 76.4 /*,0.626*/);
  addNatElem(76, 190.23977696, ScatLenFM{10.7}, 0.3, 16);//Os
  addIsotope(76, 184, 183.95249295, ScatLenFM{10}, 0, 3000 /*,0.0002*/);
  addIsotope(76, 186, 185.95383766, ScatLenFM{11.6}, 0, 80 /*,0.0158*/);
  addIsotope(76, 187, 186.95574964, ScatLenFM{10}, 0.3, 320 /*,0.016*/);
  addIsotope(76, 188, 187.95583736, ScatLenFM{7.6}, 0, 4.7 /*,0.133*/);
  addIsotope(76, 189, 188.958146, ScatLenFM{10.7}, 0.5, 25 /*,0.161*/);
  addIsotope(76, 190, 189.9584455, ScatLenFM{11}, 0, 13.1 /*,0.264*/);
  addIsotope(76, 192, 191.96147888, ScatLenFM{11.5}, 0, 2 /*,0.41*/);
  addNatElem(77, 192.21605388, ScatLenFM{10.6}, 0, 425);//Ir
  addNatElem(78, 195.0801337, ScatLenFM{9.6}, 0.13, 10.3);//Pt
  addIsotope(78, 190, 189.95994988, ScatLenFM{9}, 0, 152 /*,0.0001*/);
  addIsotope(78, 192, 191.96104274, ScatLenFM{9.9}, 0, 10 /*,0.0079*/);
  addIsotope(78, 194, 193.962683527, ScatLenFM{10.55}, 0, 1.44 /*,0.329*/);
  addIsotope(78, 195, 194.964794353, ScatLenFM{8.83}, 0.13, 27.5 /*,0.338*/);
  addIsotope(78, 196, 195.964954675, ScatLenFM{9.89}, 0, 0.72 /*,0.253*/);
  addIsotope(78, 198, 197.96789673, ScatLenFM{7.8}, 0, 3.66 /*,0.072*/);
  addNatElem(79, 196.966570114, ScatLenFM{7.63}, 0.43, 98.65);//Au
  addIsotope(79, 197, 196.966570114, ScatLenFM{7.63}, 0.43, 98.65 /*,1*/);
  addNatElem(80, 200.58545474, ScatLenFM{12.692}, 6.6, 372.3);//Hg
  addIsotope(80, 196, 195.96583344, ScatLenFM{30.3}, 0, 3080 /*,0.002*/);
  addIsotope(80, 199, 198.968280989, ScatLenFM{16.9}, 30, 2150 /*,0.17*/);
  addNatElem(81, 204.38333219, ScatLenFM{8.776}, 0.21, 3.43);//Tl
  addIsotope(81, 203, 202.97234402, ScatLenFM{6.99}, 0.14, 11.4 /*,0.29524*/);
  addIsotope(81, 205, 204.97442724, ScatLenFM{9.52}, 0.007, 0.104 /*,0.70476*/);
  addNatElem(82, 207.21690749, ScatLenFM{9.405}, 0.003, 0.171);//Pb
  addIsotope(82, 204, 203.97304342, ScatLenFM{9.9}, 0, 0.65 /*,0.014*/);
  addIsotope(82, 206, 205.97446512, ScatLenFM{9.22}, 0, 0.03 /*,0.241*/);
  addIsotope(82, 207, 206.97589673, ScatLenFM{9.28}, 0.002, 0.699 /*,0.221*/);
  addIsotope(82, 208, 207.97665192, ScatLenFM{9.5}, 0, 0.00048 /*,0.524*/);
  addNatElem(83, 208.98039852, ScatLenFM{8.532}, 0.0084, 0.0338);//Bi
  addIsotope(83, 209, 208.98039852, ScatLenFM{8.532}, 0.0084, 0.0338 /*,1*/);
  addNatElem(90, 232.03805369, ScatLenFM{10.31}, 0, 7.37);//Th
  addIsotope(90, 232, 232.03805369, ScatLenFM{10.31}, 0, 7.37 /*,1*/);
  addNatElem(92, 238.02893712, ScatLenFM{8.417}, 0.005, 7.57);//U
  addIsotope(92, 233, 233.03963437, ScatLenFM{10.1}, 0.1, 574.7);
  addIsotope(92, 234, 234.04095037, ScatLenFM{12.4}, 0, 100.1 /*,5e-05*/);
  addIsotope(92, 235, 235.04392819, ScatLenFM{10.47}, 0.2, 680.9 /*,0.0072*/);
  addIsotope(92, 238, 238.050787, ScatLenFM{8.402}, 0, 2.68 /*,0.99275*/);

  // ----- Autogenerated code end (no manual edits inside please!!!) ----- //


  //Fix 1) Four elements added here on the bottom which were available in NCrystal
  //until v2.0.0, but were not included in the auto-generated code above. We
  //manually copy them over here from the old code, to avoid elements suddenly
  //becoming unavailable with NCrystal v2.1.0:
  addNatElem(43, 98, ScatLenFM{6.8}, 0.5, 20.0);//Tc
  addNatElem(61, 145, ScatLenFM{12.6}, 1.3, 168.4);//Pm
  addNatElem(88, 226, ScatLenFM{10}, 0.0, 12.8);//Ra
  addNatElem(91, 231.03588, ScatLenFM{9.1}, 0.1, 200.6);//Pa

  //Sort the database:
  std::sort(result.begin(),result.end());//sort (and verify uniqueness in dbg builds)
  nc_assert(result.size()==numberOfEntriesCountInDB);//if this fails, please update numberOfEntriesCountInDB.

  //Fix 2) We modify the incoherent cross-section of Xenon. It appears to have
  //been generated with a value of 0 barn, even though it has a scattering
  //length around 3.04fm -> 1.16133 barn. In cif2hkl Xenon also has 1.16133
  //barn, so it seems consistent:

  auto itXenon = std::lower_bound(result.begin(),result.end(),AtomDBKey(54,0));
  nc_assert( itXenon!=result.end() && itXenon->Z()==54 && itXenon->A() == 0);
  itXenon->modifyIncXS(1.16133);

  //Fix 3) Update the coherent scattering length of Calcium with the value
  //determined in the recent (2019) paper DOI:10.1016/j.nima.2019.02.072 The
  //values above are b_coh=4.7fm (~2.77barn) and sigma_incoh=0.05barn. The 2019
  //paper reports b_coh=5.83fm+-0.03fm (4.27barn), but does not mention the
  //incoherent contribution. As a compromise we keep sigma_incoh=0.05 and fix
  //b_coh so sigma_coh+sigma_incoh=4.27barn, meaning b_coh = 5.795 (which is
  //within 1 sigma of 5.83). Concerning isotope-specific data, we update just
  //Ca40, since this isotope is predominant (~97%) in natural calcium, picking a
  //value (5.909624496225fm) which ensures that mixing all isotopes yield a
  //scattering xs of 4.27barn.
  auto itCalcium = std::lower_bound(result.begin(),result.end(),AtomDBKey(20,0));
  nc_assert( itCalcium!=result.end() && itCalcium->Z()==20 && itCalcium->A() == 0);
  itCalcium->modifyCohScatLen(ScatLenFM{5.795});//increased from 4.7, an increase of 23.3%
  auto itCalcium40 = std::lower_bound(result.begin(),result.end(),AtomDBKey(20,40));
  nc_assert( itCalcium40!=result.end() && itCalcium40->Z()==20 && itCalcium40->A() == 40);
  itCalcium40->modifyCohScatLen(ScatLenFM{5.909624496225});//increased from 4.8, an increase of 23.1%


  if (ncgetenv_bool("ATOMDB_DUMP")) {
    std::ostringstream ss;
    ss<<"NCrystal::AtomDB:BeginDump (since NCRYSTAL_ATOMDB_DUMP env var was set)\n";
    for (const auto& e : result)
      ss<<"NCrystal::AtomDB:  "<<e.getAtomDBLine()<<'\n';
    ss<<"NCrystal::AtomDB:EndDump\n";
    Msg::outputMsg(ss.str(),MsgType::RawOutput);
  }

  return result;
}

namespace NCRYSTAL_NAMESPACE {
  namespace AtomDB {
    static bool dummy_internal_force_init = [] () {
      if (ncgetenv_bool("ATOMDB_DUMP"))
        internal::internalDB();//ensure internalDB initialisation
      return true;
    }();
  }
}
