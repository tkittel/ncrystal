#ifndef NCrystal_RandUtils_hh
#define NCrystal_RandUtils_hh

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

#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //Directions or scattering angles (mu=cos(theta_scat)) under full or partial
  //isotropic distributions. Vector results are unit vectors:
  double randIsotropicScatterAngle( RNG& );//deprecated!
  NCRYSTAL_API CosineScatAngle randIsotropicScatterMu( RNG& );//marked NCRYSTAL_API since used in custom physics example
  Vector randIsotropicDirection( RNG& );
  Vector randDirectionGivenScatterMu( RNG&, double mu, const Vector& in );

  NeutronDirection randIsotropicNeutronDirection( RNG& );
  NeutronDirection randNeutronDirectionGivenScatterMu( RNG&, double mu, const Vector& in );
  inline NeutronDirection randNeutronDirectionGivenScatterMu( RNG& rng, CosineScatAngle mu, const NeutronDirection& in )
  {
    return randNeutronDirectionGivenScatterMu( rng, mu.dbl(), in.as<Vector>() );
  }

  //Sample uniformly in an interval:
  double randInterval( double rand01, double a, double b );
  double randInterval( RNG& rng, double a, double b );

  //Sample a random point on the unit circle:
  PairDD randPointOnUnitCircle( RNG& );

  //Sample a random point inside the unit circle:
  PairDD randPointInUnitCircle( RNG& );

  //Sample one or two independent values from a unit Gaussian:
  double randNorm( RNG& );
  void randNorm( RNG&, double&g1, double&g2);

  //sample gaussian tail (tail>=0!), like sampling randNorm until result is
  //>=tail (but more efficient):
  double randNormTail(double tail, RNG& rng);

  class RandIdxByWeightPicker final : private MoveOnly {
  public:

    //Utility for picking random indices with probabilities proportional to
    //weights stored in a commulative list.
    //
    //For weights which change between picks, or for small arrays, use the static
    //oneShotPick: it needs no set-up, and costs O(log n) per pick.
    //
    //For fixed weights used for many picks, construct an instance instead: it
    //costs O(n) once, after which each pick is O(1) on average however uneven
    //the weights.
    //
    //All methods give identical results for identical input, and uses a single
    //random number in case of multiple candidates
    //
    //Refer to the implementation for details.

    static std::size_t oneShotPick( RNG&, Span<const double> commulvals );
    static std::size_t oneShotPick( double rand01val,
                                    Span<const double> commulvals );

    RandIdxByWeightPicker() = default;//empty, can not be used for picking
    explicit RandIdxByWeightPicker( VectD&& commulvals );

    std::size_t pick( RNG& ) const;
    std::size_t pick( double rand01val ) const;
    std::size_t pick( RNG&, std::size_t n ) const;
    std::size_t pick( double rand01val, std::size_t n ) const;

    std::size_t size() const noexcept { return m_v.size(); }

  private:
    VectD m_v;
    std::vector<std::uint32_t> m_guide;
    double m_lo = 0.0;
    double m_step = 0.0;
    double m_invStep = 0.0;
    double guideVal( std::size_t k ) const noexcept { return m_lo + k * m_step; }
    std::size_t lowerBoundIdx( double val ) const noexcept;
  };

  //FIXME: Kept for the ncplugin-MagScat plugin (NCParamagneticScatter.cc), which
  //should be updated (PR) to use RandIdxByWeightPicker::oneShotPick before the
  //release, after which these can be removed:
  std::size_t pickRandIdxByWeight( RNG&, Span<const double> commulvals);
  std::size_t pickRandIdxByWeight( double rand01val, Span<const double> commulvals);

  //Sample exponential:
  double randExp( RNG& rng );//sample non-negative value from exp(-x)
  double randExpInterval( RNG& rng, double a, double b, double c );//sample value in [a,b] from exp(-c*x)

  //Sample a value x in [0,1] for a PDF proportional to k^x for k>0, based on a
  //value uniformly distributed in [0,1], rand01val. This can be used to sample
  //a PDF whose density is linear in log and changes by a factor of k between 0
  //and 1. To allow caching, lnk=std::log(k) must be provided by the caller.
  double randKPowX( double k, double lnk, double rand01val );

  class RandExpIntervalSampler {
  public:
    //Like randExpInterval fct, but more efficient if need to sample repeatedly
    //with same parameters.
    RandExpIntervalSampler();//invalid instance, must complete with set(..)
    void set(double a, double b, double c);
    RandExpIntervalSampler(double a, double b, double c);
    double sample(RNG&rng) const;
    double sample( double randval01 ) const;
    void invalidate();
    bool isValid() const;
  private:
    double m_a, m_c1, m_c2;
  };

  //Sample f(x) = exp(-c*x)/sqrt(x) on [a,b], a>=0 b>a, c>0:
  double randExpDivSqrt( RNG&, double c, double a, double b );

  class RandXRSRImpl final : private MoveOnly {
    //Generator implementing the xoroshiro128+ (XOR/rotate/shift/rotate) due to
    //David Blackman and Sebastiano Vigna (released into public domain / CC0
    //1.0). It has a period of 2^128-1, is very fast and passes most statistical
    //tests. The one exception is that the two lowest order bits of the directly
    //generated 64 bit integers are not of high quality, but we solve that by
    //internally generating two numbers (i.e. 128 bits) in genUInt64, and using
    //only the good quality bits in the result.
  public:
    using state_t = std::array<uint64_t,2>;
    RandXRSRImpl(uint64_t seed = 0);//NB: seed = 0 is not a special seed value.
    RandXRSRImpl( no_init_t ) {}//invalid unspecified state
#if defined(__GNUC__) && (__GNUC__*1000+__GNUC_MINOR__)<7000
    RandXRSRImpl( const state_t& st ) { m_s[0] = st[0]; m_s[1] = st[1]; }
#else
    RandXRSRImpl( const state_t& st ) : m_s{st} {}
#endif

    double generate() noexcept;// uniformly in ]0,1]
    void generateMany(std::size_t n, double* tgt)  ncnoexceptndebug;

    uint64_t genUInt64() noexcept;//uniformly over 0..uint64max (i.e. all bits randomised)
    uint32_t genUInt32() noexcept;//uniformly over 0..uint32max (i.e. all bits randomised)
    bool coinflip() noexcept;

    void seed(uint64_t seed);
    const state_t& state() const noexcept  { return m_s; }
    state_t& state() noexcept { return m_s; }
    void jump();

    //Internal functions, exposed solely for the purpose of unit tests:
    static uint64_t splitmix64(uint64_t& state);
    uint64_t genUInt64WithBadLowerBits() noexcept;
  private:
    state_t m_s;
    void genmanyimpl(int n, double* tgt)  ncnoexceptndebug;

  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline std::size_t NCrystal::RandIdxByWeightPicker::oneShotPick( RNG& rng,
                                                                 Span<const double> commulvals )
{
  return commulvals.size()==1 ? 0 : oneShotPick( rng.generate(), commulvals );
}

inline std::size_t NCrystal::RandIdxByWeightPicker::pick( RNG& rng ) const
{
  return pick( rng, m_v.size() );
}

inline std::size_t NCrystal::RandIdxByWeightPicker::pick( double rand01val ) const
{
  return pick( rand01val, m_v.size() );
}

inline std::size_t NCrystal::RandIdxByWeightPicker::pick( RNG& rng,
                                                          std::size_t n ) const
{
  return n==1 ? 0 : pick( rng.generate(), n );
}

inline std::size_t NCrystal::pickRandIdxByWeight( RNG& rng, Span<const double> commulvals)
{
  return RandIdxByWeightPicker::oneShotPick( rng, commulvals );
}

inline double NCrystal::randIsotropicScatterAngle( NCrystal::RNG& rng )
{
  return std::acos(randIsotropicScatterMu(rng).get());
}

inline NCrystal::CosineScatAngle NCrystal::randIsotropicScatterMu( RNG& rng )
{
  const double mu = -1.0+rng.generate()*2.0;
  nc_assert(mu>=-1.0&&mu<=1.0);
  return CosineScatAngle{mu};
}

inline NCrystal::NeutronDirection NCrystal::randIsotropicNeutronDirection( RNG& rng )
{
  return randIsotropicDirection(rng).as<NeutronDirection>();
}

inline NCrystal::NeutronDirection NCrystal::randNeutronDirectionGivenScatterMu( RNG& rng, double mu, const Vector& in )
{
  return randDirectionGivenScatterMu(rng,mu,in).as<NeutronDirection>();
}

inline double NCrystal::randExp( RNG& rng )
{
  return - std::log(rng.generate());
}

inline NCrystal::RandExpIntervalSampler::RandExpIntervalSampler() : m_a(0), m_c1(0), m_c2(0)
{
}

inline void NCrystal::RandExpIntervalSampler::set(double a, double b, double c)
{
  nc_assert( c > 0.0 );
  m_a = a;
  m_c1 = -1.0/c;
  m_c2 = std::expm1(-c*(b-a));
  nc_assert(isValid());
}

inline NCrystal::RandExpIntervalSampler::RandExpIntervalSampler(double a, double b, double c)
  : m_a(a), m_c1(-1.0/c), m_c2(std::expm1(-c*(b-a)))
{
  nc_assert(isValid());
}

inline void NCrystal::RandExpIntervalSampler::invalidate()
{
  m_a = m_c1 = m_c2 = 0.0;
}

inline bool NCrystal::RandExpIntervalSampler::isValid() const
{
  return m_c1 < 0.0;
}

inline double NCrystal::RandExpIntervalSampler::sample(RNG&rng) const
{
  nc_assert(isValid());
  return m_a + m_c1 * std::log( 1.0 + rng.generate() * m_c2 );
}

inline double NCrystal::RandExpIntervalSampler::sample( double randval01 ) const
{
  nc_assert(isValid());
  return m_a + m_c1 * std::log( 1.0 + randval01 * m_c2 );
}

inline double NCrystal::randExpInterval( RNG& rng, double a, double b, double c )
{
  return RandExpIntervalSampler(a,b,c).sample(rng);
}

inline uint64_t NCrystal::RandXRSRImpl::genUInt64WithBadLowerBits() noexcept
{
  const uint64_t s0 = m_s[0];
  uint64_t s1 = m_s[1];
  uint64_t result = s0 + s1;
  s1 ^= s0;
  m_s[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
  m_s[1] = (s1 << 36) | (s1 >> 28);
  return result;
}

inline double NCrystal::RandXRSRImpl::generate() noexcept
{
  return randUInt64ToFP01( genUInt64WithBadLowerBits() );
}

inline uint64_t NCrystal::RandXRSRImpl::genUInt64() noexcept
{
  //Since lower 3 bits in generator output have unwanted correlations, we simply
  //combine upper bits of two integers into one:
  const uint64_t g1 = genUInt64WithBadLowerBits();
  const uint64_t g2 = genUInt64WithBadLowerBits();
  return ( g1 >> 32 ) | ( g2 & 0xffffffff00000000 );
}

inline uint32_t NCrystal::RandXRSRImpl::genUInt32() noexcept
{
  //Since lower 3 bits in generator output have unwanted correlations, we simply
  //just use the upper bits:
  return static_cast<uint32_t>(genUInt64WithBadLowerBits() >> 32);
}

inline bool NCrystal::RandXRSRImpl::coinflip() noexcept
{
  //Test one of the high bits, stay far away from the 3 lowest:
  constexpr uint64_t onebit = 0x1000000000000000ull;
  return onebit & genUInt64WithBadLowerBits();
}

inline NCrystal::PairDD NCrystal::randPointInUnitCircle( RNG& rng )
{
  double a,b;
  do {
    a = -1.0+rng.generate()*2.0;
    b = -1.0+rng.generate()*2.0;
  } while ( a*a + b*b > 1.0 );
  return { a, b };
}


inline double NCrystal::randInterval( double rand01, double a, double b )
{
  nc_assert( rand01 >= 0.0 && rand01 <= 1.0 );
  return intervalPos01( a, b, rand01 );
}

inline double NCrystal::randInterval( RNG& rng, double a, double b )
{
  return randInterval( rng.generate(), a, b );
}

#endif
