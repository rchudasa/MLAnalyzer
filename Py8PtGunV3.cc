#include "GeneratorInterface/Core/interface/GeneratorFilter.h"
#include "GeneratorInterface/ExternalDecays/interface/ExternalDecayDriver.h"

#include "GeneratorInterface/Pythia8Interface/interface/Py8GunBase.h"

#include <numeric>

namespace gen {

class Py8PtGunV3 : public Py8GunBase {

   public:

      Py8PtGunV3( edm::ParameterSet const& );
      ~Py8PtGunV3() {}

      bool generatePartonsAndHadronize() override;
      const char* classname() const override;

   private:

      // PtGun particle(s) characteristics
      double  fMinEta;
      double  fMaxEta;
      double  fMinPt ;
      double  fMaxPt ;
      double  fMinMass ;
      double  fMaxMass ;
      bool    fAddAntiParticle;

};

// implementation
//
Py8PtGunV3::Py8PtGunV3( edm::ParameterSet const& ps )
   : Py8GunBase(ps)
{

   // ParameterSet defpset ;
   edm::ParameterSet pgun_params =
      ps.getParameter<edm::ParameterSet>("PGunParameters"); // , defpset ) ;
   fMinEta     = pgun_params.getParameter<double>("MinEta"); // ,-2.2);
   fMaxEta     = pgun_params.getParameter<double>("MaxEta"); // , 2.2);
   fMinPt      = pgun_params.getParameter<double>("MinPt"); // ,  0.);
   fMaxPt      = pgun_params.getParameter<double>("MaxPt"); // ,  0.);
   fMinMass    = pgun_params.getParameter<double>("MinMass"); // ,  0.);
   fMaxMass    = pgun_params.getParameter<double>("MaxMass"); // ,  0.);
   fAddAntiParticle = pgun_params.getParameter<bool>("AddAntiParticle"); //, false) ;

}

int sum(std::vector <int> dist) {
    return std::accumulate(dist.begin(), dist.end(), 0);
}

double max_element(std::vector <double> dist) {
    double max = 0;
    int s = dist.size();
    for (int i = 0; i < s; i++) {
        double el = dist[i];
        if (max < el){max = el;}
    }
    return max;
}

std::vector <double> get_inverse_pdf(std::vector <int> dist) {
    std::vector <double> invpdf(dist.size());
    double sum_hist = sum(dist);
    int s = dist.size();
    for (int i = 0; i < s; i++) {
        if (dist[i] != 0 ) {
            invpdf[i] = sum_hist / dist[i];
            //std::cout << "Bin " << i << " -> " << invpdf[i] << std::endl;
        }
        else {invpdf[i] = 1;}
    }
    double max_invpdf = max_element(invpdf);
    for (int i = 0; i < s; i++) {
        invpdf[i] = invpdf[i] / max_invpdf;
    }
    return invpdf;
}

double lookup_mass_invpdf(double mgen, std::vector <double> m_bins, std::vector <double> m_invpdf) {
    int im = 0;
    int s1 = m_bins.size();
    int s2 = m_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
        im = ib;
        if (ib + 1 >  s2 - 1) { break; }
        if (mgen <= m_bins[ib]) { break; }
    }
    return m_invpdf[im];
}

double lookup_pt_invpdf(double pTgen, std::vector <int> pT_bins, std::vector <double> pT_invpdf) {
    int ipt = 0;
    int s1 = pT_bins.size();
    int s2 = pT_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
        ipt = ib;
        if (ib + 1 >  s2 - 1) { break; }
        if (pTgen <= pT_bins[ib]) { break; }
    }
    return pT_invpdf[ipt];
}

double lookup_invpdf(double Mgen, std::vector <double> M_bins, int pTgen, std::vector <int> pT_bins, std::vector <double> invpdf) {
    unsigned int ibin = 0;
    unsigned int m1  = M_bins.size();
    unsigned int pt1 = pT_bins.size();
    unsigned int inv = invpdf.size();
    bool found_mass = false;
    bool found_end  = false;
    for (unsigned int ibx = 0; ibx < m1; ibx++) {
        if (found_mass || found_end) { break; }
        for (unsigned int iby = 0; iby < pt1; iby++) {
            ibin = (ibx*pt1)+ iby;
            if ( ((ibx*pt1) + iby + 1) >  (inv - 1) ) {
                found_end = true;
                break;
            }
            if ( (Mgen  <= M_bins[ibx]) && (pTgen <= pT_bins[iby]) ) {
                found_mass = true;
                break;
            }
        }
    }
    return invpdf[ibin];
}

double get_rand_el(std::vector <int> dist) {
    int randomIndex = rand() % dist.size();
      return dist[randomIndex];
}

std::vector <int> pT_bins   = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120};
std::vector <double> m_bins = {4.125, 4.65, 5.175000000000001, 5.7, 6.225, 6.75, 7.275, 7.800000000000001, 8.325000000000001, 8.85, 9.375, 9.9, 10.425, 10.950000000000001, 11.475, 12.0, 14.0};
std::vector <int> occ = {
15451, 19638, 22435, 24176, 25277, 26082, 26399, 27387, 27441, 27218, 27761, 27426, 27313, 27411, 27802, 27741, 27428, 27810, // 3.6 - 4.125
17314, 21942, 24582, 26477, 28168, 28856, 29235, 29405, 30037, 30239, 30498, 30520, 30775, 30508, 31096, 31015, 30916, 31032, // 4.125
16994, 21427, 24404, 26594, 27888, 28867, 29615, 29369, 30114, 30044, 30229, 30642, 30664, 30811, 30852, 30686, 30732, 31364, // 4.65
16788, 21219, 24544, 26520, 27865, 28559, 29063, 29879, 29835, 29523, 29974, 30281, 30452, 30620, 30885, 30785, 31114, 30828, // 5.175
16060, 21126, 23805, 26531, 27718, 28646, 28781, 29550, 29644, 30505, 30016, 30754, 30448, 30750, 30535, 30725, 30944, 31043, // 5.7
14671, 20605, 23961, 26320, 27736, 28267, 28788, 30047, 29607, 30148, 30163, 30399, 30649, 30972, 30794, 30945, 30885, 30757, // 6.225
8605,  23442, 23827, 26323, 27103, 28172, 29011, 29244, 29892, 30213, 30071, 30200, 30326, 30782, 30771, 30719, 30814, 30970, // 6.75
//8605,  19491, 23827, 26323, 27103, 28172, 29011, 29244, 29892, 30213, 30071, 30200, 30326, 30782, 30771, 30719, 30814, 30970, // 6.75
2499,  25410, 21949, 25410, 27487, 28379, 28905, 29324, 30116, 30528, 30294, 30267, 30339, 30518, 30704, 31070, 30634, 30539, // 7.275
//2499,  23442, 21949, 25410, 27487, 28379, 28905, 29324, 30116, 30528, 30294, 30267, 30339, 30518, 30704, 31070, 30634, 30539, // 7.275
//2499,  21566, 21949, 25410, 27487, 28379, 28905, 29324, 30116, 30528, 30294, 30267, 30339, 30518, 30704, 31070, 30634, 30539, // 7.275
//2499,13521, 21949, 25410, 27487, 28379, 28905, 29324, 30116, 30528, 30294, 30267, 30339, 30518, 30704, 31070, 30634, 30539, // 7.275
308,   32768, 17037, 23775, 26876, 28266, 28765, 29224, 29869, 30083, 30550, 30262, 30148, 30225, 30507, 30733, 30662, 30876,  // 7.8
//308,   23981, 17037, 23775, 26876, 28266, 28765, 29224, 29869, 30083, 30550, 30262, 30148, 30225, 30507, 30733, 30662, 30876,  // 7.8
//308,   18324, 19084, 23775, 26876, 28266, 28765, 29224, 29869, 30083, 30550, 30262, 30148, 30225, 30507, 30733, 30662, 30876,  // 7.8
//308,   13491, 21969, 23775, 26876, 28266, 28765, 29224, 29869, 30083, 30550, 30262, 30148, 30225, 30507, 30733, 30662, 30876,  // 7.8
//264,  6460, 16899, 23775, 26876, 28266, 28765, 29224, 29869, 30083, 30550, 30262, 30148, 30225, 30507, 30733, 30662, 30876,  // 7.8
149,    1726, 31253, 17619, 26508, 27531, 28430, 29130, 30123, 29497, 29875, 30348, 30656, 30422, 30367, 31253, 30821, 30571,  // 8.325
//149,    1726, 31253, 20188, 26508, 27531, 28430, 29130, 30123, 29497, 29875, 30348, 30656, 30422, 30367, 31253, 30821, 30571,  // 8.325
//149,    1646, 26262, 20188, 26508, 27531, 28430, 29130, 30123, 29497, 29875, 30348, 30656, 30422, 30367, 31253, 30821, 30571,  // 8.325
//149,    1646, 20188, 20188, 26508, 27531, 28430, 29130, 30123, 29497, 29875, 30348, 30656, 30422, 30367, 31253, 30821, 30571,  // 8.325
//149,  1351, 11015, 18867, 24717, 27531, 28430, 29130, 30123, 29497, 29875, 30348, 30656, 30422, 30367, 31253, 30821, 30571,  // 8.325
83,      181, 10094, 15507, 20656, 24972, 27825, 29257, 29837, 29746, 30002, 30449, 29924, 30677, 30240, 30792, 31195, 30579,  // 8.85 
//83,      181, 10094, 18623, 20656, 24972, 27825, 29257, 29837, 29746, 30002, 30449, 29924, 30677, 30240, 30792, 31195, 30579,  // 8.85 
//83,      181,  7674, 18623, 20656, 24972, 27825, 29257, 29837, 29746, 30002, 30449, 29924, 30677, 30240, 30792, 31195, 30579,  // 8.85 
//83,    164,  4935, 14198, 20656, 24972, 27825, 29257, 29837, 29746, 30002, 30449, 29924, 30677, 30240, 30792, 31195, 30579,  // 8.85
46,      110,   920, 44307, 15512, 21534, 25278, 27661, 29158, 29857, 30579, 30358, 30405, 30720, 30545, 31164, 31029, 30716,  // 9.375
//46,      110,   920, 31164, 15512, 21534, 25278, 27661, 29158, 29857, 30579, 30358, 30405, 30720, 30545, 31164, 31029, 30716,  // 9.375
//46,      110,   920, 31164, 15512, 21534, 25278, 27661, 29158, 29857, 30579, 30358, 30405, 30720, 30545, 31164, 31029, 30716,  // 9.375
//46,      110,   920, 18602, 18602, 21534, 25278, 27661, 29158, 29857, 30579, 30358, 30405, 30720, 30545, 31164, 31029, 30716,  // 9.375
//43,     90,   750,  9680, 16597, 21534, 25278, 27661, 29158, 29857, 30579, 30358, 30405, 30720, 30545, 31164, 31029, 30716,  // 9.375
24,       61,   106,  5120, 18499, 19519, 22644, 25139, 27920, 29514, 29978, 30362, 30714, 30367, 30574, 30753, 30583, 30533,  // 9.9
//24,       61,   106,  4846, 18499, 19519, 22644, 25139, 27920, 29514, 29978, 30362, 30714, 30367, 30574, 30753, 30583, 30533,  // 9.9
//28,     61,    98,  3707, 12921, 18356, 22644, 25139, 27920, 29514, 29978, 30362, 30714, 30367, 30574, 30753, 30583, 30533,  // 9.9
12,       35,    68,   479, 78740, 14643, 19785, 22556, 25692, 26650, 29092, 30093, 29818, 30495, 30346, 30714, 30737, 30483,  // 10.425   
//12,       35,    68,   479, 44410, 14643, 19785, 22556, 25692, 26650, 29092, 30093, 29818, 30495, 30346, 30714, 30737, 30483,  // 10.425   
//12,       35,    68,   479, 24613, 17307, 19785, 22556, 25692, 26650, 29092, 30093, 29818, 30495, 30346, 30714, 30737, 30483,  // 10.425   
//12,       35,    68,   479, 17307, 17307, 19785, 22556, 25692, 26650, 29092, 30093, 29818, 30495, 30346, 30714, 30737, 30483,  // 10.425   
//13,     27,    61,   377,  8410, 15384, 19785, 22556, 25692, 26650, 29092, 30093, 29818, 30495, 30346, 30714, 30737, 30483,  // 10.425   
7,        19,    34,    69,  3642, 18396, 17075, 20323, 22903, 25062, 26872, 28145, 29880, 30241, 30527, 30392, 30546, 30532,  // 10.95
//7,        19,    34,    69,  3642, 17075, 17075, 20323, 22903, 25062, 26872, 28145, 29880, 30241, 30527, 30392, 30546, 30532,  // 10.95
//11,     22,    25,    59,  2935, 11899, 17075, 20323, 22903, 25062, 26872, 28145, 29880, 30241, 30527, 30392, 30546, 30532,  // 10.95
4,        10,    23,    37,   220, 176481, 12807, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679,  // 11.475
//4,        10,    23,    37,   220, 26979, 15718, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679,  // 11.475
//4,        10,    23,    37,   220, 16718, 16718, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679,  // 11.475
//5,       8,    19,    31,   193,  7410, 14148, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679,  // 11.475
4,        10,    23,    37,   220, 16718, 15718, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679   //  12 - 14
//4,        10,    23,    37,   220, 71037, 12807, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679,  // 11.475
//4,        10,    23,    37,   220, 26979, 16718, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679,  // 11.475
//5,       8,    19,    31,   193,  7410, 14148, 18215, 21255, 23445, 25335, 26852, 28343, 28721, 30278, 30354, 30123, 30679  // 
};
std::vector <double> invpdf = get_inverse_pdf(occ);


bool Py8PtGunV3::generatePartonsAndHadronize()
{

   fMasterGen->event.reset();

   for ( size_t i=0; i<fPartIDs.size(); i++ )
   {

      int particleID = fPartIDs[i]; // this is PDG - need to convert to Py8 ???

      double rand_sampler = rand() / double(RAND_MAX);
      double pt           = (fMaxPt-fMinPt) * randomEngine().flat() + fMinPt;
      double mass         = (fMaxMass-fMinMass) * randomEngine().flat() + fMinMass;
      double weight       = lookup_invpdf(mass, m_bins, pt, pT_bins, invpdf);
      while ( rand_sampler > weight ) {
         rand_sampler = rand() / double(RAND_MAX);
         pt           = (fMaxPt-fMinPt) * randomEngine().flat() + fMinPt;
         mass         = (fMaxMass-fMinMass) * randomEngine().flat() + fMinMass;
         weight       = lookup_invpdf(mass, m_bins, pt, pT_bins, invpdf);
      }
      // Calculate angles
      double phi = (fMaxPhi-fMinPhi) * randomEngine().flat() + fMinPhi;
      double eta = (fMaxEta-fMinEta) * randomEngine().flat() + fMinEta;
      double the = 2.*atan(exp(-eta));

      // Calculate momenta
      double pp = pt / sin(the); // sqrt( ee*ee - mass*mass );
      double ee = sqrt( pp*pp + mass*mass );

      double px = pt * cos(phi);
      double py = pt * sin(phi);
      double pz = pp * cos(the);

      if ( !((fMasterGen->particleData).isParticle( particleID )) )
      {
         particleID = std::fabs(particleID) ;
      }

      if( 1<= fabs(particleID) && fabs(particleID) <= 6) // quarks
        (fMasterGen->event).append( particleID, 23, 101, 0, px, py, pz, ee, mass );
      else if (fabs(particleID) == 21)                   // gluons
        (fMasterGen->event).append( 21, 23, 101, 102, px, py, pz, ee, mass );
      else                                               // other
        (fMasterGen->event).append( particleID, 1, 0, 0, px, py, pz, ee, mass );

      // Here also need to add anti-particle (if any)
      // otherwise just add a 2nd particle of the same type
      // (for example, gamma)
      if ( fAddAntiParticle )
      {
        if( 1 <= fabs(particleID) && fabs(particleID) <= 6){ // quarks
          (fMasterGen->event).append( -particleID, 23, 0, 101, -px, -py, -pz, ee, mass );
        }
        else if (fabs(particleID) == 21){                   // gluons
          (fMasterGen->event).append( 21, 23, 102, 101, -px, -py, -pz, ee, mass );
        }
        else if ( (fMasterGen->particleData).isParticle( -particleID ) ){
          (fMasterGen->event).append( -particleID, 1, 0, 0, -px, -py, -pz, ee, mass );
        }
        else {
          (fMasterGen->event).append( particleID, 1, 0, 0, -px, -py, -pz, ee, mass );
        }

      } // antiparticle

   } // fPartIDs

   if ( !fMasterGen->next() ) return false;

   event().reset(new HepMC::GenEvent);
   return toHepMC.fill_next_event( fMasterGen->event, event().get() );

} // generatePartonsAndHadronize()

const char* Py8PtGunV3::classname() const
{
   return "Py8PtGunV3";
}

typedef edm::GeneratorFilter<gen::Py8PtGunV3, gen::ExternalDecayDriver> Pythia8PtGunV3;

} // end namespace

using gen::Pythia8PtGunV3;
DEFINE_FWK_MODULE(Pythia8PtGunV3);
