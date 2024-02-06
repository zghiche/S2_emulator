#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLinkUnpacking_SA.h"
#include <unistd.h> 
#include <algorithm>

using namespace std;
using namespace l1thgcfirmware;

HGCalLinkUnpacking::HGCalLinkUnpacking(const ClusterAlgoConfig& config) : config_(config) {}

void HGCalLinkUnpacking::runLinkUnpacking(const HGCalLinkTriggerCellSAPtrCollection& linksIn,
                                    HGCalTriggerCellSAPtrCollection& triggerCellsOut) const {
    HGCalLinkTriggerCellSAPtrCollection triggerCellsWork, towersWork;
    unpackLinks( linksIn, triggerCellsWork, towersWork);
    triggerCellsOut = triggerCellDistribution( triggerCellsWork );
    unpackTriggerCells( triggerCellsOut );
}

void HGCalLinkUnpacking::unpackLinks( const HGCalLinkTriggerCellSAPtrCollection& LinksIn , HGCalLinkTriggerCellSAPtrCollection& TriggerCells , HGCalLinkTriggerCellSAPtrCollection& Towers ) const
{
  const unsigned int stepLatency = config_.getStepLatency( UnpackLinks );

  for ( unsigned int i(0) ; i!=LinksIn.size() ; i+=84 ) {
    for ( int j(0) ; j!=2 ; ++j ) {
      for ( int k(0) ; k!=84 ; ++k ) {
        const auto& in = LinksIn.at( i + k );
        uint64_t val = in->data_.value_;
        uint64_t r_over_z = in->r_over_z_.value_;
        uint64_t phi = in->phi_.value_;
        for( int l(0) ; l!=3 ; ++l )
        {
          auto tc = make_unique< HGCalLinkTriggerCell >();
          tc->clock_ = in->clock_+j + stepLatency;
          tc->index_ = (3*k)+l;
          tc->data_ = val & 0x7FFF;
          val >>= 15;
          tc->r_over_z_ = r_over_z & 0x7FFF;
          r_over_z >>= 15;
          tc->phi_ = phi & 0x7FFF;
          phi >>= 15;
          tc->lastFrame_ = in->lastFrame_ and j==1;
          tc->dataValid_ = true;
          TriggerCells.emplace_back( move(tc) );
        }
      }
    }
  }  
}
// =======================================================================================================================================================

// =======================================================================================================================================================
HGCalTriggerCellSAPtrCollection HGCalLinkUnpacking::triggerCellDistribution( const HGCalLinkTriggerCellSAPtrCollection& TriggerCellsIn ) const
{
  const unsigned int stepLatency = config_.getStepLatency( TriggerCellDistribution );
  
  const size_t Nframes = 216;
  const size_t Nchannels = 84 * 3;
  
  HGCalTriggerCellSAPtrCollection TriggerCellsOut;
  for ( unsigned int frame = 0; frame != Nframes; ++frame ) {
    for ( unsigned int iColumn = 0; iColumn != config_.cColumns(); ++iColumn ) {
      auto& lut_out = config_.TriggerCellDistributionLUT( ( Nframes*iColumn ) + frame );
      int valid = ( lut_out >> 39 ) & 0x1;
   
      if( valid )
      {   
        int R_over_Z = ( lut_out >> 0 )  & 0xFFF;
        int Phi      = ( lut_out >> 12 ) & 0xFFF;
        int Layer    = ( lut_out >> 24 ) & 0x3F;
        int index    = ( lut_out >> 30 ) & 0x1FF;

        // cout << R_over_Z << ", " << Phi << endl;
        auto& in = TriggerCellsIn.at( ( Nchannels * frame ) + index );       
        TriggerCellsOut.emplace_back(
          make_unique< HGCalTriggerCell >( 
            (frame==Nframes-1),
            true,
            in->r_over_z_.value_,
            in->phi_.value_,
            Layer,
            in->data_.value_
          )
        );
        auto& tc = TriggerCellsOut.back();
        tc->setClock(in->clock_ + stepLatency);
        tc->setIndex( iColumn );
      }
    }
  }

  return TriggerCellsOut;
}
// =======================================================================================================================================================

// =======================================================================================================================================================
void HGCalLinkUnpacking::unpackTriggerCells(HGCalTriggerCellSAPtrCollection& triggerCells ) const
{
  const unsigned int stepLatency = config_.getStepLatency( UnpackTriggerCells );
  int mantissaBits_ = 3;

  for ( auto& tc : triggerCells ) {
    // Energy
    uint32_t energyMantissa = (tc->energy() >> 4) & 0x7;
    uint32_t energyExponent = tc->energy() & 0xF;
    tc->setEnergy(((1ULL << mantissaBits_) | energyMantissa) << (energyExponent - 1));
    
    // Phi
    uint32_t phiMantissa = (tc->phi() >> 4) & 0x7;
    uint32_t phiExponent = tc->phi() & 0xF;
    tc->setPhi(((1ULL << mantissaBits_) | phiMantissa) << (phiExponent - 1));
     
    // R_over_Z
    uint32_t rOverZMantissa = (tc->rOverZ() >> 4) & 0x7;
    uint32_t rOverZExponent = tc->rOverZ() & 0xF;
    tc->setROverZ(((1ULL << mantissaBits_) | rOverZMantissa) << (rOverZExponent - 1));

    // if (tc->energy() > 0) {
    //   std::cout << tc->energy()  << ", " << tc->phi() << ", " << tc->rOverZ() << std::endl;
    // }
    tc->addLatency( stepLatency );
  }

  sort(triggerCells.begin(), triggerCells.end(), [](auto& a, auto& b) {
         return a->rOverZ() > b->rOverZ();
       });
  // for ( auto& tc : triggerCells ) {
  //  if (tc->energy() != 0) {
  //    std::cout << tc->rOverZ() << ", " << tc->phi() << std::endl;
  // }}

}
