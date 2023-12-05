#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/Toolkit.h"
// #include <math.h>
// #include <algorithm>

  
using namespace l1thgcfirmware;

#define MEMBERS &HGCalCluster::clock_ , \
                &HGCalCluster::index_ , \
                &HGCalCluster::n_tc_, \
                &HGCalCluster::e_, \
                &HGCalCluster::e_em_, \
                &HGCalCluster::e_em_core_, \
                &HGCalCluster::e_h_early_, \
                &HGCalCluster::w_, \
                &HGCalCluster::n_tc_w_, \
                &HGCalCluster::w2_, \
                &HGCalCluster::wz_, \
                &HGCalCluster::wphi_, \
                &HGCalCluster::wroz_, \
                &HGCalCluster::wz2_, \
                &HGCalCluster::wphi2_, \
                &HGCalCluster::wroz2_, \
                &HGCalCluster::layerbits_, \
                &HGCalCluster::sat_tc_, \
                &HGCalCluster::shapeq_, \
                &HGCalCluster::sortKey_, \
                &HGCalCluster::sortKey2_, \
                &HGCalCluster::L_, \
                &HGCalCluster::R_, \
                &HGCalCluster::X_, \
                &HGCalCluster::lastFrame_ , \
                &HGCalCluster::dataValid_
    // &HGCalCluster::weta_, 
    // &HGCalCluster::weta2_, 

bool HGCalCluster::operator==(const HGCalCluster& other) const { 
  return ::check( *this , other , MEMBERS );
}

std::ostream& operator<< ( std::ostream& aStr , const HGCalCluster& aCell )
{
  return ::print( aStr , aCell ,  MEMBERS );
}

HGCalClusterSAPtrCollection HGCalCluster::ReadDebugFile( const std::string& aFilename )
{
  HGCalClusterSAPtrCollection lRet;
  ::ReadDebugFile( aFilename , lRet , MEMBERS );
  return lRet;
}

#undef MEMBERS







const HGCalCluster& HGCalCluster::operator+=(const HGCalCluster& c) {

  n_tc_ += c.n_tc_;
  e_ += c.e_;
  e_em_ += c.e_em_;
  e_em_core_ += c.e_em_core_;
  e_h_early_ += c.e_h_early_;
  w_ += c.w_;
  // w2_ += c.w2_;
  // wz_ += c.wz_;
  // weta_ += c.weta_;
  // wphi_ += c.wphi_;
  // wroz_ += c.wroz_;
  // wz2_ += c.wz2_;
  // weta2_ += c.weta2_;
  // wphi2_ += c.wphi2_;
  // wroz2_ += c.wroz2_;
  // n_tc_w_ += c.n_tc_w_;
  layerbits_ |= c.layerbits_;
  sat_tc_ |= c.sat_tc_;
  shapeq_ |= c.shapeq_;

  // // if ( w_ <= 52438 and shapeq_ == 1 and c.shapeq_ == 1 ) { // Magic numbers
    // // shapeq_ = 1;
    // w_ = ...
    w2_ += c.w2_;
    wz_ += c.wz_;
    // weta_ += c.weta_;
    wphi_ += c.wphi_;
    wroz_ += c.wroz_;
    wz2_ += c.wz2_;
    // weta2_ += c.weta2_;
    wphi2_ += c.wphi2_;
    wroz2_ += c.wroz2_;
    n_tc_w_ += c.n_tc_w_;      
  // // } else if ( w_ <= c.w_ ) {
    // // shapeq_ = 0;
    // // w_ = c.w_;
    // // w2_ = c.w2_;
    // // wz_ = c.wz_;
    // // // weta_ = c.weta_;
    // // wphi_ = c.wphi_;
    // // wroz_ = c.wroz_;
    // // wz2_ = c.wz2_;
    // // // weta2_ = c.weta2_;
    // // wphi2_ = c.wphi2_;
    // // wroz2_ = c.wroz2_;
    // // n_tc_w_ = c.n_tc_w_;
  // // } else {
    // // shapeq_ = 0;
  // // }

  return *this;
}


void HGCalCluster::saturate()
{
  if( n_tc_.value_ > 0x3FF ) n_tc_.value_ = 0x3FF;
  if( e_.value_ > 0x3FFFFF ) e_.value_ = 0x3FFFFF;
  if( e_em_.value_ > 0x3FFFFF ) e_em_.value_ = 0x3FFFFF;
  if( e_em_core_.value_ > 0x3FFFFF ) e_em_core_.value_ = 0x3FFFFF;
  if( e_h_early_.value_ > 0x3FFFFF ) e_h_early_.value_ = 0x3FFFFF;
  if( w_.value_ > 0xFFFF ) w_.value_ = 0xFFFF;
  if( n_tc_w_.value_ > 0x3FF ) n_tc_w_.value_ = 0x3FF;
  if( w2_.value_ > 0xFFFFFFFF ) w2_.value_ = 0xFFFFFFFF;
  if( wz_.value_ > 0x1FFFFFFF ) wz_.value_ = 0x1FFFFFFF;
  if( wphi_.value_ > 0xFFFFFFF ) wphi_.value_ = 0xFFFFFFF;
  if( wroz_.value_ > 0x1FFFFFFF ) wroz_.value_ = 0x1FFFFFFF;
  if( wz2_.value_ > 0x3FFFFFFFFFF ) wz2_.value_ = 0x3FFFFFFFFFF;
  if( wphi2_.value_ > 0xFFFFFFFFFF ) wphi2_.value_ = 0xFFFFFFFFFF;
  if( wroz2_.value_ > 0x3FFFFFFFFFF ) wroz2_.value_ = 0x3FFFFFFFFFF;      
}

