// -*- C++ -*-

#ifndef KVC_TRACK_INFO_HH
#define KVC_TRACK_INFO_HH

#include <G4VUserTrackInformation.hh>

//_____________________________________________________________________________
class KVC_TrackInfo : public G4VUserTrackInformation
{
public:
  KVC_TrackInfo(G4bool is_from_quartz) : m_is_from_quartz(is_from_quartz) {}
  ~KVC_TrackInfo() override {}

  G4bool IsFromQuartz() const { return m_is_from_quartz; }

private:
  G4bool m_is_from_quartz;
};

#endif
