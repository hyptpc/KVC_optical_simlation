#ifndef KVC_TRACK_INFO_HH
#define KVC_TRACK_INFO_HH

#include "G4VUserTrackInformation.hh"

class KVC_TrackInfo : public G4VUserTrackInformation {
public:
    KVC_TrackInfo(bool isFromQuartz) : fIsFromQuartz(isFromQuartz) {}
    virtual ~KVC_TrackInfo() {}

    bool IsFromQuartz() const { return fIsFromQuartz; }

private:
    bool fIsFromQuartz;
};

#endif
