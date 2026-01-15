#include "EventAction.hh"
#include "AnaManager.hh"

namespace
{
  auto& gAnaMan = AnaManager::GetInstance();
}

EventAction::EventAction() : G4UserEventAction(), fNCerGen(0) {
}

EventAction::~EventAction() {
}

void EventAction::BeginOfEventAction(const G4Event* anEvent) {
  gAnaMan.BeginOfEventAction(anEvent);
  fNCerGen = 0; //追加　チェレンコフ光生成数の初期化
}

void EventAction::EndOfEventAction(const G4Event* anEvent) {
  G4int eventID = anEvent->GetEventID();

  //追加
  gAnaMan.EndOfEventAction(anEvent); //AnaManagerにチェレンコフ光生成数を渡す
  gAnaMan.SetCerGen(fNCerGen);  // Anamanager側でイベントデータを保存
  //gAnaMan.FillTree();
  //追加終

  if (eventID % 100 == 0) {
    G4cout << "   Event number = " << eventID << G4endl;
  }


}