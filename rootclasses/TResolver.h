#ifndef TResolver_h
#define TResolver_h

#include <TObjArray.h>

class TObject;
class TChannel;
class TReconstruction;

class TResolver : public TObject
{

private:
  TObjArray fChannelLists;	// array of collections of TChannels

protected:
  void FindAll(EChannelKind kind); //calls the right FindAll
  void FindAll(EChannelKind kind,EChannelKind from1,EChannelKind from2);
  // finds all pairs from from1 & from2 & adds them to kind
  void FindAll(EChannelKind kind,EChannelKind from1,EChannelKind from2
                                ,EChannelKind from3);

public:
  TResolver(const TReconstruction *event); // reads gammas from event
  ~TResolver();

  void FindAll(); //finds all pairs and triplets that might make a channel

  void Display() const; //prints out the decay channel
  Int_t GetNumInChannel(EChannelKind kind)
          {return GetChannelList(kind)->GetSize();}
  const TCollection *GetChannelList(EChannelKind kind)
	  {return (TCollection*)fChannelLists[kind];}

  ClassDef(TResolver,1) // Channel resolving class
};

#endif
