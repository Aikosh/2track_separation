#ifndef RAWDATAOVERLAYPROCESSOR_H
#define RAWDATAOVERLAYPROCESSOR_H

#include <marlin/Processor.h>
#include <marlin/EventModifier.h>
#include <lcio.h>
#include <string>
#include <map>


namespace marlintpc
{
///  Short description of processor.

/**
 * Lengthy and detailed description of processor ...
 * @author Aiko. Shoji, Iwate University, Keisuke. Fujii, KEK, Oliver. Schaefer, DESY
 *
 * \par Credits:
 * The processor skeleton was generated by the script \c createProcessor.py
 *
 */

class RawDataOverlayProcessor : public marlin::Processor//, public marlin::EventModifier 
{
 public:
  
  virtual Processor*  newProcessor() { return new RawDataOverlayProcessor ; }
  
  RawDataOverlayProcessor();

  virtual ~RawDataOverlayProcessor();
 
  //virtual const std::string & name() const { return Processor::name() ; }

  //virtual void modifyEvent( EVENT::LCEvent * evt ) ;
 
  virtual void init();
  
  virtual void processRunHeader(lcio::LCRunHeader* run);
  
  virtual void processEvent(lcio::LCEvent* evt);
  
  virtual void check(lcio::LCEvent* evt);
  
  virtual void end();
  
  
 protected:
  /* Helper method.*/
    LCEvent* readNextEvent() ;

  /* the place for protected and private member data and functions */
  std::string _inputColName; ///< Name of the input collection

  StringVec _fileName;

  std::string _outputColName;

  std::map< std::pair<int,int>, std::vector<short int>* > _id2adcMap;

  IO::LCReader* _lcReader ;
  EVENT::LCEvent* _overlayEvent ;
  int _activeRunNumber;
  int _nRun ;
  int _nEvt ;
  int _nOverlayEvt ;
  IntVec _events ;

  static void _merge(EVENT::LCEvent* srcEvent, EVENT::LCEvent* destEvent, EVENT::LCCollection* out);
  static void _merge(EVENT::LCCollection* src, EVENT::LCCollection* dest, EVENT::LCCollection* out);
  static void _fillADC(EVENT::LCCollection* col, std::map<long long, std::vector<short int> *> &ADCMap);

};
} // namespace marlintpc
#endif // RAWDATAOVERLAYPROCESSOR_H
