#include "RawDataOverlayProcessor.h"
#include <iostream>
#include <string>
#include <algorithm>

#include "marlin/Processor.h"
#include "marlin/Global.h"

// include what you need from lcio
#include "IO/LCWriter.h"
#include "IO/LCReader.h"

#include <lcio.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <EVENT/LCEvent.h>
#include <EVENT/TrackerData.h>
#include <Exceptions.h>

#include <UTIL/LCTOOLS.h>
#include <UTIL/Operators.h>
#include <UTIL/BitSet32.h>

//LCCD
#include "lccd.h"
#include "lccd/DBInterface.hh"


// the common flag word definitions in MarlinTPC
#include "FlagwordDefinitions.h"

using namespace lcio;
using namespace marlin;
using namespace std;


namespace marlintpc
{

RawDataOverlayProcessor aRawDataOverlayProcessor;


RawDataOverlayProcessor::RawDataOverlayProcessor() : Processor("RawDataOverlayProcessor"),
                      				     _lcReader(0 ) ,
						     _overlayEvent(0 ) ,
						     _activeRunNumber(0 ) ,
						     _nRun(0 ) ,
						     _nEvt(0 ) ,
						     _nOverlayEvt( 0 ),
						     _events(0 ) {
  // the processor description
  _description = "RawDataOverlayProcessor is intended to ...";
  
    StringVec file ;
    file.push_back( "overlay.slcio" ) ;
    registerProcessorParameter( "InputFileName" , 
                                "Name of the lcio input file"  ,
                                _fileName ,
                                file
			      );

    registerInputCollection(LCIO::TRACKERRAWDATA,
                             "InputCollectionName" ,
                             "Name of the input collection"  ,
                             _inputColName ,
                             std::string("TPCRawData")
                            );

    registerOutputCollection(LCIO::TRACKERRAWDATA,
                             "OutputCollectionName" ,
                             "Name of the output collection"  ,
                             _outputColName ,
                             std::string("TPCRawData")
                            );

}

RawDataOverlayProcessor::~RawDataOverlayProcessor()
{
}

void RawDataOverlayProcessor::init()
{ 
  streamlog_out(DEBUG) << " init() called  " << endl;
  
  // you might want to list the used parameters by invoking
  printParameters();
  // opening background input
  _lcReader = LCFactory::getInstance()->createLCReader( LCReader::directAccess ) ;

  _lcReader->open( _fileName ) ;

  streamlog_out( DEBUG6 ) << "*** opened files for overlay - first file : " << _fileName[0]
	  << "    has " << _lcReader->getNumberOfEvents() << "  events. "  << std::endl ;


  if( _fileName.size() == 1 ) {

	  // if we have only one file, we can use direct access to random events
	  _lcReader->getEvents( _events ) ;

	  streamlog_out( DEBUG6 ) << "    will use direct access from " << _events.size() / 2 << "  events. "  << std::endl ;
  }
  //_nRun = 0 ;
  _nEvt = 0 ;
}


void RawDataOverlayProcessor::processRunHeader(LCRunHeader* run)
{
  // do something for each run
  //  _nRun++ ;
	run->parameters().setValue(_processorName + "_revision", "$Rev: $");

	for(ProcParamMap::iterator i = _map.begin(); i != _map.end(); i++) {
		if(!i->second->isOptional() || i->second->valueSet()) {
			run->parameters().setValue(_processorName + "_" + i->second->name(), i->second->value());
	        }
	}
} 


void RawDataOverlayProcessor::processEvent(LCEvent* evt)
{ 
  // do something for each event, usually access the input collection:

   if( _events.size() == 0 ) {   // in order to be reproduceable, we have to reset the file stream when using the skip mechanism

	   _lcReader->close() ;
	   _lcReader->open( _fileName ) ;
   }

   _overlayEvent = readNextEvent() ;
   ++_nOverlayEvt ;

   // merge event from storage with EVT

   LCCollection* outputCol = new LCCollectionVec(LCIO::TRACKERRAWDATA);

   // set the flags that cellID1 should be stored
   lcio::LCFlagImpl trkFlag(0) ;
   trkFlag.setBit(lcio::LCIO::TRAWBIT_ID1) ;
   outputCol->setFlag(trkFlag.getFlag());

   RawDataOverlayProcessor::_merge(_overlayEvent->getCollection(_inputColName), evt->getCollection(_inputColName), outputCol);
   
    _nEvt ++ ;

    evt->addCollection(outputCol, _outputColName);

}


void RawDataOverlayProcessor::check(LCEvent* evt)
{ 
  // check what is needed here 
}


void RawDataOverlayProcessor::end()
{ 
  // do something at the end of the processing, e.g. cleaning up

    streamlog_out( MESSAGE ) << " ------------------------------------------ "
                             << "   Overlay processor " << _nOverlayEvt << " background events on " << _nEvt << " physics events.\n"
                             << "      -> mean = " << double(_nOverlayEvt ) / double( _nEvt )
                             << "  ±  " <<   double(_nOverlayEvt ) / double( _nEvt ) / sqrt( _nEvt ) << "\n"
                             << std::endl ;
}

LCEvent* RawDataOverlayProcessor::readNextEvent() {

    LCEvent* overlayEvent = 0 ;

    // if we are reading from more than one file, we need to use the skipNEvents method,
    // otherwise we can directly access a random event

    int nEvts = _events.size() / 2 ;

    streamlog_out( DEBUG3 ) << "  --- Overlay::readNextEvent()  ---- will use "  << ( nEvts ?  " direct access " : " skipNEvents " )    << std::endl ;


    if( nEvts > 0 ){

	    int runNum = _events[ 2 * _nEvt     ] ;
	    int evtNum = _events[ 2 * _nEvt + 1 ] ;

	    streamlog_out( DEBUG3 ) << "   will read event " << evtNum << "  from  run " <<  runNum << std::endl ;

	    overlayEvent = _lcReader->readEvent(  runNum , evtNum , LCIO::UPDATE ) ;


	    if( !overlayEvent ) {

		    streamlog_out( ERROR ) << "   +++++++++ could not read event " << evtNum << "  from  run " <<  runNum << std::endl ;
		    return 0 ;
	    }

    }
    _activeRunNumber = overlayEvent->getRunNumber();

    return overlayEvent ;
}


inline long long cellID2long(int id0, int id1) { return ((long long) id0 << 32) | id1; }
inline int cellID0fromlongID(long long longID) { return (int) (longID >> 32); }
inline int cellID1fromlongID(long long longID) { return (int) (longID & 0x00000000ffffffff); }

void RawDataOverlayProcessor::_fillADC (LCCollection* col, map<long long, vector<short int> *> &ADCMap) 
{   
      TrackerRawDataImpl *rawdata;
      int nElements = col->getNumberOfElements();
      
      static const int kMaxTimeBin = 1024;

      for (int i=0; i<nElements; i++) {
        rawdata = dynamic_cast<TrackerRawDataImpl *> ( col->getElementAt(i) );
	if (!rawdata) continue;
	long long longID = cellID2long(rawdata->getCellID0(), rawdata->getCellID1());

	vector<short int> * adcvecp = ADCMap[longID];

	if (!adcvecp) {
	  adcvecp = new vector<short int>(kMaxTimeBin, 0.);
	  ADCMap[longID] = adcvecp;
	}

	vector<short int> &adcvec = *adcvecp;
	
	const EVENT::ShortVec & adcv = rawdata->getADCValues(); 	
	int nadc  = adcv.size();
	int tbin0 = rawdata->getTime();

	for (int iadc = 0; iadc < nadc; iadc++) {
	  int itime = tbin0 + iadc;
	  if (itime < 0 || itime >= kMaxTimeBin ) {
	  } else { 
	    adcvec[itime] = adcvec[itime] + adcv[iadc];
	  }
	}
      }

  }


void RawDataOverlayProcessor::_merge(LCCollection* src, LCCollection* dest, LCCollection* out) {
    int nElementsSrc = src->getNumberOfElements();
    int nElementsDest = dest->getNumberOfElements();
    const string destType = dest->getTypeName();
    
    // check if collections have the same type
    if (destType != src->getTypeName()) {
      streamlog_out( WARNING ) << "merge not possible, collections of different type" << endl;
      return;
    }
    
    streamlog_out( DEBUG4 ) << "merging collection of type: " << destType << " --- \n";
    
    if (destType == LCIO::TRACKERRAWDATA) {

      map<long long, vector<short int> *> ADCMap;

      streamlog_out( DEBUG4 ) << "merging ...  nElements(src) = "  << nElementsSrc  << endl;
      streamlog_out( DEBUG4 ) << "merging ...  nElements(dest) = " << nElementsDest << endl;
      
      cerr << "merging ...  nElements(src) = "  << nElementsSrc  << endl;
      cerr << "merging ...  nElements(dest) = " << nElementsDest << endl;

      _fillADC(dest, ADCMap);
      _fillADC( src, ADCMap);

      static const short int kThreshold = 1; 
      static const short int kPedADC = 100; 
      static const short int kMaxADC = 1024 - kPedADC; 
#if 0
      int nChannels = 0;
      cerr << "nChannels: " << ADCMap.size() << endl;
#endif
      for (auto itr = ADCMap.begin(); itr != ADCMap.end(); itr++) {
	long long longID = itr->first;
	int id0 = cellID0fromlongID(longID);
	int id1 = cellID1fromlongID(longID);
#if 0
      nChannels++;
      cerr << "-------------------------------------------- "<< endl;
      cerr << "nChannels: " << nChannels << endl;
      cerr << "longID " << (void *) longID << endl;
      cerr << "id0: " << id0 << " id1: " << id1 << endl;
#endif
        vector<short int> *adcvecp = itr->second;
	vector<short int> &adcvec  = *adcvecp;
#if 0
      cerr << "itr->second: " << (void *) itr->second << endl;
      cerr << "-------------------------------------------- "<< endl;
#endif
	EVENT::ShortVec adcvalues; 	

	int nadc = adcvec.size();
#if 0
      cerr << "nadc: " << nadc << endl;
#endif
	int itimebin0 = -9999;
	bool isfound = false;

	for (int itime = 0; itime < nadc; itime++) { //loop 1
	  short int adc = adcvec[itime];
	  
	  if (!isfound && adc < kThreshold) { 

	      continue;
	  } else if (!isfound && adc >= kThreshold) { 

	      isfound = true;
	      itimebin0 = itime;
	      adcvalues.clear(); 
	      adc = min(adc, kMaxADC);
 	      adcvalues.push_back(adc);
	  } else if (isfound && adc >= kThreshold) { 
	      adc = min(adc, kMaxADC);
 	      adcvalues.push_back(adc);
	  } else if (isfound && adc < kThreshold) { 

	      isfound = false;
	      TrackerRawDataImpl * rawdata = new TrackerRawDataImpl;
	      rawdata->setADCValues(adcvalues);
	      rawdata->setTime(itimebin0);
	      rawdata->setCellID0(id0);
	      rawdata->setCellID1(id1);
              out->addElement(rawdata);
#if 0
	      cerr << "id0: " << rawdata->getCellID0() << 
		     " id1: " << rawdata->getCellID1() ;
	      cerr << " itimebin0: " << itimebin0 << " " ;
	      cerr << "adc: " ;
	      for (unsigned int j = 0; j < adcvalues.size(); j++) {
		 cerr << adcvalues[j] << ", " ;
	      }
	      cerr << endl;
#endif
	  }
	}
#if 0
        cerr << "-------------------------------------------- "<< endl;
        cerr << "adcvecp:  " << (void *) adcvecp << endl;
#endif
	delete adcvecp;
	ADCMap[longID] = 0;

      } 
      return;
    }

}


} // close namespace brackets
