/******************************************************************************
 * @file EventDescriptor.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  EventDescriptor.h
 * ****************************************************************************/

// Decay vertex header
#include "EventDescriptor.h"

namespace AuxEvent
{
  EventDescriptor::EventDescriptor()
  {}
  EventDescriptor::~EventDescriptor()
  {}

  void EventDescriptor::Initialize(int run, int subrun, int event)
  {
    run = run;
    subrun = subrun;
    event = event;
  } // END function Initialize

} // END namespace EventDescriptor 
