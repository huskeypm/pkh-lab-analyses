"""For time stamping and time deltas"""

#Standard library
from datetime import datetime

#Standard time format
TIMEFORMAT="%a %d-%b-%Y %H:%M:%S.%f %Z"

def format_delta(tdelta):
  "Return a string in standard format for a time delta"
  hours,rem=divmod(tdelta.seconds,3600)
  mins,sec=divmod(rem,60)
  return "%dD %02d:%02d:%02d.%06d (%f sec)"%(tdelta.days,hours,mins,sec,tdelta.microseconds,tdelta.total_seconds())

class Timer(object):
  """A simple timer"""
  def __init__(self):
    self.start=datetime.now()
  @classmethod
  def start(cls):
    "Convenience method for starting a timer."
    return cls()
  def stop(self):
    "Stop the timer. Store the timedelta and a string of the elapsed time. Return the string."
    self.end=datetime.now()
    self.delta = self.end - self.start
    self.delta_str = format_delta(self.delta)
    return self.delta_str
  def split(self,when=None):
    "Store the timedelta and a string of the current elapsed time. Return the string."
    if when is None:
      when = datetime.now()
    self.lap = when - self.start
    self.lap_str = format_delta(self.lap)
    return self.lap_str

def timestamp(dt=None):
  if dt is None:
    dt=datetime.now()
  return dt.strftime(TIMEFORMAT)