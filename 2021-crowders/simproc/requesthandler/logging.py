"""For creating detailed logs

Each log entry has the following:
- timestamp: a timestamp string
- level: a logging level as defined by the python logging module
- message: a text string
- parameters: a dictionary of values to report

Log entries can also request initializing timers and reporting their current elapsed value.

Not implemented:
- entries: a list of sub-entries
"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import collections
import logging
from datetime import datetime
import sys

#Site packages
#(implement a separate yaml instance from the one in yaml_manager, to avoid conflicts)
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap
yaml=YAML(typ="rt",pure=True)
yaml.default_flow_style = False

#This package
from . import filepath
from . import timing
from . import locators
from . import yaml_manager
from . import schema

#Constants
MAX_BUFFER_MESSAGES = 100
logging.TIMING = 25 #This puts it between INFO and WARNING

#Classes for logging

class LogRecord(object):
  def __init__(self,name,level,msg,args=None,**kwargs):
    #Get creation time
    self.created = datetime.now()
    #Usual parameters, like the stdlib version
    self.name = name
    self.msg = msg
    self.args = args
    self.levelname = logging.getLevelName(level)
    self.levelno = level
    #Store the additional parameters
    self.parameters=kwargs
  @property
  def message(self):
    msg = str(self.msg)
    if self.args is not None:
      msg = msg % self.args
    return msg

class Logger(logging.Logger):
  def __init__(self, name, level=logging.NOTSET):
    #Initialization from base class
    super(Logger, self).__init__(name,level)
    #Initialize dictionary of timers
    self.timers={}
  # def wrapFindCaller(self,stack_info=False,stacklevel=1,exc_info=None):
  #   #Copied from the python standard library source, in _log
  #   #Won't work yet
  #   sinfo = None
  #   if _srcfile:
  #     #IronPython doesn't track Python frames, so findCaller raises an
  #     #exception on some versions of IronPython. We trap it here so that
  #     #IronPython can use logging.
  #     try:
  #       fn, lno, func, sinfo = self.findCaller(stack_info, stacklevel)
  #     except ValueError: # pragma: no cover
  #       fn, lno, func = "(unknown file)", 0, "(unknown function)"
  #   else: # pragma: no cover
  #     fn, lno, func = "(unknown file)", 0, "(unknown function)"
  #   if exc_info:
  #     if isinstance(exc_info, BaseException):
  #       exc_info = (type(exc_info), exc_info, exc_info.__traceback__)
  #     elif not isinstance(exc_info, tuple):
  #       exc_info = sys.exc_info()
  #   return fn,lno,exc_info,func,extra,sinfo
  def _log(self, level, msg, args=None, **kwargs): #TODO: arguments would need to be modified to let wrapFindCaller work
    record = LogRecord(self.name, level, msg, args, **kwargs)
    self.handle(record)
  def startTimer(self,timername,**kwargs):
    self.timers[timername]=timing.Timer()
    self.log(logging.TIMING,"Starting new timer.",timername=timername,**kwargs)
  def splitTimer(self,timername,**kwargs):
    delta_str=self.timers[timername].split()
    delta_totalsec=self.timers[timername].lap.total_seconds()
    self.log(logging.TIMING,"Reporting elapsed time on timer.",
            timer_name=timername,
            elapsed_sec=delta_totalsec,
            elapsed=delta_str, **kwargs)
  def stopTimer(self,timername,**kwargs):
    delta_str=self.timers[timername].stop()
    delta_totalsec=self.timers[timername].delta.total_seconds()
    self.log(logging.TIMING,"Stopping timer.",
            timer_name=timername,
            elapsed_sec=delta_totalsec,
            elapsed=delta_str,**kwargs)

class RootLogger(Logger):
  def __init__(self,level):
    #Create a logger named "root"
    Logger.__init__(self,"root",level)
    #Set up a BufferHandler for messages while other handlers are being set up
    self.bufferhandler=BufferHandler()
    self.bufferhandler.name="Root Buffer Handler"
    self.addHandler(self.bufferhandler,skipbuffer=True)
  def addHandler(self, hdlr, skipbuffer=False):
    #Use method from base class
    super(Logger, self).addHandler(hdlr)
    if not skipbuffer:
      #Send the buffered messages to the new handler
      self.bufferhandler.dump_to(hdlr)
      #Note how many messsages were sent to the new handler
      self.info("Adding new handler, and sending it previous messages.",
                handler_type=type(hdlr).__name__,
                handler_name=hdlr.name,
                num_messages_output=self.bufferhandler.num_records,
                missing_messages=self.bufferhandler.num_deleted)
  def bumpLevel(self,level):
    """Make sure the logger is enabled for this level"""
    if self.level > logging._checkLevel(level):
      self.setLevel(level)

class YAMLStreamHandler(logging.StreamHandler):
  """For output of logs to a YAML stream"""
  def __init__(self,stream=None,use_stdout=False):
    if stream is None and use_stdout:
      stream=sys.stdout
    #Initialization from base class
    super(YAMLStreamHandler, self).__init__(stream)
  def emit(self,record):
    try:
      msg=self.format(record)
      stream=self.stream
      yaml.dump(msg,stream)
      self.flush()
    except RecursionError:
      raise
    except Exception:
      self.handleError(record)

class YAMLFileHandler(logging.FileHandler):
  #We don't need to override init, even though the base class calls StreamHandler instead of YAMLStreamHandler, because YAMLStreamHandler doesn't override the method called in init.
  #The same goes for close
  #We do need to override emit.
  def emit(self,record):
    #Check if opening the stream was delayed
    if self.stream is None:
      self.stream = self._open()
    #Call the appropriate stream handler emit
    YAMLStreamHandler.emit(self, record)

class YAMLFormatter(object):
  """Intended only for use with the YAMLStreamHandler"""
  def format(self,record):
    od=CommentedMap()
    od['timestamp']=timing.timestamp(record.created)
    od['area']=record.name
    od['level']=record.levelname
    od['message']=record.message
    if getattr(record,'parameters',None):
      od['parameters']=record.parameters
    return [od]

class BufferHandler(logging.Handler):
  """A handler that just keeps a limited buffer of log records
  
  Note that no formatter is needed, because the records are stored directly."""
  def __init__(self,level=logging.NOTSET):
    #Initialization from base class
    super(BufferHandler, self).__init__(level)
    #Initalize the buffer
    self.buffer=collections.deque(maxlen=MAX_BUFFER_MESSAGES)
    #Initialize record count
    self.total_records=0
  def emit(self,record):
    self.buffer.append(record)
    self.total_records+=1
  @property
  def num_records(self):
    return len(self.buffer)
  @property
  def num_deleted(self):
    return self.total_records-len(self.buffer)
  def dump_to(self,other):
    """Ask another handler to handle all the buffered records"""
    for record in self.buffer:
      other.handle(record)

#Introduce the "timing" log level
logging.addLevelName(logging.TIMING,"TIMING")
#Tell the stdlib logging module about these new classes
# logging.setLogRecordFactory(LogRecord)
# logging.setLoggerClass(Logger)
#Set up the root logger
root=RootLogger(logging.WARNING)
Logger.root=root
#Set up the Manager
Logger.manager = logging.Manager(Logger.root)
Logger.manager.setLoggerClass(Logger)

def getLogger(name=None):
  #Redo the module-level function from stdlib
  if name:
    return Logger.manager.getLogger(name)
  else:
    return root

#Functions for configuring

def find_unique_id(stem,logdir,ext,num_digits,sepchar):
  """Create a file name that does not already exist.

  Arguments:

    - stem = log filename stem (the part before the unique ID and extension), as string
    - logdir = path to the directory to contain the log file. Absolute path. String form acceptable.
    - ext = extension to use for teh filename
    - num_digits = number of digits to use in the unique ID
    - sepchar = character to place between the stem and the unique ID

  Return: the filename calculated"""
  logdir_path=filepath.Path(logdir, isFile=False)
  assert logdir_path.is_dir(), "logdir must be a directory"
  existing_files=[c.name for c in logdir_path.iterdir()]
  id_tmpl="{0:0%dd}"%num_digits
  fname_tmpl=stem+sepchar+id_tmpl+ext
  unid=1
  trial_fname=fname_tmpl.format(unid)
  while trial_fname in existing_files:
    unid+=1
    trial_fname=fname_tmpl.format(unid)
  return trial_fname

_srcfile = find_unique_id.__code__.co_filename #So that findCaller will work, if needed

def configure_logfile(level="TIMING",stem="simproc",logdir_rel="logs",logdir_abs=None,ext=".log.yaml",num_digits=3,sepchar="."):
  """Set up a log file for the root logger, which will flow down to all children by default

  Arguments:

    - level = level below which to suppress messages to this file, as a string or integer
    - stem = log filename stem, as string
    - logdir_rel = path to the directory to contain the log file, taken relative to locators.DATAFOLDER **at the time of the call**. String form acceptable.
    - ext, num_digits, sepchar = additional arguments to ``find_unique_id``

  No return value."""
  #Get the root logger
  global root
  #Lower the level of the root logger if necessary
  root.bumpLevel(level)
  #Get the absolute path to the log directory
  global logdir
  if logdir_abs is None:
    logdir = locators.DATAFOLDER / filepath.Path(logdir_rel, isFile = False)
  else:
    logdir = filepath.Path(logdir_abs, isFile=False)
  #Make sure the directory exists
  logdir.assure_dir()
  #Get the filename and path for the log file
  logfile=find_unique_id(stem,logdir,ext,num_digits,sepchar)
  logfpath = str(logdir / logfile)
  #Create a handler to output to this file
  handler=YAMLFileHandler(logfpath,mode='w')
  handler.name=logfile
  #Set the level for the handler
  handler.setLevel(level)
  #Set the formatter for the handler
  formatter=YAMLFormatter()
  handler.setFormatter(formatter)
  #Add the handler to the root logger
  root.addHandler(handler)
  return

def configure_stdout(level="TIMING"):
  """Set up logs to sys.stdout"""
  #Get the root logger
  global root
  #Lower the level of the root logger if necessary
  root.bumpLevel(level)
  #Create the handler
  handler=YAMLStreamHandler(use_stdout=True)
  handler.setLevel(level)
  #Set the formatter for the handler
  formatter=YAMLFormatter()
  handler.setFormatter(formatter)
  #Add the handler to the root logger
  root.addHandler(handler)
  return

def configure_stderr(level="TIMING"):
  """Set up logs to sys.stderr"""
  #Get the root logger
  global root
  #Lower the level of the root logger if necessary
  root.bumpLevel(level)
  #Create the handler
  handler=YAMLStreamHandler()
  handler.setLevel(level)
  #Set the formatter for the handler
  formatter=YAMLFormatter()
  handler.setFormatter(formatter)
  #Add the handler to the root logger
  root.addHandler(handler)
  return

#Class for configuring logging from a yaml file

#Validation schema for ConfigLogging
_ConfigLogging_props_schema_yaml="""#ConfigLogging
level: {type: string}
stem: {type: string}
logdir_rel: {type: string}
ext: {type: string}
num_digits: {type: integer}
sepchar: {type: string}
"""

class ConfigLogging(schema.SelfValidating):
  """Class to configure logging parameters from within a yaml file

  This isn't really a class. It's a function.
  But when loading from yaml, there isn't a way to call a function directly.
  You can only load classes.
  Hence, this is an object that simply calls another function when initialized,
  and then does nothing else ever.
  
  This is also not a request: its action is taken at initialization,
  not when requests are run.

  Initialization arguments:
    - level: minimum level for events to log, as a string
    - stem = log filename stem, as string
    - logdir_rel = path to the directory to contain the log file, taken relative to locators.DATAFOLDER **at the time of the call**. String form acceptable.
    - ext, num_digits, sepchar = additional arguments to ``find_unique_id``
  """
  _validation_schema=schema.SelfValidating.update_schema(_ConfigLogging_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(ConfigLogging, self).__init__(**kwargs)
    #Call the configuration function
    configure_logfile(**kwargs)
 
#Register for loading from yaml
yaml_manager.register_classes([ConfigLogging])
