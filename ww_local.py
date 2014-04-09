import ctypes as ct
import numpy as np
import datetime, os

#libddww = ct.cdll.LoadLibrary('/afs/ipp-garching.mpg.de/aug/ads/lib64/@sys/libddww8.so')
libddww = ct.cdll.LoadLibrary('/afs/ipp-garching.mpg.de/aug/ads/lib64/@sys/libddww8.so.8.1.20110418')

class ObjectInfo:
  status = False
  def __init__( self , name=None,typ=None,format=None,length=None,items=None,indices=None ):
    self.status = False
    if name != None and typ != None and format!=None and length!=None and items !=None and indices!=None:
      self.Set(name,typ,format,length,items,indices)

  def __del__( self ):
    self.Unset()

  def Set( self , name , typ , format , length , items , indices ):
    self.Unset()
    vars(self)['name'] = name
    vars(self)['typ'] = typ
    vars(self)['format'] = format
    vars(self)['length'] = length
    vars(self)['items'] = items
    vars(self)['indices'] = indices
    self.status = True

  def Unset( self ):
    if self.status:
      del self.name
      del self.typ
      del self.format
      del self.length
      del self.items
      del self.indices
      self.status = False

ctList = [ct.c_bool, ct.c_char, ct.c_wchar, ct.c_byte,\
              ct.c_ubyte, ct.c_short, ct.c_ushort, ct.c_int, \
              ct.c_uint, ct.c_long, ct.c_ulong, ct.c_longlong, \
              ct.c_ulonglong , ct.c_float, ct.c_double, ct.c_longdouble, \
              ct.c_char_p, ct.c_wchar_p, ct.c_void_p]

class shotfile:
  status = False
  def __init__( self , experiment=None , diagnostic=None , shotnumber=None , edition = -1):
    self.status = False
    if experiment != None and diagnostic != None and shotnumber != None:
      self.Open( experiment , diagnostic , shotnumber , edition )
  

  def __del__( self ):
    self.Close()
  
  
  def Open( self , experiment , diagnostic , shotnumber , edition = -1, mode="new"):
    self.Close()
    if shotnumber > 0:
      self.status = True
      vars(self)['diaref'] = ct.c_int32( 0 )
      vars(self)['edition'] = ct.c_int32( edition )
      vars(self)['shotnumber'] = ct.c_uint32( shotnumber )
      error = ct.c_int( 0 )
      _error = ct.byref( error )
      _experiment = ct.c_char_p( experiment )
      _diagnostic = ct.c_char_p( diagnostic )
      _shotnumber = ct.byref( self.shotnumber )
      _mode = ct.c_char_p( mode )
      _edition = ct.byref( self.edition )
      _diaref = ct.byref( self.diaref )
      date = datetime.date.strftime(datetime.datetime.today(),'%y.%m.%d;%H:%M:%S')
      _date = ct.c_char_p( date )
      lexpr = ct.c_ulonglong( len(experiment) )
      ldiag = ct.c_ulonglong( len(diagnostic) )
      lmode = ct.c_ulonglong( len( mode ) )
      ldate = ct.c_ulonglong( len( date ) )
      print diagnostic,len(diagnostic)
      print experiment,len(experiment)
      print date,len(date)
      print edition
      result = libddww.wwopen_(_error,_experiment,_diagnostic,_shotnumber,_mode,_edition,_diaref,_date,lexpr,ldiag,lmode,ldate)
      if self.GetError( error ):
        del self.diaref
        del self.edition
        del self.shotnumber
        self.status = False
        raise Exception('ww: Error Opening Shotfile %(error)s' % {'error':error.value})
    return self.status
  
  
  def Close( self ):
    if self.status:
      self.status = False
      error = ct.c_int32( 0 )
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      disp = 'lock'
      _disp = ct.c_char_p( disp )
      ldisp = ct.c_ulonglong( len(disp ) )
      space = 'maxspace'
      _space = ct.c_char_p( space )
      lspace = ct.c_ulonglong( len( space ) )
      result = libddww.wwclose_( _error , _diaref , _disp , _space , ldisp , lspace )
      print 'Close: ',result,error.value,' Edition',self.edition.value
      if self.GetError( error ):
        raise Exception('ww: Error Closing Shotfile')
      del self.diaref
      del self.edition
      del self.shotnumber
  
  
  def GetError( self , error ):
    if error.value != 0:
      print 'error!'
      _error = ct.byref( error )
      ID = ""
      _ID = ct.c_char_p( ID )
      lID = ct.c_ulonglong( len(ID) )
      ctrl = ct.c_uint32( 3 )
      _ctrl = ct.byref( ctrl )
      result = libddww.xxerror_( _error , _ctrl , _ID , lID )
      if result != 0:
        return True
    return False
  
  
  def GetArray( self , data ):
    if type(data) == type('str'):
      return ct.c_char_p(data)
    d = data
    while np.size(d) > 1: d = d[0]
    dtyp = d.dtype
    for cty in ctList:
      if np.dtype(cty) == dtyp:
        if np.size(data) > 1:
          return cty.from_buffer(data)
        return cty(data)
    return False
  
  
  def GetType( self , data ):
    # to do at some point:
    #  0 raw
    #  5 logical
    #  6 char*n
    # ---
    # done:
    d = data
    while np.size(d) > 1: d = d[0]
    if type(d) == type('str'):
      # 6 char
      return ct.c_uint32( 6 )
    if d.dtype in ("int8", "int16", "int32"):
      #  1 integer
      return ct.c_uint32( 1 )
    if d.dtype == "float32":
      #  2 float
      return ct.c_uint32( 2 )
    if d.dtype == "float64":
      #  3 double
      return ct.c_uint32( 3 )
    if d.dtype == "bool":
      #  2 float
      return ct.c_uint32( 5 )
    if d.dtype == "int64":
      # 10 long long
      return ct.c_uint32( 10 )
    raise Exception('ww: Unknown data type. known: int8, int16, int32, int64, float32, float64, char')
    return -1
  
  
  def SetSignal( self , signame , data ):
    if self.status:
      error = ct.c_int32( 0 )
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      name = ct.c_char_p( signame )
      lname = ct.c_ulonglong( len(signame) )
      temp = self.GetArray( data )
      type = self.GetType( data[0] )
      _type = ct.byref( type )
      lbuf = ct.c_uint32( np.size( data ) )
      _lbuf = ct.byref( lbuf )
      _data = ct.byref( temp )
      stride = ct.c_uint32( 1 )
      _stride = ct.byref( stride )
      result = libddww.wwsignal_(_error,_diaref,name,_type,_lbuf,_data,_stride,lname)
      if self.GetError( error ):
        raise Exception('ww: Error Writing Signal %(signame)s' %{'signame':signame} )
  
  
  def _wwinsert(self, signame, lbuf, data, stride, indices, typ=None):
    if self.status:
      error = ct.c_int32( 0 )
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      name = ct.c_char_p( signame )
      lname = ct.c_ulonglong( len(signame) )
      temp = self.GetArray( data )
      type = self.GetType( data[0] ) if typ==None else typ
      _type = ct.byref( type )
      lbuf = ct.c_uint32( np.size( data ) )
      _lbuf = ct.byref( lbuf )
      _data = ct.byref( temp )
      stride = ct.c_uint32( 1 )
      _stride = ct.byref( stride )
      nind=len(indices)
      c_indices = (ct.c_uint32 * nind)()
      for jind in range(nind):
        c_indices[jind]=ct.c_uint32(indices[jind])
      _indices = ct.byref(c_indices)
      result = libddww.wwinsert_(_error,_diaref,name,_type,_lbuf,_data,_stride,_indices,lname)
      if self.GetError( error ):
        raise Exception('ww: Error Inserting Signal %(signame)s' %{'signame':signame} )
    

  def _wwoinfo(self, signame):
    if self.status:
      error = ct.c_int32( 0 )
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      name = ct.c_char_p( signame )
      lname = ct.c_ulonglong( len(signame) )
      typ = ct.c_uint32(0)
      _typ = ct.byref( typ )
      format = ct.c_uint16(0)
      _format = ct.byref( format )
      ntval = ct.c_uint32(0)
      _ntval = ct.byref(ntval)
      items = ct.c_uint32(0)
      _items = ct.byref(items)
      indices = self.GetArray(np.array([0,0,0,0], dtype=np.int32))
      _indices = ct.byref(indices)      
      result = libddww.wwoinfo_(_error,_diaref,name,_typ,_format,_ntval,_items,_indices,lname)
      return {'name':name.value, 'type':typ, 'format':format, 'leng':ntval, 'items':items, 'indices':indices}
      if self.GetError( error ):
        raise Exception('ww: Error Inserting Signal %(signame)s' %{'signame':signame} )
      
  
  def SetSignalGroup( self, signame, data):
    dat = np.array(data)
    # dat[time, dimension1, dimension2, dimension3]
    depth = len(dat.shape)
    print signame, depth, dat.shape
    if depth < 2:
      raise Exception('ww: SetSignalGroup requires a Signal_Group_')
      return False
    dim0 = dat.shape[0]
    dim1 = dat.shape[1]
    dim2 = dat.shape[2] if depth > 2 else 1
    dim3 = dat.shape[3] if depth > 3 else 1
    if depth == 2:
      for i in range(dim1):
        single_dat = dat[:, i]
        indices = np.array([i+1,1,1],dtype=np.int64)
        self._wwinsert(signame, dim0, np.array(single_dat), 1, indices)
    if depth == 3:
      for i in range(dim1):
        for j in range(dim2):
          single_dat = dat[:,i,j]
          indices = np.array([1+i,1+j,1],dtype=np.int64)
          self._wwinsert(signame, dim0, np.array(single_dat), 1, indices)
    
    return True
  
  def SetAreabase(self, areaname, nt, data):
    if self.status:
      error = ct.c_int32( 0 )
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      name = ct.c_char_p( areaname )
      lname = ct.c_ulonglong( len( areaname ) )
      type = self.GetType( data[0] )
      _type = ct.byref( type )
      k1  = ct.c_long(1)
      k2  = ct.c_long(nt)
      _k1 = ct.byref(k1)
      _k2 = ct.byref(k2)
      if type == -1:
        raise Exception("ww.shotfile.SetAreabase: Unsupported datatype %(type)s in areabase %(name)s." % \
          {'type':np.dtype( data[0] ),'name':areaname} )
      if nt == 1:
        sizes = (ct.c_uint32 * 3)()
        depth = len(data.shape)
        for jsiz in range(depth):
          sizes[jsiz] = data.shape[jsiz]
        print 'Size of AB '+areaname+' :',sizes[0:3]
        _sizes = ct.byref( sizes )
        temp = self.GetArray( data )
        _data = ct.byref( temp )
        libddww.wwainsert_(_error,_diaref,name,_k1,_k2,_type,_data,_sizes,lname)
        print 'Written AB '+areaname
      else:
        raise Exception('ww.shotfile.SetAreabase: not prepared yet for time dependent areabase %(name)s' % {'name':areaname})

      if self.GetError( error ):
        raise Exception('ww.shotfile.SetAreabase: Error writing areabase %(name)s' % {'name':areaname})

  def SetTimebase( self , timename , data ):
    if self.status:
      error = ct.c_int32( 0 )
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      name = ct.c_char_p( timename )
      lname = ct.c_ulonglong( len( timename ) )
      type = self.GetType( data[0] )
      _type = ct.byref( type )
      if type == -1:
        raise Exception("ww.shotfile.SetTimebase: Unsupported datatype %(type)s in timebase %(name)s." % \
          {'type':np.dtype( data[0] ),'name':timename} )
      lbuf = ct.c_uint32( np.size( data ) )
      _lbuf = ct.byref( lbuf )
      temp = self.GetArray( data )
      _data = ct.byref( temp )
      stride = ct.c_uint32( 1 )
      _stride = ct.byref( stride )
      libddww.wwtbase_(_error,_diaref,name,_type,_lbuf,_data,_stride,lname)
      if self.GetError( error ):
        raise Exception('ww.shotfile.SetTimebase: Error writing timebase %(name)s' % {'name':timename})
  
  
  def GetInfo( self , objname ):
    if self.status:
      error = ct.c_int32( 0 )
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      name = ct.c_char_p( objname )
      lname = ct.c_ulonglong( len(objname) )
      typ = ct.c_uint32( 0 )
      _typ = ct.byref( typ )
      format = ct.c_uint16( 0 )
      _format = ct.byref( format )
      leng = ct.c_uint32( 0 )
      _leng = ct.byref( leng )
      items = ct.c_uint32( 0 )
      _items = ct.byref( items )
      indices = (ct.c_uint32*3)()
      _indices = ct.byref( indices )
      libddww.wwoinfo_(_error,_diaref,name,_typ,_format,_leng,_items,_indices,lname)
      if self.GetError( error ):
        raise Exception('ww.shotfile.GetInfo: Error getting info of object %(name)s' %{'name':objname})
      index = np.frombuffer( indices , dtype = np.uint32 )
      return ObjectInfo(objname,typ.value,format.value,leng.value,items.value,index)


  def SetParameter( self , set_name , par_name , payload, typ=None, stride=1, cjust=8 ):
    """ ww.shotfile.SetParameter( SetName , ParameterName , payload )\n
    types: 0 = raw 1 = integer 2 = float 3 = double 4 = complex 5 = logical 6 = character
    """
    if self.status:
      error = ct.c_int32(0)
      _error = ct.byref( error )
      _diaref = ct.byref( self.diaref )
      pset = ct.c_char_p(set_name)
      lset = ct.c_ulonglong( len( set_name ) )
      par = ct.c_char_p(par_name)
      lpar = ct.c_ulonglong( len( par_name ) )
      mytype = self.GetType( payload ) if typ==None else ct.c_uint32(typ)
      _type = ct.byref( mytype )
      if mytype.value == 6:
        if type(payload) == type('str'):
          sbuf=payload.ljust(cjust)
        else:
          sbuf = ''
          for entry in payload:
            sbuf += entry.ljust(cjust)
        lbuf = ct.c_int32( len(sbuf) )
        _lbuf = ct.byref( lbuf )
        buffer = ct.c_char_p(sbuf)
        stride = ct.c_uint16(1)
        _stride = ct.byref(stride)
        result = libddww.wwparm_(_error,_diaref, pset, par, _type, _lbuf, buffer, _stride, lset, lpar, lbuf)
      else:
        lbuf = ct.c_int32( np.size(payload))
        _lbuf = ct.byref( lbuf )
        buffer = self.GetArray(payload)
        _buffer = ct.byref( buffer )
        stride = ct.c_uint16(1)
        _stride = ct.byref(stride)
        result = libddww.wwparm_(_error,_diaref, pset, par, _type, _lbuf, _buffer, _stride, lset, lpar)
      if self.GetError( error ):
        raise Exception('ww: Error Writing Parameters into %(setname)s -> %(parname)s' %{'setname':set_name, 'parname':par_name} )


  def write_sf(self, nshot, sfh_d, sfhdir, diag='TRA'):

    print 'Write_sf'
    exp = os.getenv('USER')
    os.chdir(sfhdir)
    if self.Open(exp,diag,nshot):

      for obj,tmp in sfh_d.iteritems():
        if tmp['devtyp'] == 8:
          print 'Writing Timebase '+obj
          self.SetTimebase( obj,np.float32(sfh_d[obj]['data']) )

      for obj,tmp in sfh_d.iteritems():
        if tmp['devtyp'] == 13:
          print 'Writing Areabase '+obj
          adim=sfh_d[obj]['data'].shape
          if len(adim) == 1:
            nta=1
            self.SetAreabase( obj,nta,np.float32(sfh_d[obj]['data']) )
          if len(adim) == 2:
            nta=adim[0]
            self.SetAreabase( obj,1,np.float32(sfh_d[obj]['data'][0,:]) )
#            self.SetAreabase( obj,nta,np.float32(sfh_d[obj]['data'][0,:]) )

      for obj,tmp in sfh_d.iteritems():
# Parameter set
        if tmp['devtyp'] == 4:
          for pn,val in tmp['info'].iteritems():
            fmt   = val[0]
            n_sfh = val[1]
            n_nml = len(tmp[pn]['data'])
            print pn, n_sfh, n_nml
            if fmt == 1794:
              if n_nml == 0:
                self.SetParameter(obj, pn, 'None')
              if n_nml >= 1:
                self.SetParameter(obj, pn, tmp[pn]['data'], typ=6)
            else:
              if n_sfh <= n_nml:
                if fmt in [3,4]:
                  par = np.zeros(n_sfh, dtype=np.int32)
                if fmt == 5:
                  par=np.zeros(n_sfh, dtype=np.float32)
                if fmt == 7:
                  par=np.zeros(n_sfh,dtype=np.bool)
                par[0:n_nml]=tmp[pn]['data']
                self.SetParameter(obj, pn, par)
# Signal
        if tmp['devtyp'] == 7:
          self.SetSignal( obj,np.float32(tmp['data']) )
# SignalGroup
        if tmp['devtyp'] == 6:
          self.SetSignalGroup( obj,np.float32(tmp['data']) )
      self.Close()        
