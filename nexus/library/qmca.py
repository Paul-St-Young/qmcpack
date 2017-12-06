#! /usr/bin/env python

import code
import os
from optparse import OptionParser
import sys
from sys import exit
import traceback
from copy import deepcopy


class Callable:                        #helper class to create 'static' methods
    def __init__(self, anycallable):
        self.__call__ = anycallable
    #end def init
#end class Callable

import copy
class Copyable:
    def _copy(self):
        return copy.deepcopy(self)
        #return copy.copy(self)
    #end def copy

    def _copies(self,count):
        copies = []
        for i in range(count):
            copies.append(self._copy())
        #end for
        return copies
    #end def copies
#end class Copyable


class Comparable:
    None
#end class Comparable


class Iterable:
    def __iter__(self):
        for item in self.__dict__:
            yield self.__dict__[item]
        #end for
    #end def __iter__

    def _iteritems(self):
        return self.__dict__.iteritems()
    #end def iteritems

    def _keys(self):
        return self.__dict__.keys()
    #end def _keys

    def _values(self):
        return self.__dict__.values()
    #end def _keys

    def _inverse(self):
        new = self.__class__()
        new.__dict__ = dict((v,k) for k, v in self._iteritems())
        return new
    #end def _inverse
#end class Iterable


import inspect
class Queryable(Iterable):
    def __str__(self,nindent=1):
        pad = '  '
        npad = nindent*pad
        s=''
        stype = type(s)
        normal = []
        qable  = []
        for k,v in self._iteritems():
            if type(k)!=stype or k[0]!='_':
                if isinstance(v,Queryable):
                    qable.append(k)
                else:
                    normal.append(k)
                #end if
            #end if
        #end for
        normal.sort()
        qable.sort()
        for k in normal:
            v = self[k]
            indent = npad+18*' '
            vstr = str(v).replace('\n','\n'+indent)
            s+=npad+'{0:<15} = '.format(k)+vstr+'\n'
        #end for
        for k in qable:
            v = self[k]
            s+=npad+str(k)+'\n'
            s+=v.__str__(nindent+1)
            if isinstance(k,str):
                s+=npad+'end '+k+'\n'
            #end if
        #end for
        return s
    #end def __str__

    def __repr__(self):
        s=''
        stype = type(s)
        mem = list(self.__dict__.keys())
        mem.sort()
        for m in mem:
            if type(m)!=stype or m[0]!='_':
#                s+='  '+m+'   '+str(type(self.__dict__[m]))+'\n'
#                s+='  '+m+'   '+self.__dict__[m].__class__.__name__+'\n'
                v=self.__dict__[m]
                if hasattr(v,'__class__'):
                    s+='  {0:<20}  {1:<20}\n'.format(m,v.__class__.__name__)
                else:
                    s+='  {0:<20}  {1:<20}\n'.format(m,type(v))
                #end if
            #end if
        #end for
        return s
    #end def __repr__

    def __len__(self):
        return len(self.__dict__)
    #end def __len__

    def __contains__(self,name):
        return name in self.__dict__
    #end def

    def __getitem__(self,name):
        return self.__dict__[name]
    #end def __getitem__

    def __setitem__(self,name,value):
        self.__dict__[name]=value
    #end def __setitem__

    def __delitem__(self,name):
        del self.__dict__[name]
    #end def __setitem__

    def _set(self,**variables):
        for name,value in variables.iteritems():
            self[name]=value
        #end for
        return self
    #end def _set

    def _clear(self):
        self.__dict__.clear()
    #end def _clear

    def _add_attribute(self,name,value=None):
        self[name] = value
    #end def _add_attribute

    def _add_attributes(self,*names,**attributes):
        for name in names:
            self[name] = None
        #end for
        for name,value in attributes.iteritems():
            self[name]=value
        #end for
    #end def _add_attributes

    def _transfer_from(self,other,copy=False):
        if not copy:
            for k,v in other._iteritems():
                self[k]=v
            #end for
        else:
            for k,v in other._iteritems():
                self[k]=deepcopy(v)
            #end for
        #end if
    #end def _transfer_from

    def _transfer_to(self,other,copy=False):
        other._transfer_from(self,copy)
    #end def _transfer_to

    def _append(self,value):
        self[len(self)] = value
    #end def _append
#end class Queryable


import sys
class Logable:
    def _open_log(self,filepath,increment=0,prefix=''):
        self._log_file=open(filepath,'w')
        self._log_pad_prefix=prefix
        self._log_pad_increment=increment
        self._log_pad_shift=0
        self._log_pad_update()
    #end def _open_log

    def _close_log(self):
        if self._log_file!=sys.stdout:
            self._log_file.close()
        #end if
    #end def _close_log

    def _set_log(self,logfile=sys.stdout,increment=0,prefix=''):
        self._log_file=logfile
        self._log_pad_prefix=prefix
        self._log_pad_increment=increment
        self._log_pad_shift=0
        self._log_pad_update()
    #end def _set_log

    def _get_log(self):
        L = Logable()
        L._log_file          = self._log_file               
        L._log_pad_prefix    = str(self._log_pad_prefix)
        L._log_pad_increment = int(self._log_pad_increment)
        L._log_pad_shift     = int(self._log_pad_shift)
        return L
    #end def _get_log

    def _transfer_log(self,L):
        self._log_file          = L._log_file               
        self._log_pad_prefix    = L._log_pad_prefix   
        self._log_pad_increment = L._log_pad_increment
        self._log_pad_shift     = L._log_pad_shift    
        self._log_pad_update()
    #end def _transfer_log

    def _increment_pad(self):
        self._log_pad_shift+=self._log_pad_increment
        self._log_pad_update()
    #end def _increment_pad

    def _decrement_pad(self):
        self._log_pad_shift-=self._log_pad_increment
        self._log_pad_update()
    #end def _decrement_pad

    def _set_pad(self,increment=0,prefix=''):
        self._log_pad_prefix=prefix
        self._log_pad_increment=increment
        self._log_pad_update()
    #end def _set_pad

    def _log(self,s=''):
        self._log_file.write(self._log_pad+str(s)+'\n')
    #end def _log

    def _exit(self,s=''):
        self._log(self.__class__.__name__+' exiting')
        self._log(s)
        sys.exit()
    #end def _exit

    def _error(self,s='',exit_now=False):
        self._log(self.__class__.__name__+' Error:  '+s)
        if exit_now:
            self._exit()
        #end if
    #end def _error

    def _warn(self,s):
        self._log(self.__class__.__name__+' Warning:  '+s)
        return
    #end def _warn


    #user should not need to use these functions
    def _log_pad_update(self):
        self._log_pad = self._log_pad_prefix+self._log_pad_shift*' '
    #end def _log_pad_update
#end class Logable


import cPickle
class Savable:
    def _save(self,fpath=None):
        #self._unlink_dynamic_methods(self)
        if fpath==None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        #self._savepath=fpath
        fobj = open(fpath,'w')
        binary = cPickle.HIGHEST_PROTOCOL
        cPickle.dump(self,fobj,binary)
        # fobj.close()
        del fobj
        del binary
        return
    #end def _save

    def _load(self,fpath=None):
        if fpath==None:
            fpath='./'+self.__class__.__name__+'.p'
        #end if
        fobj = open(fpath,'r')
        tmp = cPickle.load(fobj)
        # fobj.close()
        self.__dict__.clear()
        #self.__dict__=tmp.__dict__
        for k,v in tmp.__dict__.iteritems():
            self.__dict__[k] = v
        #end for
        del fobj
        del tmp
        #self._relink_dynamic_methods()
        return
    #end def _load

    #def _unlink_dynamic_methods(self,obj):
    #    if hasattr(obj,'__dict__'):
    #        for k,v in obj.__dict__.iteritems():
    #            if type(v)==type(self._unlink_dynamic_methods):
    #                obj.__dict__[k]=None
    #            else:
    #                self._unlink_dynamic_methods(v)
    #            #end if
    #        #end for
    #    elif hasattr(obj,'__iter__'):
    #        for a in obj:
    #            self._unlink_dynamic_methods(a)
    #        #end for
    #    #end if
    #    return
    ##end def _unlink_dynamic_methods
    #
    #def _relink_dynamic_methods(self,obj):
    #    if hasattr(obj,'_set_dynamic_methods'):
    #        obj._set_dynamic_methods()
    #    #end if
    #    if hasattr(obj,'__dict__'):
    #        for a in obj.__dict__:
    #            self._relink_dynamic_methods(a)
    #        #end for
    #    elif hasattr(obj,'__iter__'):
    #        for a in obj:
    #            self._relink_dynamic_methods(a)
    #        #end for
    #    #end if
    #    return
    ##end def _relink_dynamic_methods
    #
    #def _set_dynamic_methods(self):
    #    return
    ##end def _set_dynamic_methods
#end class Savable


class Validatable:
    def _has_type_bindings(self):
        return '_type_bindings' in self.__class__.__dict__
    #end def _has_type_bindings

    def _has_rules(self):
        return '_rules' in self.__class__.__dict__
    #end def _has_rules

    def _bind_type(self,vname,vtype):
        if not self._has_type_bindings():
            self.__class__.__dict__['_type_bindings']=dict()
        #end if
        self.__class__.__dict__['_type_bindings'][vname]=vtype
        return
    #end def _bind_type

    #rules must have an executable form, eg. vrule='self.var > 1.0'
    def _add_rule(self,vrule):
        if not self._has_rules():
            self.__class__.__dict__['_rules']=list()
        #end if
        self.__class__.__dict__['_rules'].append(vrule)
        return
    #end def

    def _is_complete(self):
        complete = True
        for k,v in self.__dict__.iteritems():
            if k[0]!='_':
                complete = complete and v!=None
            #end if
        #end for
        return complete
    #end def _is_complete

    def _is_type_valid(self,return_exceptions=False):
        type_valid=True
        exceptions = dict()
        if self._has_type_bindings():
            for k,v in self.__class__.__dict__['_type_bindings'].iteritems():
                var = self.__dict__[k]
                if hasattr(v,'__class__'):
                    vtype = var.__class__.__name__
                else:
                    vtype = type(var)
                #end if
                var_valid = vtype==v
                type_valid = type_valid and var_valid
                if not var_valid:
                    exceptions[k]=v,vtype
                #end if
            #end for
        #end if
        if not return_exceptions:
            return type_valid
        else:
            return type_valid,exceptions
        #end if
    #end def _is_type_valid

    def _is_legal(self):
        legal=True
        if self._has_rules():
            for rule in self.__class__.__dict__['_rules']:
                exec 'obeys_rule='+rule
                legal = legal and obeys_rule
            #end for
        #end if
        return legal
    #end def _is_legal

    def _is_valid(self,exit_on_fail=False,return_exceptions=False):
        complete   = self._is_complete()
        type_valid,type_exceptions = self._is_type_valid(return_exceptions=True)
        legal      = self._is_legal()
        valid = complete and type_valid and legal
        if not valid and exit_on_fail:
            name = self.__class__.__name__
            if not complete:
                print '       '+name+' has None attributes'
            #end if
            if not type_valid:
                print '       '+name+' has invalid attribute types'
                self._write_type_exceptions(type_exceptions)
            #end if
            if not legal:
                print '       '+name+' violates user-defined rules'
            #end if
            sys.exit('Error: '+name+' failed _is_valid, exiting')
        #end if
        return valid
    #end def _is_valid

    def _write_type_exceptions(self,type_exceptions,pad='         '):
        for k,v in type_exceptions.iteritems():
            print pad+'variable '+k+' should be '+v[0]+' type, is '+v[1]+' type'
        #end for
    #end def _write_type_exceptions        
#end class Validatable


class AllAbilities(Copyable,Comparable,Queryable,Savable,Validatable):
    def __init__(self,*vals,**attributes):
        names = []
        for val in vals:
            if isinstance(val,str):
                names.append(val)
            elif isinstance(val,Iterable):
                for name,value in val._iteritems():
                    self[name]=value
                #end for
            elif isinstance(val,dict):
                for name,value in val.iteritems():
                    self[name]=value
                #end for
            #end for
        #end for
        self._add_attributes(*names,**attributes)
    #end def __init__
#end class AllAbilities



exit_call = exit

class obj(AllAbilities):
    logfile = sys.stdout

    def copy(self):
        return self._copy()
    def copies(self,count):
        return self._copies(count)
    def iteritems(self):
        return self._iteritems()
    def keys(self):
        return self._keys()
    def values(self):
        return self._values()
    def inverse(self):
        return self._inverse()
    def set(self,**variables):
        return self._set(**variables)
    def clear(self):
        self._clear()
    def add_attribute(self,name,value=None):
        self._add_attribute(name,value)
    def add_attributes(self,*names,**attributes):
        self._add_attributes(*names,**attributes)
    #def transfer_from(self,other,copy=False):
    #    self._transfer_from(other,copy)
    #def transfer_to(self,other,copy=False):
    #    self._transfer_to(other,copy)
    def append(self,value):
        self._append(value)
    def save(self,fpath=None):
        self._save(fpath)
    def load(self,fpath):
        self._load(fpath)


    def list(self,*names):
        if len(names)==0:
            names = list(self.keys())
            names.sort()
        #end if
        values = []
        for name in names:
            values.append(self[name])
        #end if
        return values
    #end def list

    def tuple(self,*names):
        return tuple(self.list(*names))
    #end def tuple

    def obj(self,*names):
        o = obj()
        o.transfer_from(self,keys=names,copy=False)
        return o
    #end def obj

    def first(self):
        return self[min(self.keys())]
    #end def first

    def last(self):
        return self[max(self.keys())]
    #end def last
    
    def open_log(self,filepath):
        self.logfile = open(filepath,'w')
    #end def open_log

    def close_log(self):
        self.logfile.close()
    #end def close_log

    def _write(self,s):
        self.logfile.write(s)
    #end def _write

    def write(self,s):
        self._write(s)
    #end def write

    def logna(self,*items):
        s=''
        for item in items:
            s+=str(item)+' '
        #end for
        self.logfile.write(s)
    #end def logna

    def log(self,*items):
        s=''
        for item in items:
            s+=str(item)+' '
        #end for
        s+='\n'
        self.logfile.write(s)
    #end def log

    def error(self,message,header=None,exit=True,trace=False,post_header=' Error:'):
        pad = 4*' '
        if header==None:
            header = self.__class__.__name__
        #end if
        self.log(header+post_header)
        self.log(pad+message.replace('\n','\n'+pad))
        if exit:
            self.log('  exiting.\n')
            if trace:
                traceback.print_stack()
            #end if
            exit_call()
        #end if
    #end def error

    def warn(self,message,header=None,post_header=' Warning:'):
        pad = 4*' '
        if header==None:
            header=self.__class__.__name__
        #end if
        self.log(header+post_header)
        self.log(pad+message.replace('\n','\n'+pad))
    #end def error

    @classmethod
    def class_error(cls,message,header=None,exit=True,trace=True,post_header=' Error:'):
        pad = 4*' '
        if header==None:
            header = cls.__name__
        #end if
        cls.logfile.write(header+post_header)
        cls.logfile.write(pad+message.replace('\n','\n'+pad)+'\n')
        if exit:
            cls.logfile.write('  exiting.\n\n')
            if trace:
                traceback.print_stack()
            #end if
            exit_call()
        #end if
    #end def class_error

    def transfer_from(self,other,keys=None,copy=False):
        if keys==None:
            keys = other.keys()
        #end if
        if not copy:
            for k in keys:
                self[k]=other[k]
            #end for
        else:
            for k in keys:
                self[k]=deepcopy(other[k])
            #end for
        #end if
    #end def transfer_from

    def transfer_to(self,other,keys=None,copy=False):
        if keys==None:
            keys = self.keys()
        #end if
        if not copy:
            for k in keys:
                other[k]=self[k]
            #end for
        else:
            for k in keys:
                other[k]=deepcopy(self[k])
            #end for
        #end if
    #end def transfer_to

    def move_from(self,other,keys=None):
        if keys==None:
            keys = other.keys()
        #end if
        for k in keys:
            self[k]=other[k]
            del other[k]
        #end for
    #end def move_from

    def move_to(self,other,keys=None):
        other.move_from(self,keys)
    #end def move_to

    def copy_from(self,other,keys=None,deep=False):
        self.transfer_from(other,keys,copy=deep)
    #end def copy_from

    def copy_to(self,other,keys=None,deep=False):
        self.transfer_to(other,keys,copy=deep)
    #end def copy_to

#end class obj




class DevBase(obj):
    user_interface_data = obj()
    dev_instructions_data = obj()

    def not_implemented(self):
        #self.error('a base class function has not been implemented','Developer')
        self.error('a base class function has not been implemented',trace=True)
    #end def not_implemented

    @classmethod
    def set_user_interface(cls,class_variables=None,class_functions=None,member_variables=None,member_functions=None):
        if class_variables is None:
            class_variables = []
        if class_functions is None:
            class_functions = []
        if member_variables is None:
            member_variables = []
        if member_functions is None:
            member_functions = []
        ui = cls.user_interface_data
        ui.class_variables = class_variables
        ui.class_functions = class_functions
        ui.member_variables= member_variables
        ui.member_functions= member_functions
    #end def set_user_interface

    @classmethod
    def set_dev_instruction(cls,situation='writing a derived class',class_variables=None,class_functions=None,member_variables=None,member_functions=None):
        if class_variables is None:
            class_variables = []
        if class_functions is None:
            class_functions = []
        if member_variables is None:
            member_variables = []
        if member_functions is None:
            member_functions = []
        ins = obj()
        cls.dev_instructions_data[situation] = ins
        ins.class_variables = class_variables
        ins.class_functions = class_functions
        ins.member_variables= member_variables
        ins.member_functions= member_functions
    #end def set_dev_instruction
        
    @classmethod
    def write_dev_data(cls,s,dev_data,indent=''):
        p1 = '  '
        p2 = 2*p1
        p3 = 3*p1
        dpad = indent+p3
        for name,value in dev_data.iteritems():
            s+=indent+p1+name+'\n'
            if isinstance(value,list):
                for v in value:
                    s+=indent+p2+v+'\n'
                #end for
            else:
                for v,description in value.iteritems():
                    s+=indent+p2+v+'\n'
                    s+=dpad+description.replace('\n','\n'+dpad)
                #end for
            #end for
        #end for
        return s
    #end def write_dev_data

    @classmethod
    def user_interface(cls):
        s='User Interface for '+cls.__name__+' Class\n'
        s+= (len(s)-1)*'='+'\n'
        s+=cls.write_dev_data(cls.user_interface_data)
        return s
    #end def user_interface

    @classmethod
    def developer_instructions(cls):
        s='Developer Instructions for '+cls.__name__+' Class\n'
        s+= (len(s)-1)*'='+'\n'
        i1='  '
        i2=2*i1
        i3=3*i1
        for situation,instruction in cls.developer_instructions.iteritems():
            s+=i1+'situation: '+situation+'\n'
            s+=i2+'define the following variables and functions:'
            s+=cls.write_dev_data(instructions,i3)
        #end for
        return s
    #end def developer_instructions
#end class DevBase








class Void:
    void_items = dict()

    @classmethod
    def _unavailable(cls,self):
        sid = id(self)
        if sid in Void.void_items:
            module,item = Void.void_items[id(self)]
        else:
            module,item = None,None
        #end if
        if module is None and item is None:
            msg = 'encountered a void item from an unavailable module'
        elif module is None:
            msg = 'item '+str(item)+' is from an unavailable module'
        elif item is None:
            msg = 'encountered a void item from unavailable module '+str(module)+'  \nthis python module must be installed on your system to use this feature'
        else:
            msg = 'item '+str(item)+' is from unavailable module '+str(module)+'  \nthis python module must be installed on your system to use this feature'
        #end if
        DevBase.class_error(msg,'QMCA',trace=False)
    #end def _unavailable

        
    @classmethod
    def _class_unavailable(cls):
        msg = 'encountered a void item from an unavailable module'
        DevBase.class_error(msg,'QMCA',trace=False)
    #end def _class_unavailable
        


    def __init__(self,module=None,item=None):
        Void.void_items[id(self)] = module,item
    #end def __init__


    #list of magic functions taken from the following sources
    #  http://web.archive.org/web/20110131211638/http://diveintopython3.org/special-method-names.html
    #  http://www.ironpythoninaction.com/magic-methods.html


    #class methods
    @classmethod
    def __instancecheck__(cls,*args,**kwargs):
        Void._class_unavailable()
    def __subclasscheck__(cls,*args,**kwargs):
        Void._class_unavailable()
    def __subclasshook__(cls,*args,**kwargs):
        Void._class_unavailable()
    

    #member methods
    def __new__(self,*args,**kwargs):
        Void._unavailable(self)
    def __eq__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ne__(self,*args,**kwargs):
        Void._unavailable(self)
    def __lt__(self,*args,**kwargs):
        Void._unavailable(self)
    def __le__(self,*args,**kwargs):
        Void._unavailable(self)
    def __gt__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ge__(self,*args,**kwargs):
        Void._unavailable(self)
    def __nonzero__(self,*args,**kwargs):
        Void._unavailable(self)
    def __subclasses__(self,*args,**kwargs):
        Void._unavailable(self)
    def __call__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hash__(self,*args,**kwargs):
        Void._unavailable(self)
    #def __del__(self,*args,**kwargs):
    #    Void._unavailable(self)
    def __dir__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delitem__(self,*args,**kwargs):
        Void._unavailable(self)
    def __len__(self,*args,**kwargs):
        Void._unavailable(self)
    def __contains__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iter__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reversed__(self,*args,**kwargs):
        Void._unavailable(self)
    def __missing__(self,*args,**kwargs):
        Void._unavailable(self)
    def __length_hint__(self,*args,**kwargs):
        Void._unavailable(self)
    def __repr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __str__(self,*args,**kwargs):
        Void._unavailable(self)
    def __unicode__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delattr__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getattribute__(self,*args,**kwargs):
        Void._unavailable(self)
    def __add__(self,*args,**kwargs):
        Void._unavailable(self)
    def __sub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __mul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __floordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __div__(self,*args,**kwargs):
        Void._unavailable(self)
    def __truediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __mod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __divmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __pow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __lshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __and__(self,*args,**kwargs):
        Void._unavailable(self)
    def __xor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __or__(self,*args,**kwargs):
        Void._unavailable(self)
    def __neg__(self,*args,**kwargs):
        Void._unavailable(self)
    def __pos__(self,*args,**kwargs):
        Void._unavailable(self)
    def __abs__(self,*args,**kwargs):
        Void._unavailable(self)
    def __invert__(self,*args,**kwargs):
        Void._unavailable(self)
    def __complex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __int__(self,*args,**kwargs):
        Void._unavailable(self)
    def __long__(self,*args,**kwargs):
        Void._unavailable(self)
    def __float__(self,*args,**kwargs):
        Void._unavailable(self)
    def __oct__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __index__(self,*args,**kwargs):
        Void._unavailable(self)
    def __enter__(self,*args,**kwargs):
        Void._unavailable(self)
    def __exit__(self,*args,**kwargs):
        Void._unavailable(self)
    def __get__(self,*args,**kwargs):
        Void._unavailable(self)
    def __set__(self,*args,**kwargs):
        Void._unavailable(self)
    def __delete__(self,*args,**kwargs):
        Void._unavailable(self)
    def __doc__(self,*args,**kwargs):
        Void._unavailable(self)
    def __dict__(self,*args,**kwargs):
        Void._unavailable(self)
    def __slots__(self,*args,**kwargs):
        Void._unavailable(self)
    def __class__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bases__(self,*args,**kwargs):
        Void._unavailable(self)
    def __name__(self,*args,**kwargs):
        Void._unavailable(self)
    def __all__(self,*args,**kwargs):
        Void._unavailable(self)
    def __file__(self,*args,**kwargs):
        Void._unavailable(self)
    def __module__(self,*args,**kwargs):
        Void._unavailable(self)
    #def __metaclass__(self,*args,**kwargs):
    #    Void._unavailable(self)
    def __import__(self,*args,**kwargs):
        Void._unavailable(self)
    def __radd__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rsub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rmul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rtruediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rfloordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rdivmod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rpow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rlshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rrshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rand__(self,*args,**kwargs):
        Void._unavailable(self)
    def __rxor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ror__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iadd__(self,*args,**kwargs):
        Void._unavailable(self)
    def __isub__(self,*args,**kwargs):
        Void._unavailable(self)
    def __imul__(self,*args,**kwargs):
        Void._unavailable(self)
    def __itruediv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ifloordiv__(self,*args,**kwargs):
        Void._unavailable(self)
    def __imod__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ipow__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ilshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __irshift__(self,*args,**kwargs):
        Void._unavailable(self)
    def __iand__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ixor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ior__(self,*args,**kwargs):
        Void._unavailable(self)
    def __round__(self,*args,**kwargs):
        Void._unavailable(self)
    def __ceil__(self,*args,**kwargs):
        Void._unavailable(self)
    def __floor__(self,*args,**kwargs):
        Void._unavailable(self)
    def __trunc__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bool__(self,*args,**kwargs):
        Void._unavailable(self)
    def __copy__(self,*args,**kwargs):
        Void._unavailable(self)
    def __deepcopy__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getstate__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reduce__(self,*args,**kwargs):
        Void._unavailable(self)
    def __reduce_ex__(self,*args,**kwargs):
        Void._unavailable(self)
    def __getnewargs__(self,*args,**kwargs):
        Void._unavailable(self)
    def __setstate__(self,*args,**kwargs):
        Void._unavailable(self)
    def __hash__(self,*args,**kwargs):
        Void._unavailable(self)
    def __bytes__(self,*args,**kwargs):
        Void._unavailable(self)
    def __format__(self,*args,**kwargs):
        Void._unavailable(self)
    def __next__(self,*args,**kwargs):
        Void._unavailable(self)
#end class Void


def unavailable(module,*items):
    voids = []
    if len(items)==0:
        voids.append(Void(module))
    #end if
    for item in items:
        voids.append(Void(module,item))
    #end for
    if len(voids)==1:
        return voids[0]
    else:
        return voids
    #end if
#end def unavailable

try:
    from numpy import array,zeros,dot,loadtxt,floor,empty,sqrt,trace,savetxt,concatenate,real,imag,diag,arange,ones,identity,exp,log,log10,sign
except (ImportError,RuntimeError):
    array,zeros,dot,loadtxt,floor,empty,sqrt,trace,savetxt,concatenate,real,imag,diag,arange,ones,identity,exp,log,log10,sign = unavailable('numpy','array','zeros','dot','loadtxt','floor','empty','sqrt','trace','savetxt','concatenate','real','imag','diag','arange','ones','identity','exp','log','log10','sign')
#end try
try:
    import numpy.random as random
except:
    random = unavailable('numpy.random','random')
#end try


class Unit:
    def __init__(self,name,symbol,type,value,shift=0):
        self.name = name
        self.symbol = symbol
        self.type = type
        self.value = value
        self.shift = shift
    #end def __init__
#end class Unit


class UnitConverter:

    unassigned = None


    kb = 1.3806503e-23 #J/K

    count_set = set(['mol'])
    mol = 6.0221415e23

    distance_set = set(['m','A','B','nm','pm','fm','a','b','c'])
    m  = 1.e0
    A  = 1.e-10*m
    B  = .52917720859e-10*m
    nm = 1.e-9*m
    pm = 1.e-12*m
    fm = 1.e-15*m

    time_set = set(['s','ms','ns','ps','fs'])
    s = 1.e0
    ms = 1.e-3*s
    ns = 1.e-9*s
    ps = 1.e-12*s
    fs = 1.e-15*s

    mass_set = set(['kg','me','mp','amu','Da'])
    kg  = 1.e0
    me  = 9.10938291e-31*kg
    mp  = 1.672621777e-27*kg
    amu = 1.660538921e-27*kg
    Da  = amu

    energy_set = set(['J','eV','Ry','Ha','kJ_mol','K','degC','degF'])
    J      = 1.e0
    eV     = 1.60217646e-19*J
    Ry     = 13.6056923*eV
    Ha     = 2*Ry
    kJ_mol = 1000.*J/mol
    K      = J/kb
    degC   = K
    degF   = 9./5.*K

    degC_shift = -273.15
    degF_shift = -459.67 

    charge_set = set(['C','e'])
    C = 1.e0
    e = 1.60217646e-19*C

    pressure_set = set(['Pa','bar','Mbar','GPa','atm'])
    Pa   = 1.e0
    bar  = 1.e5*Pa
    Mbar = 1.e6*bar
    GPa  = 1.e9*Pa
    atm  = 1.01325e5*Pa

    force_set = set(['N','pN'])
    N  = 1.e0
    pN = 1e-12*N

    therm_cond_set = set(['W_mK'])
    W_mK = 1.0

    alatt = unassigned
    blatt = unassigned
    clatt = unassigned

    meter      = Unit('meter'     ,'m' ,'distance',m)
    angstrom   = Unit('Angstrom'  ,'A' ,'distance',A)
    bohr       = Unit('Bohr'      ,'B' ,'distance',B)
    nanometer  = Unit('nanometer' ,'nm','distance',nm)
    picometer  = Unit('picometer' ,'pm','distance',pm)
    femtometer = Unit('femtometer','pm','distance',fm)
    a          = Unit('a'         ,'a' ,'distance',alatt)
    b          = Unit('b'         ,'b' ,'distance',blatt)
    c          = Unit('c'         ,'c' ,'distance',clatt)

    second      = Unit('second'     ,'s' ,'time',s )
    millisecond = Unit('millisecond','ms','time',ms)
    nanosecond  = Unit('nanosecond' ,'ns','time',ns)
    picosecond  = Unit('picosecond' ,'ps','time',ps)
    femtosecond = Unit('femtosecond','fs','time',fs)

    kilogram         = Unit('kilogram'        ,'kg' ,'mass',kg)
    electron_mass    = Unit('electron mass'   ,'me' ,'mass',me)
    proton_mass      = Unit('proton mass'     ,'mp' ,'mass',mp)
    atomic_mass_unit = Unit('atomic mass unit','amu','mass',amu)
    dalton           = Unit('Dalton'          ,'Da' ,'mass',Da)
        
    joule         = Unit('Joule'        ,'J'     ,'energy',J)
    electron_volt = Unit('electron Volt','eV'    ,'energy',eV)
    rydberg       = Unit('Rydberg'      ,'Ry'    ,'energy',Ry)
    hartree       = Unit('Hartree'      ,'Ha'    ,'energy',Ha)
    kJ_mole       = Unit('kJ_mole'      ,'kJ_mol','energy',kJ_mol)
    kelvin        = Unit('Kelvin'       ,'K'     ,'energy',K)
    celcius       = Unit('Celcius'      ,'degC'  ,'energy',degC,degC_shift)
    fahrenheit    = Unit('Fahrenheit'   ,'degF'  ,'energy',degF,degF_shift)

    coulomb         = Unit('Coulomb'        ,'C','charge',C)
    electron_charge = Unit('electron charge','e','charge',e)

    pascal     = Unit('Pascal'    ,'Pa'  ,'pressure',Pa)
    bar        = Unit('bar'       ,'bar' ,'pressure',bar)
    megabar    = Unit('megabar'   ,'Mbar','pressure',Mbar)
    gigapascal = Unit('gigaPascal','Gpa' ,'pressure',GPa)
    atmosphere = Unit('atmosphere','atm' ,'pressure',atm)

    newton     = Unit('Newton'    ,'N' ,'force',N)
    piconewton = Unit('picoNewton','pN','force',pN)

    W_per_mK = Unit('W/(m*K)','W_mK','thermal_cond',W_mK)


    distance_dict = dict([('A',angstrom),\
                          ('B',bohr),\
                          ('a',a),\
                          ('b',b),\
                          ('c',c),\
                          ('m',meter),\
                          ('nm',nanometer),\
                          ('pm',picometer),\
                          ('fm',femtometer),\
                          ])

    mass_dict = dict([    ('kg',kilogram),\
                          ('me',electron_mass),\
                          ('mp',proton_mass),\
                          ('amu',atomic_mass_unit),\
                          ('Da',dalton),\
                          ])
    energy_dict = dict([\
                          ('J',joule),\
                          ('eV',electron_volt),\
                          ('Ry',rydberg),\
                          ('Ha',hartree),\
                          ('kJ_mol',kJ_mole),\
                          ('K',kelvin),\
                          ('degC',celcius),\
                          ('degF',fahrenheit),\
                          ])

    charge_dict = dict([\
                          ('C',coulomb),\
                          ('e',electron_charge),\
                          ])

    pressure_dict = dict([\
                          ('Pa',pascal),\
                          ('bar',bar),\
                          ('Mbar',megabar),\
                          ('GPa',gigapascal),\
                          ('atm',atmosphere),\
                              ])

    force_dict = dict([\
                          ('N',newton),\
                          ('pN',piconewton),\
                              ])

    therm_cond_dict = dict([\
                          ('W_mK',W_per_mK),\
                          ])


    unit_dict = dict([    ('A',angstrom),\
                          ('B',bohr),\
                          ('a',a),\
                          ('b',b),\
                          ('c',c),\
                          ('m',meter),\
                          ('nm',nanometer),\
                          ('pm',picometer),\
                          ('fm',femtometer),\
                          ('kg',kilogram),\
                          ('me',electron_mass),\
                          ('mp',proton_mass),\
                          ('amu',atomic_mass_unit),\
                          ('J',joule),\
                          ('eV',electron_volt),\
                          ('Ry',rydberg),\
                          ('Ha',hartree),\
                          ('kJ_mol',kJ_mole),\
                          ('K',kelvin),\
                          ('degC',celcius),\
                          ('degF',fahrenheit),\
                          ('C',coulomb),\
                          ('e',electron_charge),\
                          ('Pa',pascal),\
                          ('bar',bar),\
                          ('Mbar',megabar),\
                          ('GPa',gigapascal),\
                          ('atm',atmosphere),\
                          ('N',newton),\
                          ('pN',piconewton),\
                          ('W_mK',W_per_mK),\
                          ])


    def __init(self):
        None
    #def __init__

    def convert(value,source_unit,target_unit):
        ui  = UnitConverter.unit_dict[source_unit]
        uo = UnitConverter.unit_dict[target_unit]

        if(ui.type != uo.type):
            print 'ERROR: in UnitConverter.convert()'
            print '   type conversion attempted between'
            print '   '+ui.type+' and '+uo.type
        else:
            value_out = (value*ui.value+ui.shift-uo.shift)/uo.value
        #end if

        return (value_out,target_unit)

    #end def convert
    convert = Callable(convert)


    def convert_scalar_to_all(units,value_orig):
        unit_type = UnitConverter.unit_dict[units].type

        value = dict()
        value['orig'] = value_orig

        for k,u in UnitConverter.unit_dict.iteritems():
            if(u.type == unit_type and u.value!=UnitConverter.unassigned):
                (value[k],utmp) = UnitConverter.convert(value['orig'],units,k)
            #end if
        #end for

        return value
    #end def convert_scalar_to_all
    convert_scalar_to_all = Callable(convert_scalar_to_all)

    
    def convert_array_to_all(units,arr,vectors=None):

        A = dict()
        A['orig'] = arr.copy()

        (n,m) = arr.shape

        if(units=='lattice'):
            if(vectors!=None):
                A['lattice'] = arr
            #end if
            def_unit = 'A'
            arr_use = dot(arr,vectors.A[def_unit])
            ui = UnitConverter.unit_dict[def_unit]
        else:
            arr_use = arr
            ui = UnitConverter.unit_dict[units]
            if(vectors!=None):
                A['lattice'] = dot(arr,linalg.inv(vectors.A[units]))
            #end if
        #end if
            
        At = zeros((n,m));
        for k,u in UnitConverter.unit_dict.iteritems():
            if(u.type == ui.type and u.value!=UnitConverter.unassigned):
                for i in range(n):
                    for j in range(m):
                        At[i,j] = (arr_use[i,j]*ui.value+ui.shift-u.shift)/u.value
                    #end for
                #end for
                A[k] = At.copy()
            #end if
        #end for

        return A
    #end def convert_array_to_all
    convert_array_to_all = Callable(convert_array_to_all)

    def submit_unit(uo):
        UnitConverter.unit_dict[uo.symbol] = uo
    #end def submit_unit
    submit_unit = Callable(submit_unit)
#end class UnitConverter


def convert(value,source_unit,target_unit):
    return UnitConverter.convert(value,source_unit,target_unit)[0]
#end def convert




try:
    from matplotlib.pyplot import figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text

    params = {'legend.fontsize':14,'figure.facecolor':'white','figure.subplot.hspace':0.,
          'axes.labelsize':16,'xtick.labelsize':14,'ytick.labelsize':14}
    rcParams.update(params)
except (ImportError,RuntimeError):
   figure,plot,xlabel,ylabel,title,show,ylim,legend,xlim,rcParams,savefig,bar,xticks,subplot,grid,setp,errorbar,loglog,semilogx,semilogy,text = unavailable('matplotlib.pyplot','figure','plot','xlabel','ylabel','title','show','ylim','legend','xlim','rcParams','savefig','bar','xticks','subplot','grid','setp','errorbar','loglog','semilogx','semilogy','text')
#end try


def sigfig_round(some_float, nsigfig):
    if nsigfig >= 1 and nsigfig == int(nsigfig):
        return round(some_float, nsigfig - int(floor(log10(abs(some_float))))-1) 
    else:
        QMCA.class_error("Number of significant figures must be a whole number\nreceived {0}".format(nsigfig))
    #end if
#end sigfig_round


def simstats(x,dim=None):
    shape = x.shape
    ndim  = len(shape)
    if dim==None:
        dim=ndim-1
    #end if
    permute = dim!=ndim-1
    reshape = ndim>2
    nblocks = shape[dim]
    if permute:
        r = range(ndim)
        r.pop(dim)
        r.append(dim)
        permutation = tuple(r)
        r = range(ndim)
        r.pop(ndim-1)
        r.insert(dim,ndim-1)
        invperm     = tuple(r)
        x=x.transpose(permutation)
        shape = tuple(array(shape)[array(permutation)])
        dim = ndim-1
    #end if
    if reshape:        
        nvars = prod(shape[0:dim])
        x=x.reshape(nvars,nblocks)
        rdim=dim
        dim=1
    else:
        nvars = shape[0]
    #end if

    mean  = x.mean(dim)
    var   = x.var(dim)

    N=nblocks

    if ndim==1:
        i=0          
        tempC=0.5
        kappa=0.0
        mtmp=mean
        if abs(var)<1e-15:
            kappa = 1.0
        else:
            ovar=1.0/var
            while (tempC>0 and i<(N-1)):
                kappa=kappa+2.0*tempC
                i=i+1
                #tempC=corr(i,x,mean,var)
                tempC = ovar/(N-i)*sum((x[0:N-i]-mtmp)*(x[i:N]-mtmp))
            #end while
            if kappa == 0.0:
                kappa = 1.0
            #end if
        #end if
        Neff=(N+0.0)/(kappa+0.0)
        if (Neff == 0.0):
            Neff = 1.0
        #end if
        error=sqrt(var/Neff)
    else:
        error = zeros(mean.shape)
        kappa = zeros(mean.shape)
        for v in xrange(nvars):
            i=0          
            tempC=0.5
            kap=0.0
            vtmp = var[v]
            mtmp = mean[v]
            if abs(vtmp)<1e-15:
                kap = 1.0
            else:
                ovar   = 1.0/vtmp
                while (tempC>0 and i<(N-1)):
                    i += 1
                    kap += 2.0*tempC
                    tempC = ovar/(N-i)*sum((x[v,0:N-i]-mtmp)*(x[v,i:N]-mtmp))
                #end while
                if kap == 0.0:
                    kap = 1.0
                #end if
            #end if
            Neff=(N+0.0)/(kap+0.0)
            if (Neff == 0.0):
                Neff = 1.0
            #end if
            kappa[v]=kap
            error[v]=sqrt(vtmp/Neff)
        #end for    
    #end if

    if reshape:
        x     =     x.reshape(shape)
        mean  =  mean.reshape(shape[0:rdim])
        var   =   var.reshape(shape[0:rdim])
        error = error.reshape(shape[0:rdim])
        kappa = kappa.reshape(shape[0:rdim])
    #end if
    if permute:
        x=x.transpose(invperm)
    #end if

    return (mean,var,error,kappa)
#end def simstats



def simplestats(x,dim=None):
    if dim==None:
        dim=len(x.shape)-1
    #end if
    osqrtN = 1.0/sqrt(1.0*x.shape[dim])
    mean   = x.mean(dim)
    var    = x.var(dim)
    error  = var*osqrtN
    return (mean,var,error,1.0)
#end def simplestats


def equilibration_length(x,tail=.5,plot=False,xlim=None,bounces=2):
    bounces = max(1,bounces)
    eqlen = 0
    nx = len(x)
    xt = x[int((1.-tail)*nx+.5):]
    nxt = len(xt)
    if nxt<10:
        return eqlen
    #end if
    #mean  = xh.mean()
    #sigma = sqrt(xh.var())
    xs = array(xt)
    xs.sort()
    mean  = xs[int(.5*(nxt-1)+.5)]
    sigma = (abs(xs[int((.5-.341)*nxt+.5)]-mean)+abs(xs[int((.5+.341)*nxt+.5)]-mean))/2
    crossings = bounces*[0,0]
    if abs(x[0]-mean)>sigma:
        s = -sign(x[0]-mean)
        ncrossings = 0
        for i in range(nx):
            dist = s*(x[i]-mean) 
            if dist>sigma and dist<5*sigma:
                crossings[ncrossings]=i
                s*=-1
                ncrossings+=1
                if ncrossings==2*bounces:
                    break
                #end if
            #end if
        #end for
        bounce = crossings[-2:]
        bounce[1] = max(bounce[1],bounce[0])
        #print len(x),crossings,crossings[1]-crossings[0]+1
        eqlen = bounce[0]+random.randint(bounce[1]-bounce[0]+1)
    #end if
    if plot:
        xlims = xlim
        del plot,xlim
        from matplotlib.pyplot import plot,figure,show,xlim
        figure()
        ix = arange(nx)
        plot(ix,x,'b.-')
        plot([0,nx],[mean,mean],'k-')
        plot([0,nx],[mean+sigma,mean+sigma],'r-')
        plot([0,nx],[mean-sigma,mean-sigma],'r-')
        plot(ix[crossings],x[crossings],'r.')
        plot(ix[bounce],x[bounce],'ro')
        plot([ix[eqlen],ix[eqlen]],[x.min(),x.max()],'g-')
        plot(ix[eqlen],x[eqlen],'go')
        if xlims!=None:
            xlim(xlims)
        #end if
        show()
    #end if
    return eqlen
#end def equilibration_length




def ci(locs,globs):
    code.interact(local=dict(globs,**locs))
#end def ci

ls = locals
gs = globals
interact = ci


class ColorWheel(DevBase):
    def __init__(self):
        colors  = 'Black Maroon DarkOrange Green DarkBlue Purple Gray Firebrick Orange MediumSeaGreen DodgerBlue MediumOrchid'.split()
        lines  = '- -- -. :'.split()
        markers = '. v s o ^ d p'.split() 
        ls = []
        for line in lines:
            for color in colors:
                ls.append((color,line))
            #end for
        #end for
        ms = []
        for i in range(len(markers)):
            ms.append((colors[i],markers[i]))
        #end for
        mls = []
        ic=-1
        for line in lines:
            for marker in markers:
                ic = (ic+1)%len(colors)
                mls.append((colors[ic],marker+line))
            #end for
        #end for
        self.line_styles   = ls
        self.marker_styles = ms
        self.marker_line_styles = mls
        self.reset()
    #end def __init__

    def next_line(self):
        self.iline = (self.iline+1)%len(self.line_styles)
        return self.line_styles[self.iline]
    #end def next_line

    def next_marker(self):
        self.imarker = (self.imarker+1)%len(self.marker_styles)
        return self.marker_styles[self.imarker]
    #end def next_marker

    def next_marker_line(self):
        self.imarker_line = (self.imarker_line+1)%len(self.marker_line_styles)
        return self.marker_line_styles[self.imarker_line]
    #end def next_marker_line

    def reset(self):
        self.iline        = -1
        self.imarker      = -1
        self.imarker_line = -1
    #end def reset

    def reset_line(self):
        self.iline        = -1
    #end def reset_line

    def reset_marker(self):
        self.imarker      = -1
    #end def reset_marker

    def reset_marker_line(self):
        self.imarker_line = -1
    #end def reset_marker_line
#end class ColorWheel
color_wheel = ColorWheel()




class QBase(DevBase):
    options = obj()

    autocorr = True

    quantity_abbreviations = obj(
        e   = 'LocalEnergy',    v  = 'Variance',       k    = 'Kinetic',
        p   = 'LocalPotential', ee = 'ElecElec',       ii   = 'IonIon',      
        l   = 'LocalECP',       n  = 'NonLocalECP',    ce   = 'CorrectedEnergy',
        mpc = 'MPC',            kc = 'KEcorr',         bw   = 'BlockWeight',
        bc  = 'BlockCPU',       ar = 'AcceptRatio',    eff  = 'Efficiency',
        te  = 'TrialEnergy',    de = 'DiffEff',        w    = 'Weight',
        nw  = 'NumOfWalkers',   sw = 'AvgSentWalkers', tt   = 'TotalTime', 
        ts  = 'TotalSamples',   le2 = 'LocalEnergy_sq'
        )

    def log(self,*texts,**kwargs):
        if len(kwargs)>0:
            n = kwargs['n']
        else:
            n=0
        #end if
        text=''
        for t in texts:
            text+=str(t)+' '
        #end for
        pad = n*'  '
        self.logfile.write(pad+text.replace('\n','\n'+pad)+'\n')
    #end def log
#end class QBase





class DatAnalyzer(QBase):
    first = 'LocalEnergy Variance Kinetic LocalPotential ElecElec LocalECP NonLocalECP IonIon'.split()
    last = 'MPC KEcorr BlockWeight BlockCPU AcceptRatio Efficiency TotalTime TotalSamples TrialEnergy DiffEff Weight NumOfWalkers LivingFraction AvgSentWalkers CorrectedEnergy'.split()

    nonenergy = set('BlockWeight BlockCPU AcceptRatio Efficiency TotalTime TotalSamples DiffEff Weight NumOfWalkers LivingFraction AvgSentWalkers'.split())
    whole_numbers = ['TotalSamples', 'NumOfWalkers']

    energy_sq = set(['Variance','LocalEnergy_sq'])

    def __init__(self,filepath,index,type='unknown',units='Ha'):
        self.info  = obj(
            index    = index,
            type     = type,
            units    = units,
            series   = None,
            filepath = None,
            filename = None,
            prefix   = None,
            nsamples = None,
            weight   = 1.0
            )
        self.data  = obj()
        self.stats = obj()
        self.load_data(filepath)
        self.analyze()
    #end def __init__


    def extract_info(self,filepath):
        series = None
        path,filename = os.path.split(filepath)
        tokens = filename.split('.')
        for token in tokens:
            if len(token)==4 and token[0]=='s' and token[1:].isdigit():
                series=int(token[1:])
                break
            #end if
        #end for
        if series is None:
            self.error('unable to read series number for file '+filepath,'QMCA')
        #end if
        prefix = filepath.split('.s'+str(series).zfill(3))[0]
        self.info.set(
            series   = series,
            filepath = filepath,
            filename = filename,
            prefix   = prefix
            )
    #end def extract_info

        
    def get_equilibration(self):
        eq = self.options.equilibration
        if eq=='auto':
            if 'LocalEnergy' in self.data:
                nbe = equilibration_length(self.data.LocalEnergy)
            else:
                nbe = 0
            #end if
        elif isinstance(eq,int):
            nbe = eq
        elif not self.info.series in eq:
            self.error('equilibration length was not provided for series {0}\nplease use the -e option to provide an equilibration length','QMCA')
        else:
            nbe = eq[self.info.series]
        #end if
        if nbe<0:
            self.error('equilibration length cannot be negative\nyou requested an equilibration length of {0} for series {1}'.format(nbe,self.info.series),'QMCA')
        #end if
        return nbe
    #end def get_equilibration


    def load_data(self,filepath):
        self.data.clear()
        self.extract_info(filepath)
        lt = loadtxt(filepath)
        if len(lt.shape)==1:
            lt.shape = (1,len(lt))
        #end if
        data = lt[:,1:].transpose()
        self.info.nsamples = data.shape[1]
        fobj = open(filepath,'r')
        variables = fobj.readline().split()[2:]
        # fobj.close()
        Econv  = convert(1.0,'Ha',self.info.units)
        Econv2 = Econv**2
        for i in range(len(variables)):
            var = variables[i]
            if var in self.nonenergy:
                self.data[var]=data[i,:]
            elif var in self.energy_sq:
                self.data[var]=Econv2*data[i,:]
            else:
                self.data[var]=Econv*data[i,:]
            #end if
        #end for
    #end def load_data


    def stat_value(self,v):
        if self.autocorr:
            (mean,var,error,kappa)=simstats(v)
        else:
            (mean,var,error,kappa)=simplestats(v)
        #end if
        sv = obj(
            mean  = mean,
            var   = var,
            error = error,
            kappa = kappa
            )
        return sv
    #end def stat_value


    def analyze(self):
        self.stats.clear()
        nbe = self.get_equilibration()
        data  = self.data
        stats = self.stats
        for varname,samples in data.iteritems():
            if nbe>=len(samples):
                self.error('equilibration length for series {0} is greater than the amount of data in the file ({1} elements)\nyou requested an equilibration length of {2}'.format(self.info.series,len(samples),nbe),'QMCA')
            #end if
            stats[varname] = self.stat_value(samples[nbe:])
        #end for
        if 'LocalEnergy_sq' in data and 'LocalEnergy' in data:
            v = data.LocalEnergy_sq - data.LocalEnergy**2
            data.Variance  = v
            stats.Variance = self.stat_value(v[nbe:])
        #end if            
        if 'LocalEnergy' in data and 'ElecElec' in data and 'MPC' in data and 'KEcorr' in data:
            v = data.LocalEnergy - data.ElecElec + data.MPC + data.KEcorr
            data.CorrectedEnergy  = v
            stats.CorrectedEnergy = self.stat_value(v[nbe:])
        #end if
        if 'LocalEnergy' in data and 'BlockWeight' in data and 'BlockCPU' in data:
            e  = stats.LocalEnergy.error
            w  = stats.BlockWeight.mean
            t  = stats.BlockCPU.mean
            wt = (data.BlockCPU*data.BlockWeight).sum()/3600
            #tt = data.BlockCPU.sum()/3600
            tt = data.BlockCPU.sum()
            ts = data.BlockWeight.sum()

            # def. of efficiency in energy.pl 
            #stats.Efficiency = obj(mean=w/t,var=0.0,error=0.0,kappa=1.0)

            # efficiency based on total work and error bar
            stats.Efficiency = obj(mean=1.0/(e**2*wt),var=0.0,error=0.0,kappa=1.0)

            data.Efficiency  = 0*data.LocalEnergy

            # total time
            stats.TotalTime = obj(mean=tt,var=0.0,error=0.0,kappa=1.0)
            data.TotalTime = 0*data.LocalEnergy

            # total samples
            stats.TotalSamples = obj(mean=ts,var=0.0,error=0.0,kappa=1.0) ###
            data.TotalSamples = 0*data.LocalEnergy

        #end if
        if 'LocalEnergy' in data and 'NumOfWalkers' in data:
            ts = data.NumOfWalkers.sum()

            stats.TotalSamples = obj(mean=ts,var=0.0,error=0.0,kappa=1.0) ###
            data.TotalSamples = 0*data.LocalEnergy

    #end def analyze


    def zero(self):
        self.info.weight = 0.0
        for d in self.data:
            d[:] = 0.0
        #end for
        for s in self.stats:
            s.set(
                mean  = 0.0,
                error = 0.0,
                kappa = 0.0,
                var   = 0.0
                )
        #end for
    #end def zero


    def accumulate(self,other,weight):
        self.info.weight += weight
        for q,v in other.data.iteritems():
            d = self.data[q]
            minlen = min(len(d),len(v))
            self.data[q] = d[0:minlen] + weight*v[0:minlen]
        #end for
        for q,v in other.stats.iteritems():
            s = self.stats[q]
            s.mean  += weight*v.mean
            s.error += weight**2*v.error**2
            s.kappa += weight*v.kappa
            s.var   += weight*v.var
        #end for
    #end def accumulate


    def normalize(self):
        w = self.info.weight
        for d in self.data:
            d[:] /= w
        #end for
        for s in self.stats:
            s.mean  /= w
            s.error  = sqrt(s.error)/w
            s.kappa /= w
            s.var   /= w*w
        #end for
    #end def normalize

    
    def write_stats(self,quantities,first_prefix=False,first_series=False,only_series=False):
        opt   = self.options
        stats = self.stats
        prefix = self.info.prefix
        series = self.info.series
        qset  = set(quantities)
        evset = set(['LocalEnergy','Variance'])
        cvset = set(['CorrectedEnergy','Variance'])
        if first_series and not only_series:
            self.log('')
        #end if

        if opt.desired_error or opt.particle_number:
            ratio = 1.0
            current_error = sigfig_round(stats.LocalEnergy.error, 1) 
            if opt.desired_error:
                ratio = stats.LocalEnergy.error / opt.desired_error
                phrase = "{:<50} {} {}".format(
                             "To change the current error of {} {} to ".format(
                                 current_error, 
                                 opt.units
                             ),
                         opt.desired_error,
                         opt.units
                         )
                self.log("\n{:<50}".format(phrase)) 
            #end if
            if opt.particle_number:
                ratio *= opt.particle_number[0]/opt.particle_number[1]
                phrase = "{:<50} {}".format(
                            "To increase current particle number of {} to".format(
                                opt.particle_number[1]
                            ),
                         opt.particle_number[0]
                         )
                if opt.desired_error:
                    self.log("{:<50}".format(phrase)) 
                else:
                    self.log("\n{:<50}".format(phrase)) 
                    phrase = "{:<50} {} {}".format(
                                "Maintaining current error bar of ",
                                    current_error, 
                                    opt.units
                             )
                    self.log("{:<50}".format(phrase)) 
                #end if
            #end if
            required_samples = int((ratio ** 2) * stats.TotalSamples.mean)
            required_samples = int(sigfig_round(required_samples, 3))
            factor = required_samples / stats.TotalSamples.mean
            factor = sigfig_round(factor, 3)
            self.log("{:<50} {:,}".format(
                        "Change sample number to",
                        required_samples
                        )
                    )
            self.log("{:<50} {}\n".format(
                        "Multiplying current sample number by a factor of",
                        factor
                        )
                    )
            exit()
        #end if

        if len(quantities)==1 and not opt.all_quantities:
            line = self.stat_line(quantities[0],prec='8.6f')
            self.log('{0}  series {1}  {2}'.format(prefix,series,line))
        elif qset==evset:
            if first_prefix and first_series:
                self.log('                            LocalEnergy               Variance           ratio')
            #end if
            line = self.EV_line(['LocalEnergy','Variance'])
            self.log('{0}  series {1}  {2}'.format(prefix,series,line))
        elif qset==cvset:
            if first_prefix and first_series:
                self.log('                        CorrectedEnergy               Variance           ratio')
            #end if
            line = self.EV_line(['CorrectedEnergy','Variance'])
            self.log('{0}  series {1}  {2}'.format(prefix,series,line))
        else:
            if opt.all_quantities:
                squants = set(stats.keys())
            else:
                squants = set(quantities)
            #end if
            first = []
            for qn in self.first:
                if qn in squants:
                    first.append(qn)
                #end if
            #end for
            last = []
            for qn in self.last:
                if qn in squants:
                    last.append(qn)
                #end if
            #end for
            mid = sorted(list(squants-set(first)-set(last)))
            order = first+mid+last
            self.log('\n{0}  series {1}'.format(prefix,series))
            for qname in order:
                if qname in stats:
                    if qname in self.whole_numbers:
                        line = self.stat_line(qname,prec='16.0f')
                    else:
                        line = self.stat_line(qname)
                    if qname=='CorrectedEnergy':
                        self.log(len(line)*'-',n=1)
                    #end if
                    self.log(line,n=1)
                #end if
            #end for
        #end if
#AcceptRatio           =          0.492020 +/-         0.000025
    #end def write_stats

    
    def stat_line(self,qname,compress=False,prec='auto'):
        q = self.stats[qname]
        if prec=='auto':
            if (q.mean !=0.0 and abs(q.error/q.mean) > 1e-8):
                d = int(max(2,1-floor(log(q.error)/log(10.))))
            else:
                d = 2
            #end if
            d = str(d)
            line = '{0:<22}=  {1:16.'+d+'f} +/- {2:16.'+d+'f}'
        else:
            line = '{0:<22}=  {1:'+prec+'} +/- {2:'+prec+'}'
        #end if
        vals = [qname,q.mean,q.error]
        opt = self.options
        if opt.show_variance:
            line += '   {'+str(len(vals))+':4.3e}'
            vals.append(q.var)
        #end if
        if opt.show_autocorr:
            line += '   {'+str(len(vals))+':>4.1f}'
            vals.append(q.kappa)
        #end if
        line = line.format(*vals)
        if compress:
            tokens = line.split()
            line=''
            for token in tokens:
                line+=token+' '
            #end for
        #end if
        return line
    #end def stat_line


    def EV_line(self,quants):
        line = ''
        for qname in quants:
            qline = self.stat_line(qname,prec='8.6f')
            line += qline.split('=')[1].strip()+'   '
        #end for
        E = self.stats[quants[0]]
        V = self.stats[quants[1]]
        line += '{0:6.4f}'.format(abs(V.mean/E.mean))
        return line
    #end def EV_line


    def plot_series(self,quantity,style,prefix):
        if not quantity in self.stats:
            self.error('quantity '+quantity+' is not present in the data for {0}'.format(self.info.filepath),'QMCA')
        #end if
        s = self.info.series
        q = self.stats[quantity]
        #errorbar(s,q.mean,q.error,fmt=style,label=prefix)
        return (s,s),(q.mean,q.mean),s,q.mean,q.error
    #end def plot_series


    def plot_trace(self,quantity,style,prefix,shift):
        if not quantity in self.data:
            self.error('quantity '+quantity+' is not present in the data for {0}'.format(self.info.filepath),'QMCA')
        #end if
        nbe = self.get_equilibration()
        q = self.data[quantity]
        middle = int(len(q)/2)
        qmean = q[middle:].mean()
        qmax = q[middle:].max()
        qmin = q[middle:].min()
        ylims = (qmean-2*(qmean-qmin),qmean+2*(qmax-qmean))
        smean,svar = self.stats[quantity].tuple('mean','var')
        sstd = sqrt(svar)
        s = arange(len(q))+shift
        xlims = min(s),max(s)
        #plot(s,q,style,label=prefix)
        plot([nbe+shift,nbe+shift],ylims,'k-.',lw=2)
        plot(xlims,[smean,smean],'r-')
        plot(xlims,[smean+sstd,smean+sstd],'r-.')
        plot(xlims,[smean-sstd,smean-sstd],'r-.')
        ylims = qmin,qmax
        return xlims,ylims,s,q
    #end def plot_trace


    def plot_histogram(self,quantity,style,prefix):
        self.not_implemented()
        #if not quantity in self.data:
        #    self.error('quantity '+quantity+' is not present in the data for {0}'.format(self.info.filepath),'QMCA')
        ##end if
        #nbe = self.get_equilibration()
        #q = self.data[quantity]
        #middle = int(len(q)/2)
        #qmean = q[middle:].mean()
        #qmax = q[middle:].max()
        #qmin = q[middle:].min()
        #ylims = (qmean-2*(qmean-qmin),qmean+2*(qmax-qmean))
        #smean,svar = self.stats[quantity].tuple('mean','var')
        #sstd = sqrt(svar)
        #s = arange(len(q))+shift
        #xlims = min(s),max(s)
        ##plot(s,q,style,label=prefix)
        #plot([nbe+shift,nbe+shift],ylims,'k-.',lw=2)
        #plot(xlims,[smean,smean],'r-')
        #plot(xlims,[smean+sstd,smean+sstd],'r-.')
        #plot(xlims,[smean-sstd,smean-sstd],'r-.')
        #ylims = qmin,qmax
        #return xlims,ylims,s,q
    #end def plot_histogram
#end class DatAnalyzer


def comma_list(s):
    if ',' in s:
        s = s.replace(',',' ')
    #end if
    s = s.strip()
    if ' ' in s:
        tokens = s.split()
        s = ''
        for t in tokens:
            s+=t+','
        #end for
        s=s[:-1]
    #end if
    return s
#end def comma_list


class QMCA(QBase):
    def __init__(self):
        self.parser    = None
        self.file_list = []
        self.file_map  = obj()
        self.data      = obj()
        self.verbose   = False

        self.read_command_line()
        self.load_data()
        self.analyze_data()
        self.exit()
    #end def __init__


    def help(self):
        self.log('\n'+self.parser.format_help().strip(),n=1)
        self.log('\nAbbreviations and full names for quantities:\n'+str(self.quantity_abbreviations),n=1)
        self.exit()
    #end def help


    def examples(self):
        examples = '''  
QMCA examples:
=============

 all quantities single file
   qmca qmc.s000.scalar.dat
   qmca qmc.s000.dmc.dat

 single quantity (local energy) single file
   qmca -q e qmc.s000.scalar.dat
   qmca -q e qmc.s000.dmc.dat

 single quantity (local energy) multiple files
  multiple series
   qmca -q e qmc.s*.scalar.dat
  multiple prefixes
   qmca -q e *.s000.scalar.dat
  multiple series and prefixes
   qmca -q e *.s*.scalar.dat

 multiple quantities, multiple files
   qmca -q ekp *.s*.scalar.dat

 special format for energy, variance, & ratio
   qmca -q ev qmc.s*.scalar.dat

 setting equilibration lengths
  equilibration length of 10 applied to all files    
   qmca -e 10 qmc.s000.scalar.dat  qmc.s001.scalar.dat
  equilibration lengths of 10 and 20 applied series 0 and 1    
   qmca -e '10 20' qmc.s000.scalar.dat  qmc.s001.scalar.dat

 twist averaging
  equal weights
   qmca -a qmc.g000.s*.scalar.dat qmc.g001.s*.scalar.dat
  user specified weights (w=2 for g=0 and w=4 for g=1)
   qmca -a -w '2 4' qmc.g000.s*.scalar.dat qmc.g001.s*.scalar.dat

 plotting vs series
   qmca -p -q e -e 10 qmc.s*.scalar.dat

 plotting vs samples (sample traces)
   qmca -t -q e -e 10 qmc.s*.scalar.dat
   
 plotting multiple twists on the same plot
   qmca -po -q e -e 10 qmc.g*.s*.scalar.dat
   qmca -to -q e -e 10 qmc.g*.s*.scalar.dat
   
 plotting twist averages
   qmca -pa -q e -e 10 qmc.g*.s*.scalar.dat
   qmca -ta -q e -e 10 qmc.g*.s*.scalar.dat

 multiple quantities and prefixes can be used for plots as well
'''
        self.log(examples)
        self.exit()
    #end def examples


    def exit(self):
        if self.verbose:
            self.log('\nQMCA finished\n')
        #end if
        exit()
    #end def exit


    def read_command_line(self):
        usage = '''usage: %prog [options] [file(s)]'''
        parser = OptionParser(usage=usage,add_help_option=False,version='%prog 0.1')
        self.parser = parser
        # available option letters: c, f, g, k, l, m, y, z
        parser.add_option('-v','--verbose',dest='verbose',
                          action='store_true',default=False,
                          help='Print detailed information (default=%default).'
                          )
        parser.add_option('-q','--quantities',dest='quantities',
                          default='all',
                          help='Quantity or list of quantities to analyze.  See names and abbreviations below (default=%default).'
                          )
        parser.add_option('-u','--units',dest='units',
                          default='Ha',
                          help='Desired energy units.  Can be Ha (Hartree), Ry (Rydberg), eV (electron volts), kJ_mol (k. joule/mole), K (Kelvin), J (Joules) (default=%default).'
                          )
        parser.add_option('-e','--equilibration',dest='equilibration',
                          default='auto',
                          help='Equilibration length in blocks (default=%default).'
                          )
        parser.add_option('-a','--average',dest='average',
                          action='store_true',default=False,
                          help='Average over files in each series (default=%default).'
                          )
        parser.add_option('-w','--weights',dest='weights',
                          default='None',
                          help='List of weights for averaging (default=%default).'
                          )
        parser.add_option('-b','--reblock',dest='reblock',
                          action='store_true',default=False,
                          help='(pending) Use reblocking to calculate statistics (default=%default).'
                          )
        parser.add_option('-p','--plot',dest='plot',
                          action='store_true',default=False,
                          help='Plot quantities vs. series (default=%default).'
                          )
        parser.add_option('-t','--trace',dest='trace',
                          action='store_true',default=False,
                          help='Plot a trace of quantities (default=%default).'
                          )
        parser.add_option('-h','--histogram',dest='histogram',
                          action='store_true',default=False,
                          help='(pending) Plot a histogram of quantities (default=%default).'
                          )
        parser.add_option('-o','--overlay',dest='overlay',
                          action='store_true',default=False,
                          help='Overlay plots (default=%default).'
                          )
        parser.add_option('--legend',dest='legend',
                          default='upper right',
                          help='Placement of legend.  None for no legend, outside for outside legend (default=%default).'
                          )
        parser.add_option('--noautocorr',dest='noautocorr',
                          action='store_true',default=False,
                          help='Do not calculate autocorrelation. Warning: error bars are no longer valid! (default=%default).'
                          )
        parser.add_option('--noac',dest='noautocorr',
                          action='store_true',default=False,
                          help='Alias for --noautocorr (default=%default).'
                          )
        parser.add_option('--sac',dest='show_autocorr',
                          action='store_true',default=False,
                          help='Show autocorrelation of sample data (default=%default).'
                          )
        parser.add_option('--sv',dest='show_variance',
                          action='store_true',default=False,
                          help='Show variance of sample data (default=%default).'
                          )
        parser.add_option('-i','--image',dest='image',
                          action='store_true',default=False,
                          help='(pending) Save image files (default=%default).'
                          )
        parser.add_option('-r','--report',dest='report',
                          action='store_true',default=False,
                          help='(pending) Write a report (default=%default).'
                          )
        parser.add_option('-s','--show_options',dest='show_options',
                          action='store_true',default=False,
                          help='Print user provided options (default=%default).'
                          )
        parser.add_option('-x','--examples',dest='examples',
                          action='store_true',default=False,
                          help='Print examples and exit (default=%default).'
                          )
        parser.add_option('--help',dest='help',
                          action='store_true',default=False,
                          help='Print help information and exit (default=%default).'
                          )
        parser.add_option('-d', '--desired_error', dest='desired_error',
                          type="float", default=None,
                          help='Show number of samples needed for desired error bar (default=%default).'
                          )
        parser.add_option('-n', '--enlarge_system',dest='particle_number',
                          type="int", nargs=2, default=None,
                          help=('Show number of samples needed to maintain '
                                'error bar on larger system: desired particle number ' 
                                'first, current particle number second (default=%default)'
                               )
                         )

        unsupported = set('histogram image reblock report examples'.split())

        units = obj(hartree='Ha',ha='Ha',ev='eV',rydberg='Ry',kelvin='K',kj_mol='kJ_mol',joules='J')

        options,files_in = parser.parse_args()
        QBase.options.transfer_from(options.__dict__)
        self.verbose = self.options.verbose
        if self.verbose:
            self.log('\nQMCA: Qmcpack Analysis Tool')
            self.log('processing options',n=1)
        #end if
        opt = self.options
        if opt.examples:
            self.examples()
        #end if
        for optname in unsupported:
            if opt[optname]!=False:
                self.log('note: option "{0}" is not yet supported'.format(optname),n=1)
                opt[optname]=False
            #end if
        #end for
        if opt.noautocorr:
            QBase.autocorr = False
        #end if
        u = opt.units.lower()
        if not u in units:
            self.error('unrecognized unit system requested: {0}\nvalid units are: {1}'.format(u, sorted(units.keys())))
        else:
            opt.units = units[u]
        #end if
        if opt.equilibration!='auto':
            opt.equilibration = comma_list(opt.equilibration)
            equil_failed = False
            try:
                opt.equilibration = int(opt.equilibration)
            except:
                try:
                    ei = eval(opt.equilibration)
                    e = obj()
                    if isinstance(ei,int):
                        e = ei
                    elif isinstance(ei,(list,tuple)):
                        ei = array(ei,dtype=int)
                        for s in range(len(ei)):
                            e[s] = ei[s]
                        #end for
                    elif isinstance(ei,dict):
                        e.transfer_from(ei)
                    else:
                        equil_failed = True
                    #end if
                    opt.equilibration = e
                except:
                    equil_failed = True
                #end try
            #end try
            if equil_failed:
                self.error('cannot process equilibration option\nvalue must be an integer, list, or dict\nyou provided: '+str(opt.equilibration))
            #end if
        #end if
        if opt.average:
            ws = opt.weights
            if ws=='None':
                opt.weights = None
            else:
                weight_failed = False
                try:
                    w = comma_list(ws)
                    w = array(eval(w),dtype=float)
                    opt.weights = w
                    weight_failed = len(w.shape)!=1
                except:
                    weight_failed = True
                #end try
                if weight_failed:
                    self.error('weights must be a list of values\nyou provided: {0}'.format(ws))
                #end if
            #end if
        #end if
        if opt.image=='None':
            opt.image=None
        #end if
        opt.all_quantities = False
        qabbr = self.quantity_abbreviations
        if ' ' in opt.quantities:
            quants = opt.quantities.split()
        else:
            quants = [opt.quantities]
        #end if
        remove = []
        for i in range(len(quants)):
            qo = quants[i]
            q = qo.lower()
            if q=='all':
                opt.all_quantities=True
                quants.extend(qabbr.values())
                remove.append(qo)
            elif not qo in qabbr.values():
                if q in qabbr:
                    quants[i] = qabbr[q]
                else:
                    match = False
                    for qa in sorted(qabbr.keys()):
                        qname = qabbr[qa]
                        if qa in q:
                            quants.append(qname)
                            q=q.replace(qa,'')
                            match=True
                        #end if
                    #end for
                    if match and q=='':
                        remove.append(qo)
                    else:
                        self.error('{0} is not a valid quantity\nplease choose quantities from the following list:\n{1}'.format(qo,str(qabbr)))
                    #end if
                #end if
            #end if
        #end for
        opt.quantities = list(set(quants)-set(remove))
        if opt.show_options or self.verbose:
            self.log('options provided:',n=1)
            self.log(str(self.options),n=1)
        #end if
        if len(files_in)==0 or opt.help:
            self.log('no files provided, please see help info below',n=1)
            self.help()
        #end if

        self.file_list.extend(files_in)
    #end def read_command_line


    def load_data(self):
        allowed_extensions = obj(scalar='scalar.dat',dmc='dmc.dat')
        opt  = self.options
        data = self.data
        ifile = 0
        seriesset = set()
        for file in self.file_list:
            if not os.path.exists(file):
                self.error('file {0} does not exist'.format(file))
            #end if
            isdat = False
            for filetype,ext in allowed_extensions.iteritems():
                if file.endswith(ext):
                    isdat = True
                    try:
                        d=DatAnalyzer(
                            filepath = file,
                            index    = ifile,
                            type     = filetype,
                            units    = opt.units
                            )
                        self.file_map[file] = d
                        prefix = d.info.prefix
                        series = d.info.series
                        seriesset.add(series)
                        if not prefix in data:
                            data[prefix]=obj()
                        #end if
                        data[prefix][series] = d
                    except:
                        self.log('problem encountered for {0}, results unavailable '.format(file))
                    #end try
                    break
                #end if
            #end for
            if not isdat:
                self.error('{0} is not a valid source of input\nplease aim qmca at files with the following extensions: {1}'.format(file,allowed_extensions.list()))
            #end if
            ifile+=1
        #end for
        if opt.average:
            for prefix,pdata in data.iteritems():
                pset = set(pdata.keys())
                missing = seriesset-pset
                if len(missing)>1:
                    self.error('files with prefix {0} are missing series provided for other prefixes\nseries provided for prefix {0}: {1}\nseries missing: {2}'.format(prefix,sorted(list(pset)),sorted(list(missing))))
                #end if
            #end for
            serieslist = sorted(list(seriesset))
            for series in serieslist:
                filetypes = set()
                for prefix,pdata in data.iteritems():
                    filetypes.add(pdata[series].info.type)
                #end for
                if len(filetypes)>1:
                    self.error('files with mixed extensions encountered for series {0}\nextensions encountered: {1}\nplease provide files with only one of the following extensions: {2}'.format(series,sorted(list(filetypes)),allowed_extensions.list()))
                #end if
            #end for
            weights = opt.weights
            prefixes = sorted(list(data.keys()))
            if weights is None:
                weights = len(prefixes)*[1.0]
            elif len(weights)!=len(prefixes):
                self.error('one weight must be provided per file prefix for averaging\nyou provided {0} weights for {1} prefixes\nweights provided: {2}\nprefixes provided: {3}'.format(len(weights),len(prefixes),weights,prefixes))
            #end if
            adata = obj()
            file_map = obj()
            file_list = []
            ns=0
            for series in serieslist:
                np = 0
                for prefix,pdata in data.iteritems():
                    w = weights[np]
                    d = pdata[series]
                    if np==0:
                        avg = d.copy()
                        di = d.info
                        avg.info.set(
                            filename = di.filename.replace(di.prefix,'avg'),
                            filepath = di.filepath.replace(di.prefix,'avg'),
                            index    = ns,
                            prefix   = 'avg'
                            )
                        avg.zero()
                    #end if
                    avg.accumulate(d,w)
                    np+=1
                #end for
                avg.normalize()
                file_list.append(avg.info.filepath)
                file_map[avg.info.filepath] = avg
                adata[series] = avg
                ns+=1
            #end for
            self.file_list = file_list
            self.file_map  = file_map
            self.data = obj(avg=adata)
        #end if
    #end def load_data


    def analyze_data(self):
        opt = self.options
        show_text  = True
        plots = opt.plot or opt.trace or opt.histogram
        if show_text:
            first_prefix = True
            for prefix in sorted(self.data.keys()):
                pdata = self.data[prefix]
                first_series = True
                only_series  = len(pdata)==1
                for series in sorted(pdata.keys()):
                    d = pdata[series]
                    d.write_stats(opt.quantities,first_prefix,first_series,only_series)
                    first_series=False
                #end for
                first_prefix = False
            #end for
        #end if
        if plots:
            color       = 'b'
            line        = '-'
            marker      = '.'
            plotstyle   = color+marker+line
            tracestyle  = color+line
            histstyle   = color+line
            iplot,itrace,ihist = range(3)
            modes = [opt.plot,opt.trace,opt.histogram]
            quantities = []
            for q in opt.quantities:
                all_present = True
                for pdata in self.data:
                    for d in pdata:
                        all_present &= q in d.data
                    #end for
                #end for
                if all_present:
                    quantities.append(q)
                #end if
            #end for
            for quantity in quantities:
                for mode in range(len(modes)):
                    if modes[mode]:
                        np = 0
                        for prefix in sorted(self.data.keys()):
                            first_prefix = np==0
                            last_prefix  = np==len(self.data)-1
                            np+=1
                            if first_prefix or not opt.overlay:
                                fig = figure()
                                ax = fig.add_subplot(111)
                                xminp =  1e99
                                xmaxp = -1e99
                                yminp =  1e99
                                ymaxp = -1e99
                                color_wheel.reset()
                            #end if
                            lcolor,lstyle = color_wheel.next_line()
                            mlcolor,mlstyle = color_wheel.next_marker_line()
                            xvals = []
                            yvals = []
                            yerrs = []
                            shifts = [0]
                            series_list = []
                            pdata = self.data[prefix]
                            only_series  = len(pdata)==1
                            shift        = 0
                            ns = 0
                            for series in sorted(pdata.keys()):
                                first_series = ns==0
                                ns+=1
                                d = pdata[series]
                                if mode==iplot:
                                    xlims,ylims,xv,yv,ye = d.plot_series(quantity,plotstyle,prefix)
                                    xvals.append(xv)
                                    yvals.append(yv)
                                    yerrs.append(ye)
                                elif mode==itrace:
                                    xlims,ylims,xv,yv = d.plot_trace(quantity,tracestyle,prefix,shift)
                                    xvals.extend(xv)
                                    yvals.extend(yv)
                                elif mode==ihist:
                                    xlims,ylims = d.plot_histogram(quantity,histstyle,prefix)
                                #end if
                                xminp = min(xminp,xlims[0])
                                xmaxp = max(xmaxp,xlims[1])
                                yminp = min(yminp,ylims[0])
                                ymaxp = max(ymaxp,ylims[1])
                                shift += d.info.nsamples
                                shifts.append(shift)
                                series_list.append(series)
                            #end for
                            if mode==iplot:
                                errorbar(xvals,yvals,yerrs,color=mlcolor,fmt=mlstyle,label=prefix)
                            elif mode==itrace:
                                plot(xvals,yvals,lstyle,color=lcolor,label=prefix)
                            #end if
                            if last_prefix or not opt.overlay:
                                dx = max(.1*abs(xminp-xmaxp),1)
                                xminp -= dx
                                xmaxp += dx
                                dy = max(.1*abs(yminp-ymaxp),1)
                                yminp -= dy
                                ymaxp += dy
                                xlim([xminp,xmaxp])
                                ylim([yminp,ymaxp])
                                if mode==iplot:
                                    xlabel('series')
                                    ylabel(quantity)
                                    title(quantity+' vs. series')
                                elif mode==itrace:
                                    for ishift in range(0,len(shifts)-1):
                                        s1 = shifts[ishift]
                                        s2 = shifts[ishift+1]
                                        smid = (s1+s2)/2
                                        if s2-s1<dx:
                                            st = s1
                                        else:
                                            st = smid
                                        #end if
                                        text(st,ymaxp-.5*dy,'s'+str(series_list[ishift]).zfill(3))
                                        if ishift!=0:
                                            plot([s1,s1],[yminp,ymaxp],'k-')
                                        #end if
                                    #end for
                                    xlabel('samples')
                                    ylabel(quantity)
                                    title('Trace of '+quantity)
                                elif mode==ihist:
                                    xlabel(quantity)
                                    ylabel('counts')
                                    title('Histogram of '+quantity)
                                #end if
                                if opt.legend.lower()!='none':
                                    if opt.legend=='upper right':
                                        ax.legend(loc='upper right')
                                    elif opt.legend=='outside':
                                        ax.legend(loc=2,bbox_to_anchor=(1.0,1.0),borderaxespad=0.)
                                    #end if
                                #end if
                            #end if
                        #end for
                    #end if
                #end for
            #end for
            show()
        #end if
    #end def analyze_data

#end class QMCA



if __name__=='__main__':
    QMCA()
#end if

   
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  